#include "CubeVertexFit.hxx"

#include <CubeReconVertex.hxx>
#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>
#include <CubeLog.hxx>
#include <CubeUnits.hxx>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

Cube::VertexFit::VertexFit() {}
Cube::VertexFit::~VertexFit() {}

namespace {
    class VertexChi2 : public ROOT::Math::IMultiGenFunction {
        Cube::Handle<Cube::ReconObjectContainer> fTracks;
    public:
        explicit VertexChi2(Cube::Handle<Cube::ReconObjectContainer> tracks)
            : fTracks(tracks) { }
        ~VertexChi2() {}
        IBaseFunctionMultiDimTempl* Clone() const {return NULL;}
        unsigned int NDim() const {
            return 3;
        }

        // The distance will be positive if the chi2 is in the "right" part of
        // the track.  The distance is for THIS end of the track, and the
        // weight is for the OTHER end of the track.
        double OtherChi2Step(double dist) const {
            double cubeSize = 10.0*unit::mm;
            // The other chi2 is off for the first cube.
            dist = dist + cubeSize;
            if (dist > 0.0) return 0.0;
            // The step size is by cubes.
            dist = dist / cubeSize;
            // The step saturates for overlaps.
            if (dist < -1.0) return 1.0;
            // There is a quadradic edge to the step.
            return 2.0*dist*dist - dist*dist*dist*dist;
        }

        // Apply a penalty when the vertex is in the middle of the track.  The
        // vertex is in the middle of the track with frDist is positive and
        // bkDist is negative.
        double MiddlePenalty(double frDist, double bkDist) const {
            double cubeSize = 10.0*unit::mm;
            frDist = (frDist - cubeSize)/cubeSize;
            bkDist = (bkDist + cubeSize)/cubeSize;
            if (frDist < 0.0) return 0.0;
            if (bkDist > 0.0) return 0.0;
            return frDist*frDist*bkDist*bkDist;
        }

        // This calculates a chi2 for the impact parameter for an end of the
        // track.  The dist, impact, and state should be for the end being
        // calculated.  The length is for the entire track.
        double StateChi2(double dist, double impact, double length,
                         Cube::Handle<Cube::TrackState> state) const {

            double cubeSize = 10.0*unit::mm;
            double var = 0.0;

            // This is the minimum direction sigma (uncertainty).  It is based
            // on the cube size and that for a straight track, the direction
            // is really determined by the two end cubes.
            double dsig = cubeSize/length;
            var += 2.0*dist*dist*dsig*dsig/12.0;

            // This is the minimum position variance based on the cube size.
            var += cubeSize/12.0;

            // Add in the variance for the position.
            var += state->GetPositionVariance().X()/3.0;
            var += state->GetPositionVariance().Y()/3.0;
            var += state->GetPositionVariance().Z()/3.0;

            // Add in the variance for the direction.
            var += dist*dist*state->GetDirectionVariance().X();
            var += dist*dist*state->GetDirectionVariance().Y();
            var += dist*dist*state->GetDirectionVariance().Z();

            return impact*impact/var;
        }

        double DoEval(const double* par) const {
            TVector3 vtx(par[0],par[1],par[2]);
            double chi2 = 0.0;
            for (Cube::ReconObjectContainer::iterator t = fTracks->begin();
                 t != fTracks->end(); ++t) {
                Cube::Handle<Cube::ReconTrack> track = (*t);

                // Find the length of the track
                double length =
                    (track->GetFront()->GetPosition().Vect()
                     - track->GetBack()->GetPosition().Vect()).Mag();

                // Parameters for the front of the track.  The front distance
                // will be negative if the vertex is in the "right" place.
                TVector3 frDiff = vtx-track->GetFront()->GetPosition().Vect();
                double frDist = frDiff * track->GetFront()->GetDirection();
                double frImpact =
                    (frDiff - frDist*track->GetFront()->GetDirection()).Mag();

                // Parameters for the back of the track.  The back distance
                // will be positive if the vertex is in the right place.
                TVector3 bkDiff = vtx-track->GetBack()->GetPosition().Vect();
                double bkDist = bkDiff * track->GetBack()->GetDirection();
                double bkImpact
                    = (bkDiff - bkDist*track->GetBack()->GetDirection()).Mag();

                // Limit the effect of the "other" end of the track.
                double frStep = OtherChi2Step(bkDist);
                double bkStep = OtherChi2Step(-frDist);

                double frChi2 = 0.0;
                if (frStep > 0.0) {
                    frChi2 = StateChi2(frDist,frImpact,length,
                                       track->GetFront());
                }

                double bkChi2 = 0.0;
                if (bkStep > 0.0) {
                    bkChi2 = StateChi2(bkDist,bkImpact,length,
                                       track->GetBack());
                }

                // Find the penalty for being in the middle of the track.
                double midPenalty = MiddlePenalty(frDist,bkDist);

                chi2 += frChi2*frStep;
                chi2 += bkChi2*bkStep;
                chi2 += midPenalty;

            }
            return chi2;
        }
    };
}

Cube::Handle<Cube::ReconVertex>
Cube::VertexFit::Apply(Cube::Handle<Cube::ReconVertex>& input) {
    CUBE_LOG(2) << "Cube::VertexFit::Apply: Begin" << std::endl;
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input->GetConstituents();
    if (!inputObjects) {
        CUBE_ERROR << "Nothing to fit." << std::endl;
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (inputObjects->empty()) {
        CUBE_ERROR << "No objects to fit." << std::endl;
        return Cube::Handle<Cube::ReconVertex>();
    }

    Cube::Handle<Cube::ReconObjectContainer>
        tracks(new Cube::ReconObjectContainer("tracks"));

    double time = 0.0;
    double timeWeight = 0.0;
    double time2 = 0.0;
    for (Cube::ReconObjectContainer::iterator o = inputObjects->begin();
         o != inputObjects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *(o);
        if (!track) continue;
        tracks->push_back(track);
        double t = track->GetPosition().T();
        double w = track->GetPositionVariance().T();
        timeWeight += w;
        time += w*t;
        time2 += w*t*t;
    }
    if (tracks->empty()) {
        CUBE_ERROR << "No tracks to fit." << std::endl;
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (tracks->size() < 2) {
        CUBE_ERROR << "Need at least 2 tracks for fit." << std::endl;
        return Cube::Handle<Cube::ReconVertex>();
    }

    time = time/timeWeight;
    time2 = time2/timeWeight;

    VertexChi2 chi2(tracks);

    std::unique_ptr<ROOT::Math::Minimizer>
        minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad"));
    if (!minimizer) {
        CUBE_ERROR << "Minimizer not created" << std::endl;
        throw std::runtime_error("Unable to create vertex fit minimizer");
    }

    minimizer->SetFunction(chi2);
    minimizer->SetVariable(0,"X",input->GetPosition().X(),1.0);
    minimizer->SetVariable(1,"Y",input->GetPosition().Y(),1.0);
    minimizer->SetVariable(2,"Z",input->GetPosition().Z(),1.0);
    minimizer->SetPrintLevel(0);

    CUBE_LOG(3) << "Starting point is "
                << " @ (" << input->GetPosition().X()
                << ", " << input->GetPosition().Y()
                << ", " << input->GetPosition().Z()
                << ", " << input->GetPosition().T()
                << ")" << std::endl;

    minimizer->Minimize();
    double minChi2 = chi2.DoEval(minimizer->X());

    Cube::Handle<Cube::VertexState> state = input->GetState();
    for (int i=0; i<3; ++i) {
        state->SetValue(state->GetPositionIndex() + i, minimizer->X()[i]);
        for (int j=0; j<3; ++j) {
            state->SetCovarianceValue(state->GetPositionIndex() + i,
                                      state->GetPositionIndex() + j,
                                      minimizer->CovMatrix(i,j));

        }
    }
    state->SetValue(state->GetTIndex(),time);
    double timeVariance = time2 - time*time;
    if (timeVariance < 0.5) timeVariance = 0.5;
    state->SetCovarianceValue(state->GetTIndex(),
                              state->GetTIndex(),
                              timeVariance);

    input->SetQuality(minChi2);

    // Emperically check for fit failures.
    if (!std::isfinite(input->GetQuality())) {
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (!std::isfinite(input->GetPosition().Mag())) {
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (!std::isfinite(input->GetPositionVariance().Mag())) {
        return Cube::Handle<Cube::ReconVertex>();
    }

    // Protect against "zero" variances.  The variance cannot be (much)
    // smaller than the cube size.
    if (input->GetPositionVariance().X()<0.0001) {
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (input->GetPositionVariance().Y()<0.0001) {
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (input->GetPositionVariance().Z()<0.0001) {
        return Cube::Handle<Cube::ReconVertex>();
    }


    CUBE_LOG(3) << "Best fit point is "
                << " @ (" << input->GetPosition().X()
                << ", " << input->GetPosition().Y()
                << ", " << input->GetPosition().Z()
                << ", " << input->GetPosition().T()
                << ")" << std::endl;

    return input;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
