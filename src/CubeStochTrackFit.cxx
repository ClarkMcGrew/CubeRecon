#include "CubeStochTrackFit.hxx"
#include "SimpleSIR.hh"

#include <CubeReconNode.hxx>
#include <CubeReconCluster.hxx>
#include <CubeLog.hxx>
#include <CubeUnits.hxx>
#include <TUnitsTable.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>
#include <TRandom.h>
#include <TDecompChol.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <vector>
#include <cmath>

#define DEBUG_NUMERIC_PROBLEMS

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// This is a funky style of coding, but the SIR filter needs a lot of support
// functions that don't need to be seen by anything else.  They could be
// "declared" as methods, but that adds a bunch of visible complexity into the
// header file.  The ODD solution is to put them all here into the anonymous
// namespace.
namespace {
    typedef std::vector<float> FilterState;
    typedef Cube::Handle<Cube::ReconCluster> FilterMeasure;

    // Copied from the TCorrValue index definitions for Cube::TrackState.  These
    // definitions *MUST* *MATCH* the indices in Cube::TrackState.
    enum {
        kEDep = 0,
        kX = 1,
        kY = 2,
        kZ = 3,
        kT = 4,
        kDX = 5,
        kDY = 6,
        kDZ = 7,
        kCurvX = 8,
        kCurvY = 9,
        kCurvZ = 10,
        kWidth = 11,
        kSampleSize = 12,
    };

    /// A user likelihood for the SimpleSIR template.
    ///
    /// Implement a (mostly) Gaussian likelihood for the cluster position.
    /// This plays games to have more tail.
    struct FilterLikelihood {
        double operator () (const FilterState& state,
                            const FilterMeasure& meas) {
            double rr = 0.0;
            for (int i=0; i<3; ++i) {
                double r1 = meas->GetPosition()[i] - state[kX+i];
                rr += r1*r1/meas->GetPositionVariance()[i];
            }
            // Add more tail!
            if (rr > 1.0) rr = sqrt(rr);
            return std::exp(-rr/2.0);
        }
    };

    /// A user state propagator for the SimpleSIR template.
    ///
    /// Propagate the state to a measurement.  It takes a state and a
    /// measurement, and propagates the state to the point of closest approach
    /// to the measurement.  There are a handful of "input" fields that can be
    /// used to control the scattering as the state is propagated.
    struct FilterPropagate {
        // These need to be set before and between calls to UpdateState.
        FilterPropagate():  fEDepSigma(0),
                            fPosSigma(0), fTimeSigma(0), fVelocity(1.0),
                            fDirSigma(0), fCurvSigma(0) {}

        // The last position.  This is updated for each step and is used in
        // case the measurement is being completely missed (this happens when
        // there is a hard scatter).
        TVector3 fLastPosition;

        // A variation that is applied to the energy deposition.  This mostly
        // used to make sure that the energy deposition parameter doesn't
        // cause numeric problems.  The local energy deposition is calculated
        // based on the local average.
        double fEDepSigma;

        // A variation that is applied to the state position.  This is based
        // on multiple scattering.
        double fPosSigma;

        // A variation that is applied to the state time.  This mostly used to
        // make sure that the time parameter doesn't cause numeric problems.
        // The local time is calculated based on the local average.
        double fTimeSigma;

        // The velocity of the state (should usually be 1.0).
        double fVelocity; // Only 0.0 to 1.0 is physical.

        // A variation that is applied to the state direction.  This is based
        // on multiple scattering.
        double fDirSigma;

        // A variation that is applied to the state curvature.  This mostly
        // used to make sure that the curvature parameter doesn't cause
        // numeric problems.  The curvature is calculated based track after
        // fitting.
        double fCurvSigma;

        // A method that is used to propagate the state to the point of
        // closest approach to the measurement.  The update can work in either
        // the forward or backward direction.
        double operator () (FilterState& state,
                            const FilterMeasure& meas) {
            // Make an initial estimate of the distance to closest approach.
            double dist = 0.0;
            for (int i=0; i<3; ++i) {
                dist += state[kDX+i]*(meas->GetPosition()[i]-state[kX+i]);
            }
            // This is the changing part of the multiple scattering.  The
            // multiplicative constant needs to be set as part of fDirSigma,
            // and fPosSigma.  The std::abs() is so that the backward
            // filtering works OK.
            double multipleScatter = std::sqrt(std::abs(dist));
            // Add scattering to the direction and find the normalize.
            double norm = 0.0;
            for (int i=0; i<3; ++i) {
                if (fDirSigma > 0) {
                    double s = multipleScatter*fDirSigma;
                    state[kDX+i] += gRandom->Gaus(0.0,s);
                }
                norm += state[kDX+i]*state[kDX+i];
            }
            norm = std::sqrt(norm);
            // Update the distance to closest approach.
            dist = 0.0;
            for (int i=0; i<3; ++i) {
                dist += state[kDX+i]*(meas->GetPosition()[i]-state[kX+i])/norm;
            }
            // Update the direction normalization.  The direction is slightly
            // denormalized to prevent the variance from getting to small and
            // the correlations from getting to large.
            norm *= gRandom->Gaus(1.0,0.001);
            for (int i=0; i<3; ++i) {
                state[kDX+i] /= norm;
            }
            // Recalculate the multipleScatter (because, why not).
            multipleScatter = std::abs(dist);
            multipleScatter *= std::sqrt(multipleScatter);
            // Update the time
            state[kT] += dist/(30.0*fVelocity*unit::cm/unit::ns);
            if (fTimeSigma > 0) {
                state[kT] += gRandom->Gaus(0.0,fTimeSigma);
            }
            // Update the position.
            for (int i=0; i<3; ++i) {
                state[kX+i] += dist*state[kDX+i];
                if (fPosSigma > 0) {
                    double s = multipleScatter*fPosSigma;
                    state[kX+i] += gRandom->Gaus(0.0,s);
                }
            }
            // Update the energy deposit.
            if (fEDepSigma > 0) {
                state[kEDep] += gRandom->Gaus(0.0,fEDepSigma);
            }
            // Apply curvature
            if (fCurvSigma > 0) {
                state[kCurvX] += gRandom->Gaus(0.0,fCurvSigma);
                state[kCurvY] += gRandom->Gaus(0.0,fCurvSigma);
                state[kCurvZ] += gRandom->Gaus(0.0,fCurvSigma);
            }
#ifdef APPLY_CURVATURE
            // This update needs to be triple checked to make sure I did the
            // algebra properly...  Because of the multiple scattering,
            // I'm guessing that we aren't sensitive to the curvature in the
            // first place.
            TVector3 tempCurv(state[kCurvX],state[kCurvY],state[kCurvZ]);
            TVector3 tempDir(state[kDX],state[kDY],state[kDZ]);
            TVector3 tempCurv = tempCurv.Cross(tempDir);
            for (int i=0; i<3; ++i) {
                state[kX+i] += dist*dist*tempCurv(i);
                state[kDX+i] += dist*tempCurv(i);
            }
#endif
            // Fill the track width. These are not used.
            state[kWidth] = gRandom->Gaus();

#define ALLOW_HARD_SCATTERING
#ifdef ALLOW_HARD_SCATTERING
            // Use a "do" loop so it's easy to "break" out if needed.
            do {
                // Finally, check if there might be a hard scatter.  This
                // means that the state completely misses the measurement.
                // The new state is updated to be close to the measurement.
                TVector3 miss = meas->GetPosition().Vect()
                    - TVector3(state[kX], state[kY], state[kZ]);
                double missDist = miss.Mag();
                /// This next line isn't quite right... but not quite wrong.
                double missSigma = std::sqrt(meas->GetPositionVariance().X());
                /// Measurements that are more than this many standard
                /// deviations away from the state are considered to have been
                /// missed.
                double allowedMiss = 5.0;
                // Base the chance of treating this as a kink on the missed
                // distance.
                double kinkChance = 0.1; // sets the max probability.
                double sharpness = 6.0; // sets the transition speed.
                double totalMiss = missDist/missSigma - allowedMiss;
                kinkChance = kinkChance/(sharpness*std::exp(-totalMiss) + 1.0);
                if (kinkChance < gRandom->Uniform()) break;
                // Update the position to be near the measurement.
                miss = miss.Unit();
                double posCorr = gRandom->Gaus(missDist,missSigma);
                // Update the direction to be along a crude estimate of the
                // local direction, and make sure we don't reverse it.
                TVector3 progress = meas->GetPosition().Vect() - fLastPosition;
                progress = progress.Unit();
                double dirCorr = 0.0;
                for (int i=0; i<3; ++i) dirCorr += progress[i] * state[kDX+i];
                if (dirCorr < 0.0) progress = - progress;
                // Update the position, and direction.
                double norm = 0.0;
                for (int i=0; i<3; ++i) {
                    state[kX+i] += miss[i]*posCorr;
                    state[kDX + i] = gRandom->Gaus(progress[i],0.2);
                    norm += state[kDX+i]*state[kDX+i];
                }
                // Tweak the direction normalization to break the
                // correlations.
                norm = std::sqrt(norm)*gRandom->Gaus(1.0,0.001);
                for (int i=0; i<3; ++i) {
                    state[kDX+i] /= norm;
                }
            } while (false); // Never continue.
#endif
            return 1.0;
        }
    };

    // Declare the filter.
    typedef SimpleSIR<FilterState, FilterPropagate,
                          FilterMeasure, FilterLikelihood> FilterSIR;

    // Calculate the average and covariance for the current state vector.
    void MakeAverage(const FilterSIR::SampleVector& samples,
                           FilterState& stateAvg, TMatrixD& stateCov) {
        std::size_t dim = samples[0].second.size();
        if (stateAvg.size() != dim) {
            stateAvg.resize(dim);
            stateCov.ResizeTo(dim,dim);
        }
        for (std::size_t i=0; i<dim; ++i) {
            stateAvg[i] = 0.0;
            for (std::size_t j=0; j<dim; ++j) {
                stateCov(i,j) = 0.0;
            }
        }
        // Find the averages.
        double weight = 0.0;
        for (FilterSIR::SampleVector::const_iterator s = samples.begin();
             s != samples.end(); ++s) {
            weight += s->first;
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(s->first)) {
                CUBE_ERROR << "Invalid sample weight " << s->first << std::endl;
                throw std::runtime_error("Numeric problem");
            }
#endif
            for (std::size_t i=0; i<dim; ++i) {
#ifdef DEBUG_NUMERIC_PROBLEMS
                if (!std::isfinite(s->second[i])) {
                    CUBE_ERROR << "Invalid sample" << std::endl;
                    for (std::size_t j = 0; j<dim; ++j) {
                        CUBE_LOG(0) << "S[" << j << "] = " << s->second[j] << std::endl;
                    }
                    throw std::runtime_error("Numeric problem");
                }
#endif
                stateAvg[i] += s->first*s->second[i];
            }
        }
#ifdef DEBUG_NUMERIC_PROBLEMS
        if (weight < 0.099 || 1.001 < weight || !std::isfinite(weight)) {
            CUBE_ERROR << "Invalid weight " << weight << std::endl;
            throw std::runtime_error("Weight must be 1.0");
        }
#endif
        // Find the covariance (brute force!)
        for (FilterSIR::SampleVector::const_iterator s = samples.begin();
             s != samples.end(); ++s) {
            for (std::size_t i=0; i<dim; ++i) {
                double a = s->second[i]-stateAvg[i];
                for (std::size_t j=0; j<dim; ++j) {
                    double b = s->second[j]-stateAvg[j];
                    stateCov(i,j) += s->first*a*b;
                }
            }
        }
        // Check that the variances aren't getting very small.  This uses the
        // float epsilon, instead of double, to be extra careful.
        const double minimumVar
                = std::sqrt(std::numeric_limits<float>::epsilon());
        for (std::size_t i=0; i<dim; ++i) {
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(stateAvg[i])) {
                CUBE_ERROR << "Invalid average state " << i << std::endl;
                throw std::runtime_error("Numeric problem");
            }
#endif
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(stateCov(i,i))) {
                CUBE_ERROR << "Invalid average variance " << i << std::endl;
                throw std::runtime_error("Numeric problem");
            }
#endif
            if (stateCov(i,i) < minimumVar) {
                stateCov(i,i) = minimumVar;
            }
        }

        // Check that the correlations aren't getting very large.
        const double maxCorrelation = 0.95;
        for (std::size_t i=0; i<dim; ++i) {
            for (std::size_t j=0; j<dim; ++j) {
                if (i == j) continue;
#ifdef DEBUG_NUMERIC_PROBLEMS
                if (!std::isfinite(stateCov(i,j))) {
                    CUBE_ERROR << "Invalid average covariance " << i << " " << j << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
#endif
                double correlation = stateCov(i,j);
                correlation /= std::sqrt(stateCov(i,i));
                correlation /= std::sqrt(stateCov(j,j));
                if (correlation < maxCorrelation) continue;
                stateCov(i,j) = maxCorrelation*maxCorrelation;
                stateCov(i,j) *= std::sqrt(stateCov(i,i));
                stateCov(i,j) *= std::sqrt(stateCov(j,j));
                stateCov(j,i) = stateCov(i,j);
            }
        }
    }

    // Take a previosly filled state, and add the effect of a new state.  This
    // used to implement forward-backward smoothing.  Generically, the track
    // is first fit in one direction (filling all of the node states), and
    // then fit in the other direction.  At each node, the new average and
    // covariance are passed with the original state.  The state is modified
    // to combine information from both directional fits.
    void ForwardBackwardSmoothing(Cube::Handle<Cube::ReconState> state,
                                  int forwardMeasurements,
                                  const FilterState& stateAvg,
                                  const TMatrixD& stateCov,
                                  int backwardMeasurements) {
        // Copy the state covariance into a local matrix and invert it to turn
        // it into an error matrix.
        TMatrixD forwErr;
        forwErr.ResizeTo(stateAvg.size(),stateAvg.size());
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                forwErr(i,j) = state->GetCovarianceValue(i,j);
#ifdef DEBUG_NUMERIC_PROBLEMS
                if (!std::isfinite(forwErr(i,j))) {
                    CUBE_ERROR << "Forward covariance problem" << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
#endif
            }
        }
        // Invert forwErr to turn covariance into error matrix.
        forwErr.Invert();
#ifdef DEBUG_NUMERIC_PROBLEMS
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                if (!std::isfinite(forwErr(i,j))) {
                    CUBE_ERROR << "Forward error problem" << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
            }
        }
#endif

        // Copy the input covariance into a local matrix and invert it to turn
        // it into an error matrix.
        TMatrixD backErr;
        backErr.ResizeTo(stateAvg.size(),stateAvg.size());
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                backErr(i,j) = stateCov(i,j);
            }
        }
#ifdef DEBUG_NUMERIC_PROBLEMS
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                if (!std::isfinite(backErr(i,j))) {
                    CUBE_ERROR << "Backward covariance problem" << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
            }
        }
#endif
        // Invert backErr to turn covariance into error matrix
        backErr.Invert();
#ifdef DEBUG_NUMERIC_PROBLEMS
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                if (!std::isfinite(backErr(i,j))) {
                    CUBE_ERROR << "Backward error problem" << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
            }
        }
#endif

        // Estimate the influence of the prior on the current states.  This
        // defines a weight for the forward and backward states as the ends of
        // the track are approached.  The weights should be based on the
        // amount the prior is contributing to the covariance of the current
        // states, but that's expensive to estimate, so have the influence die
        // off as a Gaussian of the number of measurements that have been
        // added to each state.  The range of influence could be a tuned
        // parameter, but this is a good enough guess.  It's based on assuming
        // that the prior and measurements will have similar covariances.
        double influence = 1.0/5.0;

        // Find the weight of the forward measurement based on how much the
        // prior information will contribute to that state.
        double forwardWeight = influence*(forwardMeasurements-1.0);
        forwardWeight = 1.0 - std::exp(-0.5*forwardWeight*forwardWeight);

        // Find the weight of the backward measurement based on how much the
        // prior information will contribute to that state.
        double backwardWeight = influence*(backwardMeasurements-1.0); 1.0;
        backwardWeight = 1.0 - std::exp(-0.5*backwardWeight*backwardWeight);

        if (forwardWeight < 0.001 && backwardWeight < 0.001) {
            // Should never happen, but protect against a one measurement
            // track!
            forwardWeight = backwardWeight = 1.0;
        }

        CUBE_LOG(2) <<"Stochastic:: "
            "forw. " << forwardMeasurements
                    <<" ("<< forwardWeight << ")"
                    << "  back " << backwardMeasurements
                    <<" ("<< backwardWeight << ")" << std::endl;

        // Do a weighted average of the states and set the new state value.
        // This ignores the correlations so that we don't end up in the
        // paradox where the average of two points with a lot of correlations
        // in the covariance matrix will be far from both points.
        for (int i = 0; i < stateAvg.size(); ++i) {
            double wSum = 0.0;
            double w = forwErr(i,i)*forwardWeight;
            double r = w*state->GetValue(i);
            wSum += w;
            w = backErr(i,i)*backwardWeight;
            r += w*stateAvg[i];
            wSum += w;
            r /= wSum;
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(r)) {
                CUBE_ERROR << "Impossible average! "
                           << " " << i
                           << " " << r
                           << " " << wSum
                           << " " << forwardMeasurements
                           << " " << forwardWeight
                           << " " << backwardMeasurements
                           << " " << backwardWeight << std::endl;
                throw std::runtime_error("Numeric problem");
            }
#endif
            state->SetValue(i,r);
        }

        // Add the forward and backward error matrices to find the new
        // covariance.  Then invert to find the covariance.
        backErr = backErr + forwErr;
        backErr.Invert();
#ifdef DEBUG_NUMERIC_PROBLEMS
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                if (!std::isfinite(backErr(i,j))) {
                    CUBE_ERROR << "Combined covariance problem" << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
            }
        }
#endif

        // Set the state covariance.
        for (int i = 0; i < stateAvg.size(); ++i) {
            for (int j = 0; j < stateAvg.size(); ++j) {
                state->SetCovarianceValue(i,j,backErr(i,j));
            }
        }
    }

    // Estimate the curvature based on three nodes.  This does the right
    // thing, but we are so dominated by multiple scattering, that it doesn't
    // really matter.
    double FindCurvature(Cube::ReconNodeContainer::const_iterator begin,
                         Cube::ReconNodeContainer::const_iterator middle,
                         Cube::ReconNodeContainer::const_iterator end) {

        TVector3 dir(0,0,0);
        TVector3 curve(0,0,0);
        /////////////////////////////////////////////////
        // Find the local average direction and curvature
        /////////////////////////////////////////////////
        double weight = 0.0;
        Cube::Handle<Cube::ReconCluster> middleObject = (*middle)->GetObject();

        // Node before the middle node.
        Cube::Handle<Cube::ReconCluster> beginObject = (*begin)->GetObject();
        TVector3 diff = middleObject->GetPosition().Vect()
            - beginObject->GetPosition().Vect();
        double dist = diff.Mag();
        if (dist < 1.0) {
            CUBE_ERROR << "The begin and middle are " << dist << " apart" << std::endl;
            return 0.0;
        }
        double w = 1.0;
        dir += w*diff.Unit();
        curve += (w/dist)*diff.Unit();
        weight += w;

        // Node after the middle node.
        Cube::Handle<Cube::ReconCluster> endObject = (*end)->GetObject();
        diff = endObject->GetPosition().Vect()
            - middleObject->GetPosition().Vect();
        dist = diff.Mag();
        if (dist < 1.0) {
            CUBE_ERROR << "The end and middle are " << dist << " apart" << std::endl;
            return 0.0;
        }
        w = 1.0;
        dir += w*diff.Unit();
        curve -= (w/dist)*diff.Unit();
        weight += w;
        if (weight < 0.1) {
            CUBE_ERROR << "Weight is invalid" << std::endl;
            throw std::runtime_error("geometry problem");
        }
        curve = (1.0/weight)*curve;
        curve.SetX(0.0);
        return curve.Mag();
    }

    // Take all of the nodes and make a guess at the curvature by
    // statistically sampling triplets of nodes.
    void MakeCurvature(const Cube::ReconNodeContainer& nodes,
                       double& curvature, double& curvatureSigma) {
        curvature = 0.0;
        curvatureSigma = 0.0;
        // For a short track, just return a zero curvature and a large
        // uncertainty.
        if (nodes.size() < 10) {
            curvatureSigma = 1.0/(20.0*unit::cm);
            return;
        }
        double minOffset = 1.0 + nodes.size()/10.0;
        double maxOffset = 1.0 + nodes.size()/3.0;
        int maxTrial = 100.0;
        int trials = 0.0;
        for (int trial = 0; trial < maxTrial; ++trial) {
            int offset = (int) gRandom->Uniform(minOffset,maxOffset);
            int first = (int) gRandom->Uniform(0.0,nodes.size()-2*offset);
            Cube::ReconNodeContainer::const_iterator begin = nodes.begin()+first;
            Cube::ReconNodeContainer::const_iterator middle = begin + offset;
            Cube::ReconNodeContainer::const_iterator end = middle + offset;
            // Check the validity of begin, middle and end!
            if (begin == middle) {
                CUBE_ERROR << "Bad begin and middle "
                           << " Offset: " << offset
                           << " First: " << first
                           << " Middle: " << first + offset
                           << " End: " << first + 2*offset
                           << " Size: " << nodes.size() << std::endl;
                // This will (or should) cause a core dump, so force it.
                continue;
            }
            if (middle == end) {
                CUBE_ERROR << "Bad middle and end"
                           << " Offset: " << offset
                           << " First: " << first
                           << " Middle: " << first + offset
                           << " End: " << first + 2*offset
                           << " Size: " << nodes.size() << std::endl;
                // This will (or should) cause a core dump, so force it.
                continue;
            }
            if (nodes.size() < (middle - nodes.begin())) {
                CUBE_ERROR << "Bad end"
                           << " Offset: " << offset
                           << " First: " << first
                           << " Middle: " << first + offset
                           << " End: " << first + 2*offset
                           << " Size: " << nodes.size() << std::endl;
                // This will (or should) cause a core dump, so force it.
                continue;
            }
            double c = FindCurvature(begin, middle, end);
            curvature += c;
            curvatureSigma += c*c;
            trials += 1.0;
        }
        if (trials < 1.0) {
            curvature = 0.0;
            curvatureSigma = 0.005;
            return;
        }
        curvature /= maxTrial;
        curvatureSigma /= maxTrial;
        curvatureSigma = curvatureSigma - curvature*curvature;
        if (curvatureSigma > 0.0) curvatureSigma = std::sqrt(curvatureSigma);
        else curvatureSigma = 0.001 + std::abs(curvature);
    }

    // Fill the sample vector with prior states.  The vector needs to be
    // resized before this is called.  There need to be at least two
    // measurements (usually, the first few measurements should be used).  The
    // directions go from the first measurement toward the last measurement.
    void MakePrior(FilterSIR::SampleVector& samples,
                   std::vector<FilterMeasure> meas,
                   double curv, double curvSigma) {
        CUBE_LOG(2) <<"Stochastic::" <<
                       "Initial curvature " << curv << "+/-" << curvSigma << std::endl;
        double mSize = meas.size();
        // Estimate energy deposition
        double avgEDep = 0.0;
        for (int i=0; i<meas.size(); ++i) {
            avgEDep += meas[i]->GetEDeposit();
        }
        double length =
            (meas.front()->GetPosition().Vect()
             - meas.back()->GetPosition().Vect()).Mag();
        avgEDep /= length;
        for (FilterSIR::SampleVector::iterator s = samples.begin();
             s != samples.end(); ++s) {
            s->first = 1.0/samples.size();
            if (s->second.size() != kSampleSize) {
                CUBE_ERROR << "Sample size is wrong.  Problem someplace"
                           << " " << s->second.size()
                           << " vs " << kSampleSize << std::endl;
                throw std::runtime_error("Bad sample size");
            }
            int m1 = (int) gRandom->Uniform(0.0, mSize);
            int m2 = m1;
            while (m1 == m2) m2 = (int) gRandom->Uniform(0.0, mSize);
            // Fill the position.  Assume a 1cm cube size.
            for (int i=0; i<3; ++i) {
                s->second[kX+i] = meas[m1]->GetPosition()[i]
                    + gRandom->Uniform(-5.0*unit::mm,5.0*unit::mm);
            }
            s->second[kT] = gRandom->Gaus(meas[m1]->GetPosition().T(),
                                          meas[m1]->GetPositionVariance().T());
            // Fill the direction.  It will be normalized when propagated.
            double dirSign = 1.0;
            if (m2 < m1) dirSign = -1.0;
            for (int i=0; i<3; ++i) {
                s->second[kDX+i] = meas[m2]->GetPosition()[i]
                    + gRandom->Uniform(-5.0*unit::mm,5.0*unit::mm)
                    - s->second[kX+i];
                s->second[kDX+i] *= dirSign;
            }
            // Fill the energy deposition and curvature.
            s->second[kEDep] =  gRandom->Gaus(avgEDep,std::sqrt(avgEDep));
            s->second[kCurvX] =  curv;
            s->second[kCurvY] =  curv;
            s->second[kCurvZ] =  curv;
            if (curvSigma > 0.0) {
                s->second[kCurvX] += gRandom->Gaus(0.0,curvSigma);
                s->second[kCurvY] += gRandom->Gaus(0.0,curvSigma);
                s->second[kCurvZ] += gRandom->Gaus(0.0,curvSigma);
            }
            // Fill the widths (not used).  Fill with gaus to prevent a
            // singular matrix.
            s->second[kWidth] = gRandom->Gaus();
        }
    }

    // Just print things to see how they are evolving...
    void DebugState(std::string text, Cube::Handle<Cube::ReconState> state) {
        Cube::Handle<Cube::TrackState> ts = state;
        if (!ts) {
            CUBE_LOG(2) <<"Stochastic::" <<"Not a state " << text << std::endl;
            return;
        }
        CUBE_LOG(2) <<"Stochastic::" <<"Print state " << text << std::endl;
        CUBE_LOG(2) <<"Stochastic::" <<"E: "
                 << unit::AsString(ts->GetEDeposit(),
                                   std::sqrt(ts->GetEDepositVariance()),
                                   "pe") << std::endl;
        CUBE_LOG(2) <<"Stochastic::" <<"P: "
                 << "("
                 << unit::AsString(ts->GetPosition().X(),
                                   std::sqrt(ts->GetPositionVariance().X()),
                                   "length")
                 << ","
                 << unit::AsString(ts->GetPosition().Y(),
                                   std::sqrt(ts->GetPositionVariance().Y()),
                                   "length")
                 << ","
                 << unit::AsString(ts->GetPosition().Z(),
                                   std::sqrt(ts->GetPositionVariance().Z()),
                                   "length") << std::endl;
        CUBE_LOG(2) <<"Stochastic::" <<"T: "
                 << unit::AsString(ts->GetPosition().T(),
                                   std::sqrt(ts->GetPositionVariance().T()),
                                   "time") << std::endl;
        CUBE_LOG(2) <<"Stochastic::" <<"D: "
                 << "("
                 << unit::AsString(ts->GetDirection().X(),
                                   std::sqrt(ts->GetDirectionVariance().X()),
                                   "direction")
                 << ","
                 << unit::AsString(ts->GetDirection().Y(),
                                   std::sqrt(ts->GetDirectionVariance().Y()),
                                   "direction")
                 << ","
                 << unit::AsString(ts->GetDirection().Z(),
                                   std::sqrt(ts->GetDirectionVariance().Z()),
                                   "direction") << std::endl;
        CUBE_LOG(2) <<"Stochastic::" <<"D: "
                 << "("
                 << unit::AsString(ts->GetCurvature().X(),
                                   std::sqrt(ts->GetCurvatureVariance().X()),
                                   "direction")
                 << ","
                 << unit::AsString(ts->GetCurvature().Y(),
                                   std::sqrt(ts->GetCurvatureVariance().Y()),
                                   "direction")
                 << ","
                 << unit::AsString(ts->GetCurvature().Z(),
                                   std::sqrt(ts->GetCurvatureVariance().Z()),
                                   "direction") << std::endl;
    }

    // Do a resampling where all of the new states are drawn from a Gaussian
    // distribution described by the stateAvg and stateCov.  This is used to
    // make sure that the resampling doesn't collapse all of the states into a
    // single point (a typical failing of a SIR filter).
    void GaussianResample(FilterSIR::SampleVector& samples,
                          const FilterState& stateAvg,
                          const TMatrixD& stateCov) {

        TMatrixD newCov(stateCov);
        TMatrixD decomposition(stateAvg.size(),stateAvg.size());

        double minimumCov
                = std::sqrt(std::numeric_limits<float>::epsilon());
        double maxCorrelation = 0.95;
        for (int trial=0; trial<10; ++trial) {
            TDecompChol chol(newCov);
            if (chol.Decompose()) {
                decomposition = chol.GetU();
                break;
            }

            CUBE_LOG(0) << "Bad decomposition" << std::endl;
            // The Cholesky decomposition has failed.  That usually means that
            // the current estimate of the covariance has a one or more pairs
            // of variables that are too correlated, or that the variance of
            // one of the variables has become to small.  This can happen when
            // the samples gets "stuck" in one point.  A solution is to look
            // variables that have a very small variance and to increase it.
            // It also looks for pairs of variables that have a large
            // correlation and "manually" reduces the correlation.

            // Check for tiny variances.
            for (std::size_t i=0; i<stateAvg.size(); ++i) {
                if (newCov(i,i) < minimumCov) newCov(i,i) = 2.0*minimumCov;
            }

            // Check for large correlations
            for (std::size_t i=0; i<stateAvg.size(); ++i) {
                for (std::size_t j=0; j<stateAvg.size(); ++j) {
                    if (i == j) continue;
                    double correlation = newCov(i,j);
                    correlation /= std::sqrt(newCov(i,i));
                    correlation /= std::sqrt(newCov(j,j));
                    if (correlation < maxCorrelation) continue;
                    newCov(i,j) = maxCorrelation*maxCorrelation;
                    newCov(i,j) *= std::sqrt(newCov(i,i));
                    newCov(i,j) *= std::sqrt(newCov(j,j));
                    newCov(j,i) = newCov(i,j);
                }
            }
            minimumCov *= 4.0;
            maxCorrelation *= 0.9;
        }

        // Update the samples using the Cholesky decomposition of the
        // covariance.  This makes sure that we don't have any duplicated
        // samples.
        for (FilterSIR::SampleVector::iterator s = samples.begin();
             s != samples.end(); ++s) {
            // Set all the samples to have a uniform weight (normalized to 1).
            s->first = 1.0/samples.size();
            // Set the average value for the state;
            for (std::size_t i = 0; i < stateAvg.size(); ++i) {
                s->second[i] = stateAvg[i];
            }
            // Add a fluctuation around the average using the Cholesky
            // decomposition.
            for (std::size_t i = 0; i < stateAvg.size(); ++i) {
                double r = gRandom->Gaus(0.0,1.0);
                for (std::size_t j = 0; j < stateAvg.size(); ++j) {
                    s->second[j] += r*decomposition(i,j);
                }
            }
        }
    }

    // Find the multiple scattering constants for a particular mass and
    // momentum.  The mass, momentum and length are HEP units.
    double MultipleScatteringAngle(double mass, double mom, double dist) {
        if (mass < 1.0*unit::MeV) {
            // Electrons don't multiple scatter, so return a big value to take
            // into account how "flexible" the track is.
            return 0.10;
        }
        double enr = std::sqrt(mass*mass + mom*mom);
        double beta = std::sqrt(1.0-mass*mass/enr/enr);
        double sqrtRadLen = std::sqrt(41.31 * unit::cm);
        double logCorr = 1.0 - 0.038*std::log(dist/sqrtRadLen/sqrtRadLen);
        return (13.6*unit::MeV*logCorr)/(beta*mom*sqrtRadLen);
    }

    // This must be destroyed before the event is freed.
    struct NodeTimeFit : public ROOT::Math::IMultiGenFunction {
        typedef std::vector<
            std::pair<double,Cube::Handle<Cube::ReconNode>>> TravelVector;
        TravelVector fTravel;

        explicit NodeTimeFit(Cube::ReconNodeContainer& nodes) {
            fTravel.clear();
            double totalDist = 0.0;
            for (int i = 0; i < (int) nodes.size(); ++i) {
                Cube::Handle<Cube::ReconCluster> object = nodes[i]->GetObject();
                if (i < 1) {
                    fTravel.push_back(std::make_pair(totalDist,nodes[i]));
                    continue;
                }
                Cube::Handle<Cube::TrackState> state1 = nodes[i]->GetState();
                Cube::Handle<Cube::TrackState> state0 = nodes[i-1]->GetState();
                TVector3 diff
                    = state1->GetPosition().Vect()-state0->GetPosition().Vect();
                TVector3 dir = state1->GetDirection() + state0->GetDirection();
                totalDist += diff*dir.Unit();
                fTravel.push_back(std::make_pair(totalDist,nodes[i]));
            }
            // Correct so that the zero distance is close to the center of the
            // track.
            int iOffset = fTravel.size()/2;
            double offsetDist = fTravel[iOffset].first;
            for (int i = 0; i < (int) nodes.size(); ++i) {
                fTravel[i].first = fTravel[i].first - offsetDist;
            }
        }

        virtual ~NodeTimeFit() {fTravel.clear();}

        IBaseFunctionMultiDimTempl* Clone() const {return NULL;}

        unsigned int NDim() const {return 2;}

        double DoEval(const double* par) const {
            double invBeta = par[0];
            double t0 = par[1];
            // Cut off the superluminal penalty (velocities more that
            // C are OK).
            double timeGoodness = std::abs(invBeta) - 1.0;
            timeGoodness = std::max(0.0,invBeta);
            // Apply a baseline penalty for velocities less than the speed of
            // light.
            timeGoodness = 4.0*timeGoodness*timeGoodness;
            for (int i=0; i<fTravel.size(); ++i) {
                Cube::Handle<Cube::ReconCluster> cluster
                    = fTravel[i].second->GetObject();
                // The speed of light.
                const double ccc = (30.0*unit::cm)/(1.0*unit::ns);
                // Expected, measured and the variance.
                double tExp = invBeta*fTravel[i].first/ccc + t0;
                double tMea = cluster->GetPosition().T();
                double tVar = cluster->GetPositionVariance().T();
                // The timeGoodness for this measurement.  Truncated to reject
                // outliers.
                double cc = (tExp-tMea);
                cc = cc*cc/tVar;
                const double chi2Epsilon = 1E-6; // Divide by zero is bad...
                const double chi2Cut = 9.0;
                const double chi2Offset = 1.0/chi2Epsilon + 1.0/chi2Cut;
                cc = 1.0/(cc + chi2Epsilon) + 1.0/chi2Cut;
                cc = 1.0/cc - 1.0/chi2Offset;
                timeGoodness += cc;
                // Set the time.  Done this way because the position is
                // returned by value (intentionally).
                Cube::Handle<Cube::TrackState> state
                    = fTravel[i].second->GetState();
                TLorentzVector pos = state->GetPosition();
                pos.SetT(tExp);
                state->SetPosition(pos);
            }
            return timeGoodness;
        }
    };
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// End of namespace
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

Cube::StochTrackFit::StochTrackFit(int nSamples)
    : fSampleCount(nSamples) {}
Cube::StochTrackFit::~StochTrackFit() {}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Apply the actual fit here!
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Cube::Handle<Cube::ReconTrack>
Cube::StochTrackFit::Apply(Cube::Handle<Cube::ReconTrack>& input) {

    Cube::ReconNodeContainer& nodes = input->GetNodes();
    if (nodes.size() < 2) {
        CUBE_ERROR << "Not enough nodes to do a fit." << std::endl;
        return Cube::Handle<Cube::ReconTrack>();
    }

    CUBE_LOG(2) << "Stochastic::" << "Fit with " << nodes.size()
                << " nodes" << std::endl;

    int dim = nodes.front()->GetState()->GetDimensions();
    if (kSampleSize != dim) {
        std::runtime_error("Stochastic fitter state definitions are wrong");
    }

    // Define the samples that will describe the state PDF.
    FilterSIR::SampleVector samples(fSampleCount);
    for (FilterSIR::SampleVector::iterator s = samples.begin();
         s != samples.end(); ++s) {
        s->first = 1.0/dim;
        s->second.resize(dim);
    }

    CUBE_LOG(2) <<"Stochastic::" <<"Start a stochastic fit with "
              << nodes.size() << " nodes "
              << " sampled with " << fSampleCount << " samples" << std::endl;

    // Make some variables for calculating the averages and covariance from
    // the samples.
    FilterState stateAvg;
    TMatrixD stateCov;

    // Create the filter
    FilterSIR filter;
    filter.SetResampleFraction(0.0);

    // Estimate the curvature
    double priorCurvature;
    double priorCurvatureSigma;
    MakeCurvature(nodes,priorCurvature,priorCurvatureSigma);
    CUBE_LOG(2) << "Stochastic::" << "Estimated curvature prior: "
                << priorCurvature
                << " +/- " << priorCurvatureSigma << std::endl;
    // Make a prior near the front of the track.
    std::vector<FilterMeasure> priorMeasurements;
    if (nodes.size() > 5) priorMeasurements.resize(5);
    else priorMeasurements.resize(nodes.size());
    for (int i=0; i<priorMeasurements.size(); ++i) {
        priorMeasurements[i] = nodes[i]->GetObject();
    }
    MakePrior(samples,priorMeasurements,priorCurvature,priorCurvatureSigma);
    for (FilterSIR::SampleVector::iterator s = samples.begin();
         s != samples.end(); ++s) {
        filter.Propagator(s->second,priorMeasurements[0]);
    }
    MakeAverage(samples,stateAvg,stateCov);
#ifdef DEBUG_PRIOR
    for (int i=0; i<stateAvg.size(); ++i) {
        nodes.front()->GetState()->SetValue(i,stateAvg[i]);
        for (int j=0; j < stateAvg.size(); ++j) {
            nodes.front()->GetState()->SetCovarianceValue(i,j,stateCov(i,j));
        }
    }
    DebugState("Prior", nodes.front()->GetState());
#endif

    // Set the noise terms (this should be done for each step when real
    // multiple scattering is considered).  The constant values work well
    // enough for a simple fitter.
    filter.Propagator.fDirSigma
        = MultipleScatteringAngle(105*unit::MeV,500*unit::MeV,1*unit::cm);
    filter.Propagator.fPosSigma = filter.Propagator.fDirSigma/sqrt(3.0);
    filter.Propagator.fCurvSigma = 0.0001;
    filter.Propagator.fEDepSigma = 0.01;
    filter.Propagator.fTimeSigma = 0.01;
    filter.Propagator.fVelocity = 1.0;

    CUBE_LOG(2) << "Stochastic::" << "Direction sigma per sqrt(mm): "
                   << filter.Propagator.fDirSigma << std::endl;

    // Forward filter from the front to the back of the track.  This sets the
    // forward estimate of the states.
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        filter.UpdateSamples(samples,(*n)->GetObject());

        MakeAverage(samples,stateAvg,stateCov);

        // Check if we need to resample.
        double eff = filter.GetEffectiveSamples(samples.begin(), samples.end());
        if (eff < 0.5*samples.size()) {
            GaussianResample(samples,stateAvg,stateCov);
            MakeAverage(samples,stateAvg,stateCov);
        }
        filter.Propagator.fLastPosition.SetXYZ(
            stateAvg[kX],stateAvg[kY],stateAvg[kZ]);

        // Save the forward going state.
        for (int i=0; i<stateAvg.size(); ++i) {
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(stateAvg[i])) {
                CUBE_ERROR << "Invalid forward going average value"
                           << std::endl;
                throw std::runtime_error("Numeric problem");
            }
#endif
            (*n)->GetState()->SetValue(i,stateAvg[i]);
            for (int j=0; j < stateAvg.size(); ++j) {
#ifdef DEBUG_NUMERIC_PROBLEMS
                if (!std::isfinite(stateCov(i,j))) {
                    stateCov.Print();
                    CUBE_ERROR << "Invalid forward going covariance value"
                               << std::endl;
                    throw std::runtime_error("Numeric problem");
                }
#endif
                (*n)->GetState()->SetCovarianceValue(i,j,stateCov(i,j));
            }
        }
    }

    // Make a prior near the back of the track.
    for (int i = 0 ; i<priorMeasurements.size(); ++i) {
        priorMeasurements[i]
            = nodes[i+nodes.size()-priorMeasurements.size()]->GetObject();
    }
    MakePrior(samples,priorMeasurements,priorCurvature,priorCurvatureSigma);

    double energyDeposit = 0.0;
    double energyVariance = 0.0;
    double chiSquared = 0.0;

    // Backward filter from the back to the front of the track.  This then
    // applies forward-backward smoothing.
    int forwardMeas = nodes.size();
    int backwardMeas = 0;
    for (Cube::ReconNodeContainer::reverse_iterator n = nodes.rbegin();
         n != nodes.rend(); ++n) {
        filter.UpdateSamples(samples,(*n)->GetObject());

        MakeAverage(samples,stateAvg,stateCov);

        // Check if we need to resample.
        double eff = filter.GetEffectiveSamples(samples.begin(), samples.end());
        if (eff < 0.5*samples.size()) {
            GaussianResample(samples,stateAvg,stateCov);
            MakeAverage(samples,stateAvg,stateCov);
        }
        filter.Propagator.fLastPosition.SetXYZ(
            stateAvg[kX],stateAvg[kY],stateAvg[kZ]);

        // Combine states also updates the current state.  This is also
        // updating the number of measurements contributing to the forward and
        // backward states.
        ForwardBackwardSmoothing((*n)->GetState(),forwardMeas--,
                                 stateAvg,stateCov,++backwardMeas);

        // Sum up the total track energy.
        Cube::Handle<Cube::ReconCluster> cluster = (*n)->GetObject();
        Cube::Handle<Cube::ClusterState> clusterState = cluster->GetState();
        energyDeposit += clusterState->GetEDeposit();
        energyVariance += clusterState->GetEDepositVariance();

        // Calculate the goodness.  This depends on the idea that the cluster
        // covariance is always diagonal.
        Cube::Handle<Cube::TrackState> trackState = (*n)->GetState();
        TVector3 nodeDiff = trackState->GetPosition().Vect()
            - clusterState->GetPosition().Vect();
        for (int i=0; i<3; ++i) {
            chiSquared +=
                nodeDiff[i]*nodeDiff[i]/cluster->GetPositionVariance()[i];
        }
    }

    // Fill the state energy deposits.  This is implicitly assuming that the
    // nodes are uniformly spaced along the track, or put another way, each
    // node samples a uniform length of track.  Because of how SFG pattern
    // recognition works, this should USUALLY be true.
    double eSum = 0.0;
    for (int i = 0; i < (int) nodes.size(); ++i) {
        double eDep = 0.0;
        double eWght = 0.0;
        for (int j = (int) -2*fWidth-1; j < (int) 2*fWidth + 2; ++j) {
            int k = i + j;
            if (k < 1) continue;
            if (nodes.size() <= k + 1) continue;
            Cube::Handle<Cube::ReconCluster> object = nodes[k]->GetObject();
            double r = std::exp(-0.5*j*j/fWidth);
            eDep += r*object->GetEDeposit();
            eWght += r;
        }
        Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
        Cube::Handle<Cube::ReconCluster> object = nodes[i]->GetObject();
        if (eWght < 0.01) {
            state->SetEDeposit(object->GetEDeposit());
            continue;
        }
        eDep /= eWght;
        state->SetEDeposit(eDep);
        eSum += eDep;
    }
    // Force the sum of the states to equal the total deposit.
    for (int i = 0; i < (int) nodes.size(); ++i) {
        Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
        state->SetEDeposit(state->GetEDeposit()*energyDeposit/eSum);
    }

#ifndef SIMPLE_TIME_LINE_FIT
    {
        // Fit an inverse velocity (ns/mm) to the track.  The range of the
        // inverse velocity is limited so there is a minimum track velocity.
        // Superluminal velocities are allowed so that timing resolution can
        // be handled.
        std::unique_ptr<ROOT::Math::Minimizer>
            minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad"));
        if (!minimizer) {
            CUBE_ERROR << "Minimizer not created" << std::endl;
            throw std::runtime_error("Unable to create vertex fit minimizer");
        }

        NodeTimeFit timeGoodness(nodes);
        minimizer->SetFunction(timeGoodness);
        minimizer->SetVariable(0,"invBeta",0.0,1.0);
        minimizer->SetVariable(1,"T0",input->GetMedian().T(),1.0);
        minimizer->SetPrintLevel(0);

        minimizer->Minimize();
        // Make sure the nodes are filled at the best fit position.  They are
        // filled in DoEval.
        double bestTimeGoodness = timeGoodness.DoEval(minimizer->X());

        // Since the NodeTimeFit is setup as a robust estimator (not an actual
        // chi2), the estimated uncertainties are not very appropriate.  The
        // slope uncertainty is particularly poorly estimated since the slope
        // has an extra non-linear penalty term.  Based on a quick look, the
        // return Error on the fitted time offset seems to be a reasonable
        // uncertainty estimate, even if it isn't formally correct.
        double tUnc = minimizer->Errors()[1];
        for (int i = 0; i < (int) nodes.size(); ++i) {
            Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
            // Eliminate all time correlations.
            for (int j = 0; j<4; ++j)  state->SetPositionCovariance(3,j,0.0);
            // Set the time variance
            state->SetPositionCovariance(3,3,tUnc);
        }
    }
#else
    {
        // Fill the state times.  This uses a line fit of hit times and
        // (incorrectly) assumes the hits are uniformly spaced along the
        // track!!  That should usually be OK (I hope).
        double xy = 0.0;
        double xx = 0.0;
        double x = 0.0;
        double y = 0.0;
        double s = 0.0;
        for (int i = 0; i < (int) nodes.size(); ++i) {
            Cube::Handle<Cube::ReconCluster> object = nodes[i]->GetObject();
            // Never let the variance become less than this minimum.  It
            // should be a parameter.
            double v = std::max(object->GetPositionVariance().T(),
                                0.7*unit::ns);
            xx += i*i/v;
            x += i/v;
            xy += i*object->GetPosition().T()/v;
            y += object->GetPosition().T()/v;
            s += 1.0/v;
        }
        double d = x*x - s*xx;
        double m = (y*x - s*xy)/d;  // slope
        double b = (x*xy - xx*y)/d; // intercept
        double tUnc = 0.0;
        // Set the state times.
        for (int i = 0; i < (int) nodes.size(); ++i) {
            double t = m*i + b;         // Apply line
#ifdef DEBUG_NUMERIC_PROBLEMS
            if (!std::isfinite(t)) throw std::runtime_error("Bug in time");
#endif
            Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
            Cube::Handle<Cube::ReconCluster> object = nodes[i]->GetObject();
            TLorentzVector pos = state->GetPosition(); // Not by reference!!
            pos.SetT(t);
            double dt = object->GetPosition().T() - t;
            double tt = dt*dt;
            tUnc += tt*tt;
            state->SetPosition(pos);
        }
        // Set the time variances (use the same value for all states).
        if (nodes.size() < 2) {
            tUnc = 0.9*unit::ns;
        }
        else {
            tUnc /= nodes.size()-1.0;
            if (tUnc > 0) tUnc = std::sqrt(tUnc);
            tUnc = std::max(tUnc,0.9*unit::ns/std::sqrt(3.0*nodes.size()));
        }
        for (int i = 0; i < (int) nodes.size(); ++i) {
            Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
            // Eliminate all time correlations.
            for (int j = 0; j<4; ++j)  state->SetPositionCovariance(3,j,0.0);
            // Set the time variance
            state->SetPositionCovariance(3,3,tUnc);
        }
    }
#endif

    // Find the curvature and set it to be the same for all nodes.
    double curv = 0.0;
    if (nodes.size() > 5) {
        Cube::ReconNodeContainer::iterator middle = nodes.begin()
            + nodes.size()/2;
        curv = FindCurvature(nodes.begin(), middle, nodes.end()-1);
    }
    for (int i = 0; i < (int) nodes.size(); ++i) {
        Cube::Handle<Cube::TrackState> state = nodes[i]->GetState();
        state->SetCurvature(curv,curv,curv);
        state->SetCurvatureVariance(0.005,0.005,0.005);
    }

    // Fill the overall track state at the front.
    Cube::Handle<Cube::TrackState> trackState = input->GetFront();
    Cube::Handle<Cube::TrackState> nodeState = nodes.front()->GetState();
    *trackState = *nodeState;
    trackState->SetEDeposit(energyDeposit);
    trackState->SetEDepositVariance(energyVariance);

    // Fill the overall track state at the back.
    trackState = input->GetBack();
    nodeState = nodes.back()->GetState();
    *trackState = *nodeState;
    trackState->SetEDeposit(energyDeposit);
    trackState->SetEDepositVariance(energyVariance);

    DebugState("Track",input->GetFront());
    DebugState("Back",input->GetBack());

    // Setup the track information and status fields.
    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(Cube::ReconObject::kSuccess);
    input->SetStatus(Cube::ReconObject::kRan);
    input->SetAlgorithmName("StochTrackFit");
    input->SetQuality(chiSquared);
    input->SetNDOF(trackDOF);

    // Check the track direction based on timing and reverse if necessary.
    // This could be replaced by sorting the nodes according to the state
    // time.  If a sorting is done, it would handle track configurations where
    // a node is out of order.
    double dt = input->GetBack()->GetPosition().T()
        - input->GetFront()->GetPosition().T();
    if (dt < 0.0) input->ReverseTrack();

    CUBE_LOG(1) << "Stochastic::"
                << " " << nodes.size() << " nodes"
                << " with " << fSampleCount << " samples"
                << ", Chi-Squared: "  << chiSquared
                << "/" << trackDOF << " d.o.f." << std::endl;

    return input;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
