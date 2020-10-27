#include "CubeBuildPairwiseVertices.hxx"

#include <CubeHandle.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconVertex.hxx>
#include <CubeLog.hxx>
#include <CubeMakeUsed.hxx>
#include <CubeUnits.hxx>
#include <CubeVertexFit.hxx>

#include <TVector3.h>
#include <TMath.h>

#include <list>
#include <set>

Cube::BuildPairwiseVertices::BuildPairwiseVertices(const char* name)
    : Cube::Algorithm(name,"Build vertices for tracks that connect") {
    fLikelihoodCut = 0.05;
    fMaxOverlap = 10.0*unit::mm;
    fMaxDistance = 50.0*unit::cm;
    fMaxApproach = 40.0*unit::mm;
    fMinTrackLength = 5.0*unit::cm;
}

Cube::BuildPairwiseVertices::~BuildPairwiseVertices() {}

Cube::Handle<Cube::AlgorithmResult>
Cube::BuildPairwiseVertices::Process(const Cube::AlgorithmResult& input,
                             const Cube::AlgorithmResult&,
                             const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "BuildPairwiseVertices:: Process "
                << GetName() <<  std::endl;
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input.GetObjectContainer();

    if (!inputObjects || !inputHits) {
        CUBE_ERROR << "No input objects or hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    // Create the output containers.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // Collect all of the tracks that can be used to build the vertices, and
    // copy all the objects to the output.
    Cube::ReconObjectContainer allTracks;
    for (Cube::ReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        Cube::Handle<Cube::ReconTrack> track = *t;
        if (!track) continue;
        allTracks.push_back(track);
    }

    CUBE_ERROR << "BuildPairwiseVertices:: Have tracks "
               << allTracks.size()<< std::endl;

    // Build all of the plausible vertices.
    std::list<Cube::Handle<Cube::ReconVertex>> allVertices;
    for (Cube::ReconObjectContainer::iterator t1 = allTracks.begin();
         t1 != allTracks.end(); ++t1) {
        if (TrackLength(*t1) < fMinTrackLength) continue;
        for (Cube::ReconObjectContainer::iterator t2 = t1 + 1;
             t2 != allTracks.end(); ++t2) {
            if (TrackLength(*t2) < fMinTrackLength) continue;
            Cube::Handle<Cube::ReconVertex> vtx = MakePairVertex((*t1),(*t2));
            if (!vtx) continue;
            allVertices.push_back(vtx);
        }
    }

    CUBE_LOG(0) << "BuildPairwiseVertices:: Have raw vertices "
                << allVertices.size() << std::endl;

    // Combine the vertices
    while (!allVertices.empty()) {
        // Find the best pair;
        std::list<Cube::Handle<Cube::ReconVertex>>::iterator vtx1
            = allVertices.end();
        std::list<Cube::Handle<Cube::ReconVertex>>::iterator vtx2 = vtx1;

        Cube::Handle<Cube::ReconVertex> bestVertex;
        double bestProb = 0.0;
        double bestDist = 1E+6;
        for (std::list<Cube::Handle<Cube::ReconVertex>>::iterator v1
                 = allVertices.begin();
             v1 != allVertices.end(); ++v1) {
            for (std::list<Cube::Handle<Cube::ReconVertex>>::iterator v2=v1;
                 v2 != allVertices.end(); ++v2) {
                if (v1 == v2) continue;
                double d = ((*v1)->GetPosition().Vect()
                            - (*v2)->GetPosition().Vect()).Mag();
                if (d > fMaxApproach) continue;
                Cube::Handle<Cube::ReconVertex> checkVertex
                    = CombineVertices(*v1,*v2);
                if (!checkVertex) continue;
                double prob = TMath::Prob(checkVertex->GetQuality(),
                                          checkVertex->GetNDOF());
                if (bestProb > prob) continue;
                vtx1 = v1;
                vtx2 = v2;
                bestDist = d;
                bestVertex = checkVertex;
                bestProb = prob;
            }
        }
        // Nothing was found (belt and suspenders)
        if (vtx1 == vtx2) break;
        if (vtx1 == allVertices.end()) break;
        if (vtx2 == allVertices.end()) break;
        if (!bestVertex) break;

        double variance = bestVertex->GetPositionVariance().Vect().Mag();
        CUBE_LOG(0) << "BuildPairwiseVertices:: Combination vertex "
                    << bestDist
                    << " T: " << bestVertex->GetConstituents()->size()
                    << " C: " << bestVertex->GetQuality()
                    << "/" << bestVertex->GetNDOF()
                    << " P: " << bestProb
                    << " S: " << std::sqrt(variance)
                    << std::endl;

        if (bestProb < fLikelihoodCut) {
            CUBE_LOG(0) << "BuildPairwiseVertices:: Vertex rejected"
                        << std::endl;
            break;
        }

        // Combine them and put it back onto the list.
        allVertices.erase(vtx1);
        allVertices.erase(vtx2);
        allVertices.push_back(bestVertex);
    }

    std::copy(allVertices.begin(), allVertices.end(),
              std::back_inserter(*finalObjects));

    CUBE_LOG(0) << "BuildPairwiseVertices:: Total vertices saved "
                << finalObjects->size() << std::endl;

    result->AddObjectContainer(finalObjects);

    // Build the hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

double Cube::BuildPairwiseVertices::TrackLength(
    Cube::Handle<Cube::ReconTrack> track) {
    return (track->GetFront()->GetPosition().Vect()
            - track->GetBack()->GetPosition().Vect()).Mag();
}

double Cube::BuildPairwiseVertices::ClosestApproach(
    const TVector3& A, const TVector3& Ad,
    const TVector3& B, const TVector3& Bd) {
    TVector3 ABcrs = Ad.Cross(Bd);
    double mag = ABcrs.Mag();
    if (mag>0) return std::abs(ABcrs*(A-B)/mag);
    return 10000;
}

double Cube::BuildPairwiseVertices::TravelDistance(
    const TVector3& A, const TVector3& Ad,
    const TVector3& B, const TVector3& Bd) {
    double impact = ClosestApproach(B, Bd, A, Ad);
    TVector3 BAcrs = Bd.Cross(Ad).Unit();
    TVector3 Ap = A + BAcrs*impact;
    double t = ((B-Ap)*Ad - ((B-Ap)*Bd)*(Bd*Ad))/(1.0-(Bd*Ad)*(Bd*Ad));
    return t;
}

TLorentzVector Cube::BuildPairwiseVertices::PairVertex(
    const TLorentzVector& A, const TVector3& Ad, double Avar,
    const TLorentzVector& B, const TVector3& Bd, double Bvar) {
    double cLight = 300.0*unit::cm/unit::ns;
    double t1 = TravelDistance(A.Vect(),Ad.Unit(), B.Vect(),Bd.Unit());
    TVector3 Ia = A.Vect() + t1*Ad.Unit();
    double t2 = TravelDistance(B.Vect(),Bd.Unit(), A.Vect(),Ad.Unit());
    TVector3 Ib = B.Vect() + t2*Bd.Unit();
    double Itime = A.T() + t1/cLight;
    Itime += B.T() + t2/cLight;
    Itime /= 2.0;
    TVector3 vtx = Ia*(1.0/Avar) + Ib*(1.0/Bvar);
    double w = 1.0/Avar + 1.0/Bvar;
    vtx = vtx * (1.0/w);
    return TLorentzVector(vtx,Itime);
}

Cube::Handle<Cube::ReconVertex>
Cube::BuildPairwiseVertices::CombineVertices(
    Cube::Handle<Cube::ReconVertex> vertex1,
    Cube::Handle<Cube::ReconVertex> vertex2) {

    Cube::Handle<Cube::ReconVertex> vtx(new Cube::ReconVertex());
    vtx->SetAlgorithmName("BuildPairwiseVertices::CombineVertices");
    vtx->SetStatus(Cube::ReconObject::kSuccess|Cube::ReconObject::kRan);
    vtx->AddDetector(Cube::ReconObject::kDST);

    double var1 = 0.0;
    var1 += vertex1->GetPositionVariance().X()/3.0;
    var1 += vertex1->GetPositionVariance().Y()/3.0;
    var1 += vertex1->GetPositionVariance().Z()/3.0;

    double var2 = 0.0;
    var2 += vertex2->GetPositionVariance().X()/3.0;
    var2 += vertex2->GetPositionVariance().Y()/3.0;
    var2 += vertex2->GetPositionVariance().Z()/3.0;

    // Find the average position for the new vertex.  This will be the
    // starting position for the fit.
    TVector3 pos = (1.0/var1)*vertex1->GetPosition().Vect()
        + (1.0/var2)*vertex2->GetPosition().Vect();
    pos = pos*(1.0/(1.0/var1 + 1.0/var2));
    double t = 0.5*(vertex1->GetPosition().T() + vertex2->GetPosition().T());
    vtx->GetState()->SetPosition(pos.X(),pos.Y(),pos.Z(),t);

    Cube::Handle<Cube::ReconObjectContainer> constituents;
    std::set<Cube::Handle<Cube::ReconTrack>> tracks;

    constituents = vertex1->GetConstituents();
    if (constituents) {
        for (Cube::ReconObjectContainer::iterator c = constituents->begin();
             c != constituents->end(); ++c) {
            tracks.insert(*c);
        }
    }

    constituents = vertex2->GetConstituents();
    if (constituents) {
        for (Cube::ReconObjectContainer::iterator c = constituents->begin();
             c != constituents->end(); ++c) {
            tracks.insert(*c);
        }
    }

    double dof = 0;
    for (std::set<Cube::Handle<Cube::ReconTrack>>::iterator t = tracks.begin();
         t != tracks.end(); ++t) {
        if (!(*t)) {
            CUBE_ERROR << "Illegal object in vertex:: Must be a track"
                       << std::endl;
            throw std::runtime_error("Illegal object in vertex");
        }
        vtx->AddConstituent(*t);
        dof += 1.0;
    }
    dof = dof - 1.0;
    vtx->SetNDOF(dof);

    if (dof < 1) return Cube::Handle<Cube::ReconVertex>();

    Cube::VertexFit vtxFit;
    vtx = vtxFit(vtx);

    if (!vtx) return Cube::Handle<Cube::ReconVertex>();
    // Protect against crazy vertices
    if (vtx->GetPositionVariance().Vect().Mag() > fMaxApproach*fMaxApproach) {
        return Cube::Handle<Cube::ReconVertex>();
    }

    return vtx;
}

Cube::Handle<Cube::ReconVertex>
Cube::BuildPairwiseVertices::MakePairVertex(
    Cube::Handle<Cube::ReconTrack> track1,
    Cube::Handle<Cube::ReconTrack> track2) {
    double bestApproach = 100000;
    double bestDist1 = 0.0;
    double bestDist2 = 0.0;
    TLorentzVector bestPos1;
    TLorentzVector bestPos2;
    TVector3 bestDir1;
    TVector3 bestDir2;

    // front front
    double t1 = TravelDistance(track1->GetFront()->GetPosition().Vect(),
                               track1->GetFront()->GetDirection(),
                               track2->GetFront()->GetPosition().Vect(),
                               track2->GetFront()->GetDirection());
    double t2 = TravelDistance(track2->GetFront()->GetPosition().Vect(),
                               track2->GetFront()->GetDirection(),
                               track1->GetFront()->GetPosition().Vect(),
                               track1->GetFront()->GetDirection());
    double b = ClosestApproach(track1->GetFront()->GetPosition().Vect(),
                               track1->GetFront()->GetDirection(),
                               track2->GetFront()->GetPosition().Vect(),
                               track2->GetFront()->GetDirection());
    if (t1 < fMaxOverlap && t2 < fMaxOverlap) {
        if (b < bestApproach) {
            bestApproach = b;
            bestDist1 = t1;
            bestDist2 = t2;
            bestPos1 = track1->GetFront()->GetPosition();
            bestDir1 = track1->GetFront()->GetDirection();
            bestPos2 = track2->GetFront()->GetPosition();
            bestDir2 = track2->GetFront()->GetDirection();
        }
    }

    // front back
    t1 = TravelDistance(track1->GetFront()->GetPosition().Vect(),
                        track1->GetFront()->GetDirection(),
                        track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection());
    t2 = TravelDistance(track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection(),
                        track1->GetFront()->GetPosition().Vect(),
                        track1->GetFront()->GetDirection());
    b = ClosestApproach(track1->GetFront()->GetPosition().Vect(),
                        track1->GetFront()->GetDirection(),
                        track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection());
    if (t1 < fMaxOverlap && t2 > - fMaxOverlap) {
        if (b < bestApproach) {
            bestApproach = b;
            bestDist1 = t1;
            bestDist2 = t2;
            bestPos1 = track1->GetFront()->GetPosition();
            bestDir1 = track1->GetFront()->GetDirection();
            bestPos2 = track2->GetBack()->GetPosition();
            bestDir2 = track2->GetBack()->GetDirection();
        }
    }

    // back front
    t1 = TravelDistance(track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection(),
                        track2->GetFront()->GetPosition().Vect(),
                        track2->GetFront()->GetDirection());
    t2 = TravelDistance(track2->GetFront()->GetPosition().Vect(),
                        track2->GetFront()->GetDirection(),
                        track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection());
    b = ClosestApproach(track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection(),
                        track2->GetFront()->GetPosition().Vect(),
                        track2->GetFront()->GetDirection());
    if (t1 > -fMaxOverlap && t2 < fMaxOverlap) {
        if (b < bestApproach) {
            bestApproach = b;
            bestDist1 = t1;
            bestDist2 = t2;
            bestPos1 = track1->GetBack()->GetPosition();
            bestDir1 = track1->GetBack()->GetDirection();
            bestPos2 = track2->GetFront()->GetPosition();
            bestDir2 = track2->GetFront()->GetDirection();
        }
    }

    // back back
    t1 = TravelDistance(track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection(),
                        track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection());
    t2 = TravelDistance(track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection(),
                        track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection());
    b = ClosestApproach(track1->GetBack()->GetPosition().Vect(),
                        track1->GetBack()->GetDirection(),
                        track2->GetBack()->GetPosition().Vect(),
                        track2->GetBack()->GetDirection());
    if (t1 > -fMaxOverlap && t2 > -fMaxOverlap) {
        if (b < bestApproach) {
            bestApproach = b;
            bestDist1 = t1;
            bestDist2 = t2;
            bestPos1 = track1->GetBack()->GetPosition();
            bestDir1 = track1->GetBack()->GetDirection();
            bestPos2 = track2->GetBack()->GetPosition();
            bestDir2 = track2->GetBack()->GetDirection();
        }
    }

    // Check if this is a decent vertex pair.
    if (bestApproach > fMaxApproach) return Cube::Handle<Cube::ReconVertex>();
    if (std::abs(bestDist1) > fMaxDistance) {
        return Cube::Handle<Cube::ReconVertex>();
    }
    if (std::abs(bestDist2) > fMaxDistance) {
        return Cube::Handle<Cube::ReconVertex>();
    }

    double t1Len = TrackLength(track1);
    // Add Uncertainty for the discreetness of the cubes.
    double t1Var = 10*unit::mm/t1Len;  //Temp: Min dir sigma (1 cube/length).
    t1Var = bestDist1*bestDist1*t1Var*t1Var + 100.0*unit::mm*unit::mm/12.0;
    TVector3 t1DirVar = track1->GetState()->GetDirectionVariance();
    t1Var += track1->GetPositionVariance().X()/3.0;
    t1Var += track1->GetPositionVariance().Y()/3.0;
    t1Var += track1->GetPositionVariance().Z()/3.0;
    t1Var += bestDist1*bestDist1*t1DirVar.X();
    t1Var += bestDist1*bestDist1*t1DirVar.Y();
    t1Var += bestDist1*bestDist1*t1DirVar.Z();

    double t2Len = TrackLength(track2);
    // Add Uncertainty for the discreetness of the cubes.
    double t2Var = 10*unit::mm/t2Len;  //Temp: Min dir sigma (1 cube/length).
    t2Var = bestDist2*bestDist2*t2Var*t2Var + 100.0*unit::mm*unit::mm/12.0;
    TVector3 t2DirVar = track2->GetState()->GetDirectionVariance();
    t2Var += track2->GetPositionVariance().X()/3.0;
    t2Var += track2->GetPositionVariance().Y()/3.0;
    t2Var += track2->GetPositionVariance().Z()/3.0;
    t2Var += bestDist2*bestDist2*t2DirVar.X();
    t2Var += bestDist2*bestDist2*t2DirVar.Y();
    t2Var += bestDist2*bestDist2*t2DirVar.Z();

    TLorentzVector bestVertex = PairVertex(bestPos1,bestDir1,t1Var,
                                           bestPos2,bestDir2,t2Var);

    Cube::Handle<Cube::ReconVertex> vtx(new Cube::ReconVertex());
    vtx->SetAlgorithmName("BuildPairwiseVertices::MakePairVertex");
    vtx->SetStatus(Cube::ReconObject::kSuccess|Cube::ReconObject::kRan);
    vtx->AddDetector(Cube::ReconObject::kDST);
    vtx->SetNDOF(1.0);
    vtx->SetQuality(1.0);

    vtx->GetState()->SetPosition(bestVertex);
    vtx->GetState()->SetPositionVariance(10.0,10.0,10.0,10.0);

    vtx->AddConstituent(track1);
    vtx->AddConstituent(track2);

    Cube::VertexFit vtxFit;
    vtx = vtxFit(vtx);

    if (!vtx) return Cube::Handle<Cube::ReconVertex>();
    if (vtx->GetPositionVariance().Vect().Mag() > fMaxApproach*fMaxApproach) {
        return Cube::Handle<Cube::ReconVertex>();
    }

    return vtx;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
