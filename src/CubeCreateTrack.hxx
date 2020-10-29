#ifndef CubeCreateTrack_hxx_seen
#define CubeCreateTrack_hxx_seen

#include "CubeClusterManagement.hxx"

#include <CubeReconTrack.hxx>
#include <CubeReconNode.hxx>
#include <CubeHandle.hxx>
#include <CubeUnits.hxx>

namespace Cube {
    /// Take iterators from a container holding TReconCluster objects in
    /// the right order for the track nodes, and construct a track.  The
    /// first argument becomes the name of the algorithm that created the
    /// track.  The track will be "fitted" by interpolating along the
    /// clusters.  The resulting track should be refit.
    template<typename clusterIterator>
    Cube::Handle<Cube::ReconTrack>
    CreateTrackFromClusters(const char* name,
                            clusterIterator begin, clusterIterator end,
                            bool verify=true);

    /// Take iterators from a container holding THandle<THit> objects, and
    /// create a track.  The hits must be in the order that they should
    /// appear in the track.  The first argument becomes the name of the
    /// algorithm that created the track.  The track will be "fitted" by
    /// interpolating along the list of hits.  The resulting track should
    /// be refit!
    template<typename hitIterator>
    Cube::Handle<Cube::ReconTrack>
    CreateTrackFromHits(const char* name,
                        hitIterator begin, hitIterator end);

    /// Combine two tracks to form a single new track with track1 at the
    /// front of the track and track2 at the back.  The resulting will
    /// need to be refitted.
    ///
    /// For efficiency within "intermediate" algorithms, the front and
    /// back states are copied from the input tracks.  This means that the
    /// new track will have the estimated front state from track1, and the
    /// estimated back state as track2.  Be careful about your assumptions
    /// when either track1 or track2 needed to be reversed, but the
    /// resulting track gets the appropriate estimate.
    Cube::Handle<Cube::ReconTrack>
    CombineTracks(const Cube::ReconTrack& track1,
                  const Cube::ReconTrack& track2);
}

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

template<typename clusterIterator>
Cube::Handle<Cube::ReconTrack>
Cube::CreateTrackFromClusters(const char* name,
                                 clusterIterator begin, clusterIterator end,
                                 bool verify) {
    if (verify) {
        for (clusterIterator i = begin; i!=end; ++i) {
            clusterIterator j = i;
            while ((++j) != end) {
                if (Cube::GetPointer(*i) != Cube::GetPointer(*j)) continue;
                CUBE_ERROR << "Invalid track: multiple copies of object"
                           << std::endl;
            }
        }
    }

    Cube::Handle<Cube::ReconTrack> track(new Cube::ReconTrack);
    track->SetAlgorithmName(name);
    track->SetStatus(Cube::ReconObject::kSuccess);
    // track->AddDetector(Cube::ReconObject::kTPC);
    track->SetName("track");

    Cube::Handle<Cube::HitSelection>
        trackHits(new Cube::HitSelection("trackHits"));
    for (clusterIterator c = begin; c != end; ++c) {
        std::copy((*c)->GetHitSelection()->begin(),
                  (*c)->GetHitSelection()->end(),
                  std::back_inserter(*trackHits));
    }
    Cube::HitSelection::iterator h = std::unique(trackHits->begin(),
                                                 trackHits->end());
    trackHits->erase(h,trackHits->end());
    track->SetHitSelection(trackHits);

    Cube::ReconNodeContainer& nodes = track->GetNodes();
    double totalCharge = 0.0;
    while (begin != end) {
        Cube::Handle<Cube::ReconCluster> cluster = *begin;
        if (!cluster) {
            CUBE_ERROR << "Invalid track: object not a cluster" << std::endl;
            return Cube::Handle<Cube::ReconTrack>();
        }
        Cube::Handle<Cube::ReconNode> node(new Cube::ReconNode);
        Cube::Handle<Cube::ReconState> state(new Cube::TrackState);
        Cube::Handle<Cube::ReconObject> object = cluster;
        node->SetState(state);
        node->SetObject(object);
        nodes.push_back(node);
        totalCharge += cluster->GetEDeposit();
        ++begin;
    }

    if (nodes.size() < 2) {
        CUBE_ERROR << "Not enough nodes for a track (" << nodes.size() << ")"
                   << std::endl;
        return Cube::Handle<Cube::ReconTrack>();
    }

    const double dirScale = 60.0*unit::mm;
    const double posScale = 30.0*unit::mm;
    const double chargeScale = 15.0*unit::mm;

    // Average the positions and directions along the track.
    for (Cube::ReconNodeContainer::iterator node = nodes.begin();
         node != nodes.end(); ++node) {
        Cube::Handle<Cube::TrackState> state = (*node)->GetState();
        Cube::Handle<Cube::ReconCluster> nodeObject = (*node)->GetObject();
        TVector3 dir(0,0,0);
        TVector3 curve(0,0,0);
        /////////////////////////////////////////////////
        // Find the local average direction and curvature
        /////////////////////////////////////////////////
        double weight = 0.0;
        Cube::ReconNodeContainer::iterator other = nodes.begin();
        while (other != node) {
            // Nodes before the current node.
            Cube::Handle<Cube::ReconCluster> otherObject
                = (*other)->GetObject();
            TVector3 diff = nodeObject->GetPosition().Vect()
                - otherObject->GetPosition().Vect();
            double dist = diff.Mag();
            double w = dist/dirScale; w = w*std::exp(-w);
            dir += w*diff.Unit();
            curve += (w/dist)*diff.Unit();
            weight += w;
            ++other;
        }
        other = node+1;
        while (other != nodes.end()) {
            // Nodes after the current node.
            Cube::Handle<Cube::ReconCluster> otherObject
                = (*other)->GetObject();
            TVector3 diff = otherObject->GetPosition().Vect()
                - nodeObject->GetPosition().Vect();
            double dist = diff.Mag();
            double w = dist/dirScale; w = w*std::exp(-w);
            dir += w*diff.Unit();
            curve -= (w/dist)*diff.Unit();
            weight += w;
            ++other;
        }
        dir = dir.Unit();
        curve = (1.0/weight)*curve;
        curve = curve - (curve*dir)*dir; // Make curve perp to dir.
        /////////////////////////////////////////////////////////////////
        // Now find the local position based on the average direction and
        // curvature.
        /////////////////////////////////////////////////////////////////
        const double objWeight = 0.1;
        TVector3 pos = nodeObject->GetPosition().Vect();
        TVector3 posVar(objWeight*pos.X()*pos.X(),
                        objWeight*pos.Y()*pos.Y(),
                        objWeight*pos.Z()*pos.Z());
        pos = objWeight*pos;
        double timeState = objWeight*nodeObject->GetPosition().T();
        double timeVar = 0.0;
        weight = objWeight;
        other = nodes.begin();
        while (other != node) {
            // Nodes before the current node.
            Cube::Handle<Cube::ReconCluster> otherObject
                = (*other)->GetObject();
            TVector3 diff = otherObject->GetPosition().Vect()
                - nodeObject->GetPosition().Vect();
            double dist = diff.Mag();
            double w = dist/posScale; w = w*std::exp(-w*w);
            TVector3 temp = otherObject->GetPosition().Vect();
            temp += dist*dir;
            pos += w*temp;
            posVar.SetX(w*temp.X()*temp.X() + posVar.X());
            posVar.SetY(w*temp.Y()*temp.Y() + posVar.Y());
            posVar.SetZ(w*temp.X()*temp.Z() + posVar.Z());
            double t = (otherObject->GetPosition().T());
            timeState += w*t;
            timeVar += w*t*t;
            weight += w;
            ++other;
        }
        other = node+1;
        while (other != nodes.end()) {
            // Nodes after the current node.
            Cube::Handle<Cube::ReconCluster> otherObject
                = (*other)->GetObject();
            TVector3 diff = otherObject->GetPosition().Vect()
                - nodeObject->GetPosition().Vect();
            double dist = diff.Mag();
            double w = dist/posScale; w = w*std::exp(-w*w);
            TVector3 temp = otherObject->GetPosition().Vect();
            temp -= dist*dir;
            pos += w*temp;
            posVar.SetX(w*temp.X()*temp.X() + posVar.X());
            posVar.SetY(w*temp.Y()*temp.Y() + posVar.Y());
            posVar.SetZ(w*temp.X()*temp.Z() + posVar.Z());
            double t = (otherObject->GetPosition().T());
            timeState += w*t;
            timeVar += w*t*t;
            weight += w;
            ++other;
        }
        pos = (1.0/weight)*pos;
        posVar = (1.0/weight)*posVar;
        posVar.SetX(posVar.X() - pos.X()*pos.X());
        posVar.SetY(posVar.Y() - pos.Y()*pos.Y());
        posVar.SetZ(posVar.Z() - pos.Z()*pos.Z());
        timeState = timeState/weight;
        timeVar = timeVar/weight - timeState*timeState;

        /////////////////////////////////////////////////////
        // Smooth the charge for this node.
        /////////////////////////////////////////////////////
        other = nodes.begin();
        weight = 0.0;
        double charge = 0.0;
        while (other != nodes.end()) {
            Cube::Handle<Cube::ReconCluster> otherObject
                = (*other)->GetObject();
            TVector3 diff = otherObject->GetPosition().Vect()
                - nodeObject->GetPosition().Vect();
            double dist = diff.Mag();
            // if (diff*dir < 0.0) dist = -dist;  // in case look at slope.
            double w = dist/chargeScale; w = std::exp(-0.5*w*w);
            charge += w*otherObject->GetEDeposit();
            weight += w;
            ++other;
        }
        charge /= weight;

        /////////////////////////////////////////////////////
        // Set the node state.
        /////////////////////////////////////////////////////
        state->SetPosition(pos.X(), pos.Y(), pos.Z(), timeState);
        state->SetPositionVariance(posVar.X(), posVar.Y(), posVar.Z(), timeVar);
        state->SetDirection(dir.X(), dir.Y(), dir.Z());
        state->SetDirectionVariance(0.001, 0.001, 0.001);
        state->SetCurvature(curve);
        state->SetCurvatureVariance(curve.X(),curve.Y(),curve.Z());
        state->SetEDeposit(charge);
        state->SetEDepositVariance(0.04*charge*charge);
    }

    // Set the total track state at the front of the track.
    Cube::Handle<Cube::TrackState> trackState = track->GetFront();
    Cube::Handle<Cube::TrackState> nodeState = nodes.front()->GetState();
    *trackState = *nodeState;

    trackState->SetEDeposit(totalCharge);
    trackState->SetEDepositVariance(0.04*totalCharge*totalCharge);

    // Set the total track state at the back of the track.
    trackState = track->GetBack();
    nodeState = nodes.back()->GetState();
    *trackState = *nodeState;

    trackState->SetEDeposit(totalCharge);
    trackState->SetEDepositVariance(0.04*totalCharge*totalCharge);

    // The states have been estimated, but to make sure the fields are
    // initialized, but mark the track as not fit anyway.
    track->ClearStatus(Cube::ReconObject::kStatusMask);

    return track;
}

template<typename hitIterator>
Cube::Handle<Cube::ReconTrack>
Cube::CreateTrackFromHits(const char* name,
                             hitIterator begin, hitIterator end) {
    // This creates one cluster per hit (EXPENSIVE!)
    Cube::Handle<Cube::ReconObjectContainer> c = CreateHitClusters(begin, end);
    return CreateTrackFromClusters(name, c->begin(), c->end(), false);
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
