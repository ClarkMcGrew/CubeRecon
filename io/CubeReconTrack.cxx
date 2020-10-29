#include "CubeReconTrack.hxx"
#include "CubeCorrValues.hxx"
#include "CubeReconNode.hxx"

#include <TROOT.h>

ClassImp(Cube::ReconTrack);

Cube::ReconTrack::ReconTrack()
    : ReconObject("track","Reconstructed SFG Track") {
    fState = new Cube::TrackState;
    fBackState = new Cube::TrackState;
    fNodes = new Cube::ReconNodeContainerImpl<Cube::TrackState>;
}

Cube::ReconTrack::ReconTrack(const Cube::ReconTrack& track)
    : Cube::ReconObject(track) {
    fNodes = new Cube::ReconNodeContainerImpl<Cube::TrackState>;

    // Copy the nodes
    // Create new nodes with Cube::TrackState's
    Cube::ReconNodeContainer::const_iterator in;
    for (in=track.GetNodes().begin(); in!=track.GetNodes().end(); ++in){
        Cube::Handle<Cube::ReconNode> node(new Cube::ReconNode);
        Cube::Handle<Cube::ReconObject> object = (*in)->GetObject();
        node->SetObject(object);
        Cube::Handle<Cube::TrackState> tstate = (*in)->GetState();
        if (tstate) {
            Cube::Handle<Cube::ReconState> pstate(
                new Cube::TrackState(*tstate));
            node->SetState(pstate);
        }
        node->SetQuality((*in)->GetQuality());

        fNodes->push_back(node);
    }

    if (track.GetState()) {
        Cube::Handle<Cube::TrackState> state = track.GetState();
        fState = new Cube::TrackState(*state);
    }
    else {
        fState = new Cube::TrackState;
    }

    if (track.GetBack()) {
        Cube::Handle<Cube::TrackState> state = track.GetBack();
        fBackState = new Cube::TrackState(*state);
    }
    else {
        fBackState = new Cube::TrackState;
    }

}

Cube::ReconTrack::~ReconTrack() {}

double Cube::ReconTrack::GetEDeposit() const {
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetEDeposit();
}

TLorentzVector Cube::ReconTrack::GetPosition() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetPosition();
}

TLorentzVector Cube::ReconTrack::GetPositionVariance() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetPositionVariance();
}

bool Cube::ReconTrack::IsXTrack() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.X())) return false;
    return true;
}

bool Cube::ReconTrack::IsYTrack() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Y())) return false;
    return true;
}

bool Cube::ReconTrack::IsZTrack() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Z())) return false;
    return true;
}

int Cube::ReconTrack::GetDimensions() const{
    TLorentzVector var = GetPositionVariance();
    int dim = 0;
    if (IsXTrack()) ++dim;
    if (IsYTrack()) ++dim;
    if (IsZTrack()) ++dim;
    return dim;
}

TVector3 Cube::ReconTrack::GetDirection() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetDirection();
}

TVector3 Cube::ReconTrack::GetCurvature() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetCurvature();
}

double Cube::ReconTrack::GetWidth() const {
    Cube::Handle<Cube::TrackState> state = GetState();
    if (!state) throw std::runtime_error("Track state messing");
    return state->GetWidth();
}

void Cube::ReconTrack::ReverseTrack() {
    // Reverse the order of the nodes.
    std::reverse(GetNodes().begin(), GetNodes().end());

    // Reverse the state directions.
    for (Cube::ReconNodeContainer::iterator n = GetNodes().begin();
         n != GetNodes().end(); ++n) {
        Cube::Handle<Cube::TrackState> state = (*n)->GetState();
        state->SetDirection(-state->GetDirection());
    }

    // Swap the front and back states.
    Cube::TrackState tempState(*(GetFront()));
    *GetFront() = *GetBack();
    *GetBack() = tempState;

    // Reverse the direction at the front and back.
    GetFront()->SetDirection(-GetFront()->GetDirection());
    GetBack()->SetDirection(-GetBack()->GetDirection());

    // Swap the energy deposit at the front and back.
    double tempDeposit = GetFront()->GetEDeposit();
    GetFront()->SetEDeposit(GetBack()->GetEDeposit());
    GetBack()->SetEDeposit(tempDeposit);

}

void Cube::ReconTrack::ls(Option_t *opt) const {
    ls_base(opt);

    TROOT::IncreaseDirLevel();
    std::string option(opt);
    if (fState) {
        TROOT::IncreaseDirLevel();
        fState->ls(opt);
        TROOT::DecreaseDirLevel();
    }
    if (fBackState) {
        TROOT::IncreaseDirLevel();
        fBackState->ls(opt);
        TROOT::DecreaseDirLevel();
    }
    if (fNodes && (option.find("dump") != std::string::npos
                   || option.find("recon") != std::string::npos)) {
        TROOT::IncreaseDirLevel();
        fNodes->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
