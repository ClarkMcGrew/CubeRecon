#include "CubeReconVertex.hxx"
#include "CubeReconNode.hxx"

ClassImp(Cube::ReconVertex);

Cube::ReconVertex::ReconVertex()
    : ReconObject("vertex","Reconstructed SFG Vertex") {
    fState = new Cube::VertexState;
    fNodes = new Cube::ReconNodeContainerImpl<Cube::VertexState>;
}

Cube::ReconVertex::~ReconVertex() {}

TLorentzVector Cube::ReconVertex::GetPosition() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::VertexState> state = GetState();
    if (!state) throw std::runtime_error("Vertex state is missing");
    return state->GetPosition();
}

TLorentzVector Cube::ReconVertex::GetPositionVariance() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::VertexState> state = GetState();
    if (!state) throw std::runtime_error("Vertex state is missing");
    return state->GetPositionVariance();
}

bool Cube::ReconVertex::IsXVertex() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.X())) return false;
    return true;
}

bool Cube::ReconVertex::IsYVertex() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Y())) return false;
    return true;
}

bool Cube::ReconVertex::IsZVertex() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Z())) return false;
    return true;
}

int Cube::ReconVertex::GetDimensions() const{
    TLorentzVector var = GetPositionVariance();
    int dim = 0;
    if (IsXVertex()) ++dim;
    if (IsYVertex()) ++dim;
    if (IsZVertex()) ++dim;
    return dim;
}
