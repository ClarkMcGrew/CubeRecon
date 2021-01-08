#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <TROOT.h>

#include "CubeReconCluster.hxx"
#include "CubeReconNode.hxx"

ClassImp(Cube::ReconCluster);

Cube::ReconCluster::ReconCluster()
    : ReconObject("cluster","Reconstructed SFG Cluster"),
      fMoments(3), fTemporariesInitialized(false) {
    fState = new Cube::ClusterState;
    fNodes = new Cube::ReconNodeContainerImpl<Cube::ClusterState>;
}

Cube::ReconCluster::ReconCluster(const Cube::ReconCluster& cluster)
    : Cube::ReconObject(cluster), fMoments(3), fTemporariesInitialized(false) {
    fNodes = new Cube::ReconNodeContainerImpl<Cube::ClusterState>;

    // Copy the nodes.  Create new nodes with Cube::ClusterState's
    Cube::ReconNodeContainer::const_iterator in;
    for (in=cluster.GetNodes().begin(); in!=cluster.GetNodes().end();in++){
        Cube::Handle<Cube::ReconNode> node(new Cube::ReconNode);
        Cube::Handle<Cube::ReconObject> object = (*in)->GetObject();
        node->SetObject(object);
        Cube::Handle<Cube::ClusterState> tstate = (*in)->GetState();
        if (tstate){
            Cube::Handle<Cube::ReconState> pstate(
                new Cube::ClusterState(*tstate));
            node->SetState(pstate);
        }
        node->SetQuality((*in)->GetQuality());

        fNodes->push_back(node);
    }


    if (cluster.GetState()){
        Cube::Handle<Cube::ClusterState> state = cluster.GetState();
        fState = new Cube::ClusterState(*state);
    }
    else {
        fState = new Cube::ClusterState;
    }
}

Cube::ReconCluster::~ReconCluster() {}

double Cube::ReconCluster::GetEDeposit() const {
    // I'm being a bit pedantic and casting to the base mix-in class.  This
    // could just as well cast to a Cube::ClusterState.
    const Cube::ClusterState* state
        = dynamic_cast<const Cube::ClusterState*>(fState);
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetEDeposit();
}

double Cube::ReconCluster::GetEDepositVariance() const {
    // I'm being a bit pedantic and casting to the base mix-in class.  This
    // could just as well cast to a Cube::ClusterState.
    const Cube::ClusterState* state
        = dynamic_cast<const Cube::ClusterState*>(fState);
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetEDepositVariance();
}

TLorentzVector Cube::ReconCluster::GetPosition() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ClusterState> state = GetState();
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetPosition();
}

TLorentzVector Cube::ReconCluster::GetPositionVariance() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ClusterState> state = GetState();
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetPositionVariance();
}

bool Cube::ReconCluster::IsXCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.X())) return false;
    return true;
}

bool Cube::ReconCluster::IsYCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Y())) return false;
    return true;
}

bool Cube::ReconCluster::IsZCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Z())) return false;
    return true;
}

int Cube::ReconCluster::GetDimensions() const{
    TLorentzVector var = GetPositionVariance();
    int dim = 0;
    if (IsXCluster()) ++dim;
    if (IsYCluster()) ++dim;
    if (IsZCluster()) ++dim;
    return dim;
}

const Cube::ReconCluster::MomentMatrix& Cube::ReconCluster::GetMoments() const {
    return fMoments;
}

void Cube::ReconCluster::SetMoments(double xx, double yy, double zz,
                                   double xy, double xz, double yz) {
    fTemporariesInitialized = false;
    fMoments(0,0) = xx;
    fMoments(1,1) = yy;
    fMoments(2,2) = zz;
    fMoments(0,1) = xy;
    fMoments(1,0) = xy;
    fMoments(0,2) = xz;
    fMoments(2,0) = xz;
    fMoments(1,2) = yz;
    fMoments(2,1) = yz;
}

void Cube::ReconCluster::SetMoments(const TMatrixT<double>& moments) {
    if (moments.GetNrows() != fMoments.GetNrows()
        || moments.GetNcols() != fMoments.GetNcols()) {
        throw std::runtime_error("Moment matrix sizes do not match.");
    }
    fTemporariesInitialized = false;
    for (int row=0; row<3; ++row) {
        for (int col=0; col<3; ++col) {
            fMoments(row,col) = moments(row,col);
        }
    }
}

void Cube::ReconCluster::UpdateFromHits() {
    fTemporariesInitialized = false;
    // Make sure there is a hit container.
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    if (!hits) return;

    // Make sure the hit container isn't empty.
    Cube::HitSelection::const_iterator beg = hits->begin();
    Cube::HitSelection::const_iterator end = hits->end();
    if (end-beg < 1) return;

    fStatus = Cube::ReconObject::kSuccess;
    fQuality = 1.0;
    fNDOF = std::max(1,int(end-beg-1));
    Cube::Handle<Cube::ClusterState> state = GetState();

    int dim = state->GetDimensions();
    TVectorT<double> vals(dim);
    TVectorT<double> sigs(dim);
    TVectorT<double> rms(dim);

    // Save the index into the state for each of the values.
    int eDep = state->GetEDepositIndex();
    int posX = state->GetXIndex();
    int posY = state->GetYIndex();
    int posZ = state->GetZIndex();
    int posT = state->GetTIndex();

    // Find the energy deposit and the average position.  The position
    // averages are energy deposition weighted.
    TVectorT<double> stateValues(dim);
    TVectorT<double> stateNorms(dim);
    for (Cube::HitSelection::const_iterator h = beg;
         h != end;
         ++h) {
        vals(eDep) = (*h)->GetCharge();
        sigs(eDep) = 1.0;

        vals(posX) = (*h)->GetPosition().X();
        sigs(posX) = (*h)->GetUncertainty().X();

        vals(posY) = (*h)->GetPosition().Y();
        sigs(posY) = (*h)->GetUncertainty().Y();

        vals(posZ) = (*h)->GetPosition().Z();
        sigs(posZ) = (*h)->GetUncertainty().Z();

        vals(posT) = (*h)->GetTime();
        sigs(posT) = (*h)->GetTimeUncertainty();

        // Form the sums needed to find averages.  Note that the energy
        // deposit is not averaged,  It's taken out of the sum.
        stateValues(eDep) += vals(eDep);
        stateNorms(eDep) = 1.0;
        for (int i=0; i<dim; ++i) {
            if (i == eDep) continue;
#define ENERGY_WEIGHTED
#ifdef ENERGY_WEIGHTED
            stateValues(i) += vals(eDep)*vals(i)/(sigs(i)*sigs(i));
            stateNorms(i) += vals(eDep)*1.0/(sigs(i)*sigs(i));
#else
            stateValues(i) += vals(i)/(sigs(i)*sigs(i));
            stateNorms(i) += 1.0/(sigs(i)*sigs(i));
#endif
        }
    }

    // Reset the edep norm so it is handled with everything else, then form
    // the averages.
    stateNorms(eDep) = 1.0;
    for (int i=0; i<dim; ++i) {
        if (stateNorms(i) > 0) stateValues(i) /= stateNorms(i);
    }

    // Estimate the covariance of cluster state values.
    TMatrixTSym<double> stateCov(dim);
    // This counts the weight for each bin of the covariance.
    TMatrixTSym<double> weights(dim);
    TMatrixTSym<double> weights2(dim);
    // This counts the number of degrees of freedom contributing to each bin.
    TMatrixTSym<double> dof(dim);
    for (Cube::HitSelection::const_iterator h = beg;
         h != end;
         ++h) {
        vals(eDep) = (*h)->GetCharge() - stateValues(eDep);
        sigs(eDep) = std::sqrt((*h)->GetCharge());
        rms(eDep) = 0.0;

        vals(posX) = (*h)->GetPosition().X() - stateValues(posX);
        sigs(posX) = (*h)->GetUncertainty().X();
        rms(posX) = (*h)->GetSize().X();

        vals(posY) = (*h)->GetPosition().Y() - stateValues(posY);
        sigs(posY) = (*h)->GetUncertainty().Y();
        rms(posY) = (*h)->GetSize().Y();

        vals(posZ) = (*h)->GetPosition().Z() - stateValues(posZ);
        sigs(posZ) = (*h)->GetUncertainty().Z();
        rms(posZ) = (*h)->GetSize().Z();

        vals(posT) = (*h)->GetTime() - stateValues(posT);
        sigs(posT) = (*h)->GetTimeUncertainty();
        rms(posT) = (*h)->GetTimeUncertainty();

        for (int row = 0; row<dim; ++row) {
            for (int col = row; col<dim; ++col) {
                double weight = 1.0/(sigs(row)*sigs(col));
                stateCov(row,col) += weight*vals(row)*vals(col);
                weights(row,col) += weight;
                weights2(row,col) += weight*weight;
                double degrees
                    = 4*rms(row)*rms(col)/(12*sigs(row)*sigs(col));
                if (row == posT && col == posT) degrees = 1.0;
                dof(row,col) += degrees;
            }
        }
    }

    // Turn the "stateCov" variable into the RMS.
    for (int row = 0; row<dim; ++row) {
        for (int col = row; col<dim; ++col) {
            if (weights(row,col)>0) {
                stateCov(row,col) /= weights(row,col);
            }
            else stateCov(row,col) = 0.0;
            stateCov(col,row) = stateCov(row,col);
        }
    }

    // Turn the RMS into a covariance of the mean (This is almost a repeat of
    // turning the value into an RMS.
    for (int row = 0; row<dim; ++row) {
        for (int col = row; col<dim; ++col) {
            if (dof(row,col)>0.9) {
#ifdef NEW_IMPLEMENTATION
                double w = weights(row,col);
                stateCov(row,col) *= w*w / (w*w - weights2(row,col));
#else
                stateCov(row,col) /= std::sqrt(dof(row,col));
#endif
            }
            else if (row==col) stateCov(row,col) = Cube::CorrValues::kFreeValue;
            else stateCov(row,col) = 0.0;
            stateCov(col,row) = stateCov(row,col);
        }
    }

#ifndef NO_FINITE_SIZE_CORRECTION
    // Add the correction for finite size of the hits.
    TVectorT<double> hitWeights(dim);
    for (Cube::HitSelection::const_iterator h = beg;
         h != end;
         ++h) {
        sigs(eDep) = std::sqrt((*h)->GetCharge());
        sigs(posX) = (*h)->GetUncertainty().X();
        sigs(posY) = (*h)->GetUncertainty().Y();
        sigs(posZ) = (*h)->GetUncertainty().Z();
        sigs(posT) = (*h)->GetTimeUncertainty();

        for (int idx = 0; idx<dim; ++idx) {
            double weight = 1.0/(sigs(idx)*sigs(idx));
            hitWeights(idx) += weight;
        }
    }

    for (int idx = 0; idx<dim; ++idx) {
        if (hitWeights(idx)<1E-8) continue;
        stateCov(idx,idx) += 1.0/hitWeights(idx);
    }
#endif

    // Fix the variance of the deposited energy.  This assumes it's Poisson
    // distributed.
    stateCov(eDep,eDep) = stateValues(eDep);

    // Set the state value and covariance.
    for (int row=0; row<dim; ++row) {
        state->SetValue(row,stateValues(row));
        for (int col=0; col<dim; ++col) {
            state->SetCovarianceValue(row,col,stateCov(row,col));
            state->SetCovarianceValue(col,row,stateCov(row,col));
        }
    }

    // Find the moments of the cluster.  This could be done as part of the
    // covariance calculation, but is seperated so it's easier to see what is
    // going on.
    TMatrixTSym<double> moments(3);
    TMatrixTSym<double> chargeSum(3);
    TVector3 center(state->GetValue(posX),
                    state->GetValue(posY),
                    state->GetValue(posZ));

    // Calculate the charge weighted "sum of squared differences".
    for (Cube::HitSelection::const_iterator h = beg;
         h != end;
         ++h) {
        TVector3 diff = (*h)->GetPosition() - center;
        for (int row = 0; row<3; ++row) {
            for (int col = row; col<3; ++col) {
                moments(row,col) += diff(row)*diff(col)*(*h)->GetCharge();
                // Correct for the finite size of the hit.
                if (row==col) {
                    double r = (*h)->GetSize()[row];
                    moments(row,col) += r*r*(*h)->GetCharge();
                }
                chargeSum(row,col) += (*h)->GetCharge();
            }
        }
    }

    // Turn this into the average square difference (i.e. the moments)
    for (int row = 0; row<3; ++row) {
        for (int col = row; col<3; ++col) {
            // No divide by zero please!
            if (chargeSum(row,col) > 1E-6) {
                moments(row,col) = moments(row,col)/chargeSum(row,col);
            }
            else if (row == col) {
                // There were no measurements on this axis so spread the
                // charge over the entire range of the detector
                // (i.e. +-"10 m + epsilon")
                moments(row,col) = 1E+9;
            }
            else {
                // There were no measurements so the aren't any off-axis
                // correlations.
                moments(row,col) = 0;
            }
            // Keep it symmetric.
            moments(col,row) = moments(row,col);
        }
    }

    SetMoments(moments);
}

void Cube::ReconCluster::ls(Option_t *opt) const {
    Cube::ReconObject::ls_base(opt);
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();

    if (fState) {
        TROOT::IncreaseDirLevel();
        fState->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    std::ios::fmtflags save = std::cout.flags();
    for (int i = 0; i<3; ++i) {
        TROOT::IndentLevel();
        if (i == 0) std::cout << "Position Cov: ";
        else        std::cout << "              ";
        Cube::Handle<Cube::ClusterState> state = GetState();
        for (int j = 0; j<3; ++j) {
            std::cout << std::setw(12) << state->GetPositionCovariance(i,j);
        }
        std::cout << std::setw(0) << std::endl;
    }
    for (int i = 0; i<3; ++i) {
        TROOT::IndentLevel();
        if (i == 0) std::cout << "Moments:      ";
        else        std::cout << "              ";
        for (int j = 0; j<3; ++j) {
            std::cout << std::setw(12) << fMoments(i,j);
        }
        std::cout << std::setw(0) << std::endl;
    }
    std::cout.flags(save);

    if (fNodes) {
        TROOT::IncreaseDirLevel();
        fNodes->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    TROOT::DecreaseDirLevel();
}

void Cube::ReconCluster::FillTemporaries() const {
    if (fTemporariesInitialized) return;
    fTemporariesInitialized = true;

    TVectorF eigenValues;
    TMatrixF eigenVectors(GetMoments().EigenVectors(eigenValues));

    // Long axis comes first.
    fLongAxis.SetXYZ(eigenVectors(0,0), eigenVectors(1,0), eigenVectors(2,0));
    fLongAxis = std::sqrt(eigenValues(0))*fLongAxis;

    // Major axis comes second.
    fMajorAxis.SetXYZ(eigenVectors(0,1), eigenVectors(1,1), eigenVectors(2,1));
    fMajorAxis = std::sqrt(eigenValues(1))*fMajorAxis;

    // Minor axis comes last.
    fMinorAxis.SetXYZ(eigenVectors(0,2), eigenVectors(1,2), eigenVectors(2,2));
    fMinorAxis = std::sqrt(eigenValues(2))*fMinorAxis;

    const double epsilon = 1E-6;
    if (fLongAxis.X() < -epsilon) {
        fLongAxis = -fLongAxis;
    }
    else if (fLongAxis.X() < epsilon && fLongAxis.Y() < -epsilon) {
        fLongAxis = -fLongAxis;
    }
    else if (fLongAxis.Y() < epsilon && fLongAxis.Z() < 0.0) {
        fLongAxis = -fLongAxis;
    }

    if (fMajorAxis.Y() < -epsilon) {
        fMajorAxis = -fMajorAxis;
    }
    else if (fMajorAxis.Y() < epsilon && fMajorAxis.Z() < -epsilon) {
        fMajorAxis = -fMajorAxis;
    }
    else if (fMajorAxis.Z() < epsilon && fMajorAxis.X() < 0.0) {
        fMajorAxis = -fMajorAxis;
    }

    TVector3 temp = fLongAxis.Cross(fMajorAxis);
    if (temp*fMinorAxis < 0) {
        fMinorAxis = -fMinorAxis;
    }

}

const TVector3& Cube::ReconCluster::GetLongAxis() const {
    FillTemporaries();
    return fLongAxis;
}

const TVector3& Cube::ReconCluster::GetMajorAxis() const {
    FillTemporaries();
    return fMajorAxis;
}

const TVector3& Cube::ReconCluster::GetMinorAxis() const {
    FillTemporaries();
    return fMinorAxis;
}

double Cube::ReconCluster::GetLongExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fLongAxis));
    }
    maxLen /= fLongAxis.Mag();
    return maxLen;

}

double Cube::ReconCluster::GetMajorExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fMajorAxis));
    }
    maxLen /= fMajorAxis.Mag();
    return maxLen;

}

double Cube::ReconCluster::GetMinorExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fMinorAxis));
    }
    maxLen /= fMinorAxis.Mag();
    return maxLen;

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
