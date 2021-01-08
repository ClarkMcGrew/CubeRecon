#include "TReconClusterElement.hxx"
#include "TMatrixElement.hxx"
#include "TEventDisplay.hxx"
#include "TGUIManager.hxx"

#include <CubeUnits.hxx>
#include <CubeHandle.hxx>
#include <TUnitsTable.hxx>

#include <CubeReconCluster.hxx>
#include <CubeClusterState.hxx>
#include <CubeHitSelection.hxx>
#include <CubeHandle.hxx>

#include <TEveLine.h>

#include <sstream>

Cube::TReconClusterElement::~TReconClusterElement() {}

Cube::TReconClusterElement::TReconClusterElement(Cube::ReconCluster& cluster,
                                                 bool showUncertainty,
                                                 bool printInfo)
    : TEveElementList(), fValid(true) {

    double minEnergy = 10.0;
    double maxEnergy = 200.0;

    Cube::Handle<Cube::ClusterState> state = cluster.GetState();
    TLorentzVector var = state->GetPositionVariance();
    TLorentzVector pos = state->GetPosition();
    if (printInfo) {
        std::cout << "Cluster(" << cluster.GetUniqueID() << ") @ "
                  << unit::AsString(pos.X(),std::sqrt(var.X()),"length")
                  << ", "
                  << unit::AsString(pos.Y(),std::sqrt(var.Y()),"length")
                  << ", "
                  << unit::AsString(pos.Z(),std::sqrt(var.Z()),"length")
                  << std::endl;

        std::cout << "  Time "
                  << unit::AsString(pos.T(),std::sqrt(var.T()),"time")
                  << "  Charge "
                  << unit::AsString(cluster.GetEDeposit(),"pe")
                  << std::endl;
    }
    // Find the rotation of the object to be displayed.  A cluster is
    // represented as a tube with the long axis along the local Z direction,
    // and the major and minor in the XY plane.
    TMatrixD tubeRot(3,3);
    if (showUncertainty) {
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                tubeRot(i,j) = state->GetPositionCovariance(i,j);
            }
        }
    }
    else {
        const Cube::ReconCluster::MomentMatrix& moments = cluster.GetMoments();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                tubeRot(i,j) = moments(i,j);
            }
        }
    }

    ///////////////////////////////////////////
    // Fill the object to be drawn.
    ///////////////////////////////////////////

    std::unique_ptr<Cube::TMatrixElement> eveCluster(
        new Cube::TMatrixElement("cluster",
                               cluster.GetPosition().Vect(),
                               tubeRot,
                               true));

    // Set the name.
    std::ostringstream name;
    name << "cluster(" << cluster.GetUniqueID() << ")";
    eveCluster->SetName(name.str().c_str());

    // Set the color.
    int color = kCyan-9;
    // EDeposit is in measured charge.
    double energy = cluster.GetEDeposit()/cluster.GetHitSelection()->size();
    color = TEventDisplay::Get().LogColor(energy,
                                          minEnergy, maxEnergy,2.0);
    eveCluster->SetMainColor(color);

    bool transparentClusters = true;
    if (transparentClusters) eveCluster->SetMainTransparency(60);
    else eveCluster->SetMainTransparency(0);

    // Build the cluster title.
    std::ostringstream title;
    title << "Cluster(" << cluster.GetUniqueID() << ") @ ";
    title << unit::AsString(pos.T(),"time");
#ifdef ADD_POSITION_TO_TITLE
    title << unit::AsString(pos.X(),std::sqrt(var.X()),"length")
          << ", " << unit::AsString(pos.Y(),std::sqrt(var.Y()),"length")
          << ", " << unit::AsString(pos.Z(),std::sqrt(var.Z()),"length");
#endif
    title << " w/ " << unit::AsString(energy,-1,"pe") << "/Cube";
    title << " in " << cluster.GetHitSelection()->size() << " cubes";

    SetName(name.str().c_str());
    SetTitle(title.str().c_str());

    eveCluster->SetTitle(title.str().c_str());
    eveCluster->SetSourceObject(&cluster);
    if (eveCluster->Valid()) AddElement(eveCluster.release());

    // Add the connections between the hits!
    Cube::Handle<Cube::HitSelection> hits = cluster.GetHitSelection();
    if (hits->size() < 2) return;

    bool connectEnds
        = Cube::TEventDisplay::Get().GUI().GetShowFitEndsButton()->IsOn();
    bool connectHits
        = Cube::TEventDisplay::Get().GUI().GetShowClusterConnectionsButton()
        ->IsOn();
    if (!connectEnds && !connectHits) return;

#ifdef CLUSTERS_AS_TRACKS
    int points = 1;
    if (connectEnds) ++points;
    if (connectHits) points += hits->size() - 1;

    TEveLine* eveLine = new TEveLine(points);
    eveLine->SetName(name.str().c_str());
    eveLine->SetTitle(title.str().c_str());
    eveLine->SetSourceObject(&cluster);
    eveLine->SetLineColor(color);
    eveLine->SetLineStyle(1);
    eveLine->SetLineWidth(2);

    int p = 0;
    if (connectEnds) {
        eveLine->SetPoint(p++,
                          hits->back()->GetPosition().X(),
                          hits->back()->GetPosition().Y(),
                          hits->back()->GetPosition().Z());
    }

    eveLine->SetPoint(p++,
                      hits->front()->GetPosition().X(),
                      hits->front()->GetPosition().Y(),
                      hits->front()->GetPosition().Z());

    if (connectHits) {
        for (Cube::HitSelection::iterator h = hits->begin()+1;
             h != hits->end(); ++h) {
            eveLine->SetPoint(p++,
                              (*h)->GetPosition().X(),
                             (*h)->GetPosition().Y(),
                              (*h)->GetPosition().Z());
        }
    }

    AddElement(eveLine);
#else
    if (connectEnds) {
            TEveLine* eveLine = new TEveLine(2);
            eveLine->SetName(name.str().c_str());
            eveLine->SetTitle(title.str().c_str());
            eveLine->SetSourceObject(&cluster);
            eveLine->SetLineColor(color);
            eveLine->SetLineStyle(1);
            eveLine->SetLineWidth(2);
            eveLine->SetPoint(0,
                              hits->front()->GetPosition().X(),
                              hits->front()->GetPosition().Y(),
                              hits->front()->GetPosition().Z());
            eveLine->SetPoint(1,
                              hits->back()->GetPosition().X(),
                              hits->back()->GetPosition().Y(),
                              hits->back()->GetPosition().Z());
            AddElement(eveLine);
    }
    if (connectHits) {
        for (Cube::HitSelection::iterator h = hits->begin();
             h != hits->end(); ++h) {
            TEveLine* eveLine = new TEveLine(2);
            eveLine->SetName(name.str().c_str());
            eveLine->SetTitle(title.str().c_str());
            eveLine->SetSourceObject(&cluster);
            eveLine->SetLineColor(color);
            eveLine->SetLineStyle(1);
            eveLine->SetLineWidth(2);
            eveLine->SetPoint(0, pos.X(), pos.Y(), pos.Z());
            eveLine->SetPoint(1,
                              (*h)->GetPosition().X(),
                              (*h)->GetPosition().Y(),
                              (*h)->GetPosition().Z());
            AddElement(eveLine);
        }
    }
#endif
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
