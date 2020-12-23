#include "TReconTrackElement.hxx"
#include "TReconClusterElement.hxx"
#include "TMatrixElement.hxx"
#include "TEventDisplay.hxx"
#include "TGUIManager.hxx"

#include <CubeUnits.hxx>
#include <CubeHandle.hxx>
#include <TUnitsTable.hxx>

#include <CubeReconTrack.hxx>
#include <CubeReconNode.hxx>
#include <CubeReconCluster.hxx>
#include <CubeTrackState.hxx>
#include <CubeHandle.hxx>

#include <TEveLine.h>

#include <sstream>

Cube::TReconTrackElement::~TReconTrackElement() {}

Cube::TReconTrackElement::TReconTrackElement(Cube::ReconTrack& track,
                                             bool showStates,
                                             bool showObjects,
                                             bool showDirection)
    : TEveElementList(), fValid(true) {

    Cube::Handle<Cube::TrackState> frontState = track.GetState();
    if (!frontState) {
        CUBE_ERROR << "TrackState missing!" << std::endl;
        fValid = false;
        return;
    }

#define USE_BACKSTATE
#ifdef USE_BACKSTATE
    Cube::Handle<Cube::TrackState> backState
        = track.GetNodes().back()->GetState();
#endif

    // The unit length for drawing the line directions.  The direction vectors
    // are this long.
    const double dirLineLen = 7.0*unit::mm;

    //////////////////////////////////////////////////////////
    // Generate the title for the track.
    //////////////////////////////////////////////////////////
    std::ostringstream title;
    {
        TLorentzVector pos = frontState->GetPosition();
        TLorentzVector var = frontState->GetPositionVariance();
        TVector3 dir = frontState->GetDirection().Unit();
        TVector3 dvar = frontState->GetDirectionVariance();

        double length = GetLength(track);
        double dEdX = 0.0;
        if (length > 0) {
            dEdX = track.GetEDeposit()/length;
        }
        title << "Track(" << track.GetUniqueID() << ") --"
              << " Nodes: " << track.GetNodes().size()
              << " Hits: " << track.GetHitSelection()->size()
              << " Length: " << unit::AsString(length,"length");
        title << std::endl
              << "   Energy Deposit: "
              << unit::AsString(track.GetEDeposit(),"pe")
              << ", dEdX: " << unit::AsString(dEdX,"pe") << "/mm";
        title << ", Qual: " << track.GetQuality() << " / " << track.GetNDOF();
        title << std::endl
              << "   Front Pos: ("
              << unit::AsString(pos.X(),std::sqrt(var.X()),"length")
              << ", "<<unit::AsString(pos.Y(),std::sqrt(var.Y()),"length")
              << ", "<<unit::AsString(pos.Z(),std::sqrt(var.Z()),"length")
              << ")";
        title << std::endl
              << "   Front Dir: ("
              << unit::AsString(dir.X(), dvar.X(),"direction")
              << ", " << unit::AsString(dir.Y(), dvar.Y(),"direction")
              << ", " << unit::AsString(dir.Z(), dvar.Z(),"direction") << ")"
              << " T: " << unit::AsString(pos.T(),"time");

    }

#ifdef USE_BACKSTATE
    if (backState) {
        TLorentzVector v = backState->GetPositionVariance();
        TLorentzVector p = backState->GetPosition();
        TVector3 d = backState->GetDirection().Unit();
        TVector3 dv = backState->GetDirectionVariance();
        title << std::endl
              << "   Back Pos: ("
              << unit::AsString(p.X(),std::sqrt(v.X()),"length")
              <<", "<<unit::AsString(p.Y(),std::sqrt(v.Y()),"length")
              <<", "<<unit::AsString(p.Z(),std::sqrt(v.Z()),"length") << ")";
        title << std::endl
              << "   Back Dir: ("
              << unit::AsString(d.X(), dv.X(),"direction")
              << ", " << unit::AsString(d.Y(), dv.Y(),"direction")
              << ", " << unit::AsString(d.Z(), dv.Z(),"direction") << ")"
              << " T: " << unit::AsString(p.T(),"time");
    }
    else {
        title << std::endl
              << "      BACK STATE IS MISSING";
    }
#endif

    Cube::ReconNodeContainer& nodes = track.GetNodes();
#define DUMP_INFO_TO_TITLE
#ifdef DUMP_INFO_TO_TITLE
    int maxNodes = 10;
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        if (maxNodes-- < 1) break;
        Cube::Handle<Cube::TrackState> nodeState = (*n)->GetState();
        Cube::Handle<Cube::ReconCluster> nodeObject = (*n)->GetObject();
        title << std::endl << "N: " << nodeState->GetPosition().T();
        if (nodeObject) {
            title << " M: " << nodeObject->GetPosition().T()
                  << "+/-" << std::sqrt(nodeObject->GetPositionVariance().T());
            Cube::Handle<Cube::HitSelection> hits
                = nodeObject->GetHitSelection();
            if (!hits) continue;
            int maxHits = 6;
            title << " --";
            for (Cube::HitSelection::iterator h = hits->begin();
                 h != hits->end(); ++h) {
                if (maxHits-- <1) break;
                int t = 2.0*(*h)->GetTime() + 0.5;
                int q = (*h)->GetCharge() + 0.5;
                title << " "  << 0.5*t << "(" << q
                      << "," <<(*h)->GetTimeUncertainty()<< ")";
            }
        }
    }
#endif

    CUBE_LOG(0) << "track:: " << title.str() << std::endl;

    std::ostringstream objName;
    objName << track.GetName() << "(" << track.GetUniqueID() << ")";

    SetMainColor(kBlue);
    SetName(objName.str().c_str());
    SetTitle(title.str().c_str());

    ////////////////////////////////////////////////////////////////
    // Draw the line connecting all of the nodes.
    ////////////////////////////////////////////////////////////////

    std::unique_ptr<TEveLine> trackLine(new TEveLine(nodes.size()));

    trackLine->SetName(objName.str().c_str());

    trackLine->SetTitle(title.str().c_str());
    trackLine->SetSourceObject(&track);
    trackLine->SetLineStyle(1);
    trackLine->SetLineWidth(4);

    int color = kBlue;
    if (GetLength(track) > 0.5*unit::mm) {
        double dEdX = track.GetEDeposit()/GetLength(track);
        double minEnergy = 5.0;
        double maxEnergy = 50.0;
        color = TEventDisplay::Get().LinearColor(dEdX,
                                              minEnergy,
                                              maxEnergy);
    }
    trackLine->SetLineColor(color);

    int p=0;

    // Start at the front state of the track
    if (frontState) {
        TLorentzVector frontPos = frontState->GetPosition();
        TLorentzVector frontVar = frontState->GetPositionVariance();
        if (std::isfinite(frontPos.X())
            && std::isfinite(frontPos.Y())
            && std::isfinite(frontPos.Z())) {
            trackLine->SetPoint(p++, frontPos.X(), frontPos.Y(), frontPos.Z());
        }
        else fValid = false;
        CUBE_DEBUG(2) << "nodes::" << "Front:"
                  << unit::AsString(frontPos.X(),
                                    std::sqrt(frontVar.X()),"length")
                  <<", "<<unit::AsString(frontPos.Y(),
                                         std::sqrt(frontVar.Y()),"length")
                  <<", "<<unit::AsString(frontPos.Z(),
                                         std::sqrt(frontVar.Z()),"length")
                  << std::endl;
        CUBE_DEBUG(2) << "nodes"
                  << "Front Dir: "
                  << unit::AsString(frontState->GetDirection())
                  << std::endl;
    }

    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        Cube::Handle<Cube::TrackState> nodeState = (*n)->GetState();
        Cube::Handle<Cube::ReconObject> nodeObject = (*n)->GetObject();
        if (!nodeState) continue;
        TLorentzVector nodePos = nodeState->GetPosition();
        TLorentzVector nodeVar = nodeState->GetPositionVariance();
        if (std::isfinite(nodePos.X())
            && std::isfinite(nodePos.Y())
            && std::isfinite(nodePos.Z())) {
            trackLine->SetPoint(p++, nodePos.X(), nodePos.Y(), nodePos.Z());
        }
        else {
            fValid = false;
        }
        CUBE_DEBUG(2) << "nodes"
                  << "Pos:"
                  << unit::AsString(nodePos.X(),
                                    std::sqrt(nodeVar.X()),"length")
                  <<", "<<unit::AsString(nodePos.Y(),
                                         std::sqrt(nodeVar.Y()),"length")
                  <<", "<<unit::AsString(nodePos.Z(),
                                         std::sqrt(nodeVar.Z()),"length")
                  << std::endl;
        Cube::Handle<Cube::ReconCluster> cluster = nodeObject;
        if (cluster) {
            double delta = (cluster->GetPosition().Vect()-nodePos.Vect()).Mag();
            CUBE_DEBUG(2) << "nodes::" << "Cluster: "
                      << unit::AsString(cluster->GetPosition().Vect(),
                                        "length")
                      << "  diff: " << unit::AsString(delta,"length")
                      << std::endl;
        }
        CUBE_DEBUG(2) << "nodes::"
                  << "Dir: " << unit::AsString(nodeState->GetDirection())
                  << std::endl;
    }

#ifdef USE_BACKSTATE
    // finish at the back state of the track
    if (backState) {
        TLorentzVector backPos = backState->GetPosition();
        TLorentzVector backVar = backState->GetPositionVariance();
        if (std::isfinite(backPos.X())
            && std::isfinite(backPos.Y())
            && std::isfinite(backPos.Z())) {
            trackLine->SetPoint(p++, backPos.X(), backPos.Y(), backPos.Z());
        }
        else {
            fValid = false;
        }
        CUBE_DEBUG(2) << "nodes::" << "Back:"
                  << unit::AsString(backPos.X(),
                                    std::sqrt(backVar.X()),"length")
                  <<", "<<unit::AsString(backPos.Y(),
                                         std::sqrt(backVar.Y()),"length")
                  <<", "<<unit::AsString(backPos.Z(),
                                         std::sqrt(backVar.Z()),"length")
                  << std::endl;
        CUBE_DEBUG(2) << "nodes::"
                  << "Back Dir: "
                  << unit::AsString(backState->GetDirection())
                  << std::endl;
    }
#endif
    AddElement(trackLine.release());

    //////////////////////////////////////////////////////////
    // Draw the state for the track and nodes.  (Show both the position and
    //  uncertainty)
    //////////////////////////////////////////////////////////

    if (showStates && frontState) {
        TLorentzVector nodePos = frontState->GetPosition();
        TMatrixD nodeVar(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                nodeVar(i,j) = frontState->GetPositionCovariance(i,j);
            }
        }
        std::unique_ptr<Cube::TMatrixElement> uncertainty(
            new Cube::TMatrixElement("Uncertainty",
                                   nodePos.Vect(),
                                   nodeVar,
                                   false));
        uncertainty->SetMainColor(kSpring-9);
        uncertainty->SetSourceObject(&(*frontState));
        if (uncertainty->Valid()) AddElement(uncertainty.release());
    }

    // Handle the states for the nodes.
    if (showStates) {
        int cnt = 0;
        for (Cube::ReconNodeContainer::iterator n = nodes.begin();
             n != nodes.end(); ++n) {
            ++cnt;
            Cube::Handle<Cube::TrackState> nodeState = (*n)->GetState();
            Cube::Handle<Cube::ReconObject> nodeObject = (*n)->GetObject();
            if (!nodeState) {
                CUBE_ERROR << "Node state is missing" << std::endl;
                continue;
            }

            // Create a "tube" that shows the position and position
            // uncertainty for the state.
            TLorentzVector nodePos = nodeState->GetPosition();
            TMatrixD nodeVar(3,3);
            for (int i=0; i<3; ++i) {
                for (int j=0; j<3; ++j) {
                    nodeVar(i,j) = nodeState->GetPositionCovariance(i,j);
                }
            }
            std::unique_ptr<Cube::TMatrixElement> stateElement(
                new Cube::TMatrixElement("Uncertainty",
                                       nodePos.Vect(),
                                       nodeVar,
                                       false));

            // Guestimate the length of track being described by this node.
            double length = 1.0;
            if (n == nodes.begin()) {
                Cube::Handle<Cube::TrackState> b = (*(n+1))->GetState();
                length = 0.75*(frontState->GetPosition().Vect()
                               - b->GetPosition().Vect()).Mag();
            }
            else if ((n+1) == nodes.end()) {
                Cube::Handle<Cube::TrackState> b = (*(n-1))->GetState();
                Cube::Handle<Cube::TrackState> e = (*(n))->GetState();
                length = 0.75*(e->GetPosition().Vect()
                               - b->GetPosition().Vect()).Mag();
            }
            else {
                Cube::Handle<Cube::TrackState> a = (*(n-1))->GetState();
                Cube::Handle<Cube::TrackState> b = (*(n+1))->GetState();
                length = 0.5*(a->GetPosition().Vect()
                              - b->GetPosition().Vect()).Mag();
            }

            // Figure out a color to use based on the measured charge.
            int color = kBlue;
            double energy
                = Cube::TEventDisplay::Get().CrudeEnergy(
                    nodeState->GetEDeposit());
            double dEdX = energy;
            if (length > 0.5*unit::mm) {
                dEdX /= length;              // Get charge per length;
                double minEnergy = 5.0;
                double maxEnergy = 30.0;
                color = TEventDisplay::Get().LogColor(dEdX,
                                                      minEnergy,
                                                      maxEnergy,2.0);
            }
            stateElement->SetMainColor(color);
            stateElement->SetSourceObject(&(*nodeState));
            if (stateElement->Valid()) AddElement(stateElement.release());

            // Show the object that is attached to the node.
            Cube::Handle<Cube::ReconCluster> cluster = nodeObject;
            if (showObjects) {
                Cube::TReconClusterElement* objectElement
                    = new Cube::TReconClusterElement((*cluster), true);
                AddElement(objectElement);
            }
            // Add a line connecting the state and the object.
            TEveLine* objectLine = new TEveLine(2);
            objectLine->SetLineColor(kBlue-7);
            objectLine->SetLineStyle(1);
            objectLine->SetLineWidth(2);
            objectLine->SetPoint(0,
                                 nodeState->GetPosition().X(),
                                 nodeState->GetPosition().Y(),
                                 nodeState->GetPosition().Z());
            objectLine->SetPoint(1,
                                 cluster->GetPosition().X(),
                                 cluster->GetPosition().Y(),
                                 cluster->GetPosition().Z());
            AddElement(objectLine);
        }
    }

#ifdef USE_BACKSTATE
    // Add the back state position and position uncertainty.
    if (showStates && backState) {
        TLorentzVector nodePos = backState->GetPosition();
        TMatrixD nodeVar(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                nodeVar(i,j) = backState->GetPositionCovariance(i,j);
            }
        }
        std::unique_ptr<Cube::TMatrixElement> uncertainty(
            new Cube::TMatrixElement("Uncertainty",
                                   nodePos.Vect(),
                                   nodeVar,
                                   false));
        uncertainty->SetMainColor(kSpring-9);
        uncertainty->SetSourceObject(&(*backState));
        if (uncertainty->Valid()) AddElement(uncertainty.release());
    }
#endif

    /////////////////////////////////////////////////////////////
    // Draw the direction at each node.
    /////////////////////////////////////////////////////////////
#define NODE_DIRECTION
#ifdef NODE_DIRECTION
    if (showDirection) {
        // Add the node direction information.
        for (Cube::ReconNodeContainer::iterator n = nodes.begin();
             n != nodes.end(); ++n) {
            Cube::Handle<Cube::TrackState> nodeState = (*n)->GetState();
            if (!nodeState) {
                CUBE_ERROR << "Node state is missing" << std::endl;
                continue;
            }
            TLorentzVector nodePos = nodeState->GetPosition();
            TVector3 nodeDir = nodeState->GetDirection();
            TEveLine* directionLine = new TEveLine(2);
            directionLine->SetLineColor(kRed);
            directionLine->SetLineStyle(1);
            directionLine->SetLineWidth(1);
            if (std::isfinite(nodePos.X())
                && std::isfinite(nodePos.Y())
                && std::isfinite(nodePos.Z())
                && std::isfinite(nodeDir.X())
                && std::isfinite(nodeDir.Y())
                && std::isfinite(nodeDir.Z())) {
                directionLine->SetPoint(0,
                                        nodePos.X(),
                                        nodePos.Y(),
                                        nodePos.Z());
                directionLine->SetPoint(1,
                                        nodePos.X()+dirLineLen*nodeDir.X(),
                                        nodePos.Y()+dirLineLen*nodeDir.Y(),
                                        nodePos.Z()+dirLineLen*nodeDir.Z());
                AddElement(directionLine);
            }
            else fValid = false;
        }
    }
#endif

    /////////////////////////////////////////////////////////////
    // Draw the direction at the front of the track.
    /////////////////////////////////////////////////////////////
#define FRONT_DIRECTION
#ifdef FRONT_DIRECTION
    if (showDirection && frontState) {
        TLorentzVector frontPos = frontState->GetPosition();
        TVector3 frontDir = frontState->GetDirection();
        TEveLine* directionLine = new TEveLine(2);
        directionLine->SetLineColor(kMagenta-9);
        directionLine->SetLineStyle(1);
        directionLine->SetLineWidth(1);
        TVector3 tipPos = frontPos.Vect() - dirLineLen*frontDir;
        if (std::isfinite(frontPos.X())
            && std::isfinite(frontPos.Y())
            && std::isfinite(frontPos.Z())
            && std::isfinite(tipPos.X())
            && std::isfinite(tipPos.Y())
            && std::isfinite(tipPos.Z())) {
            directionLine->SetPoint(0, frontPos.X(), frontPos.Y(), frontPos.Z());
            directionLine->SetPoint(1, tipPos.X(), tipPos.Y(), tipPos.Z());
            AddElement(directionLine);
        }
        else fValid = false;

        TMatrixD tipVar(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                tipVar(i,j)
                    = dirLineLen*frontState->GetDirectionCovariance(i,j);
            }
        }
        std::unique_ptr<Cube::TMatrixElement> uncertainty(
            new Cube::TMatrixElement("Uncertainty",
                                   tipPos,
                                   tipVar,
                                   false));
        uncertainty->SetMainColor(kMagenta-9);
        uncertainty->SetSourceObject(&(*frontState));
        if (uncertainty->Valid()) AddElement(uncertainty.release());
    }
#endif

#define BACK_DIRECTION
#ifdef BACK_DIRECTION
    /////////////////////////////////////////////////////////////
    // Draw the direction at the back of the track (if it exists).
    /////////////////////////////////////////////////////////////
    if (showDirection && backState) {
        TLorentzVector backPos = backState->GetPosition();
        TVector3 backDir = backState->GetDirection();
        TEveLine* directionLine = new TEveLine(2);
        directionLine->SetLineColor(kMagenta-9);
        directionLine->SetLineStyle(1);
        directionLine->SetLineWidth(1);
        TVector3 tipPos = backPos.Vect() + dirLineLen*backDir;
        if (std::isfinite(backPos.X())
            && std::isfinite(backPos.Y())
            && std::isfinite(backPos.Z())
            && std::isfinite(tipPos.X())
            && std::isfinite(tipPos.Y())
            && std::isfinite(tipPos.Z())) {
            directionLine->SetPoint(0, backPos.X(), backPos.Y(), backPos.Z());
            directionLine->SetPoint(1, tipPos.X(), tipPos.Y(), tipPos.Z());
            AddElement(directionLine);
        }
        else fValid = false;

        TMatrixD tipVar(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                tipVar(i,j) = backState->GetDirectionCovariance(i,j);
            }
        }
        std::unique_ptr<Cube::TMatrixElement> uncertainty(
            new Cube::TMatrixElement("Uncertainty",
                                   tipPos,
                                   tipVar,
                                   false));
        uncertainty->SetMainColor(kMagenta-9);
        uncertainty->SetSourceObject(&(*backState));
        if (uncertainty->Valid()) AddElement(uncertainty.release());
    }
#endif

}

double Cube::TReconTrackElement::GetLength(Cube::ReconTrack& track) const {
    double length = 0.0;
    Cube::ReconNodeContainer::iterator n = track.GetNodes().begin();
    Cube::Handle<Cube::TrackState> lastState = (*(n++))->GetState();
    while (n != track.GetNodes().end()) {
        Cube::Handle<Cube::TrackState> nodeState = (*(n++))->GetState();
        length += (nodeState->GetPosition().Vect()
                   - lastState->GetPosition().Vect()).Mag();
        lastState = nodeState;
    }
    return length;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
