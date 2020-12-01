#include "TFitChangeHandler.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TGUIManager.hxx"
#include "TShowHits.hxx"
#include "TReconTrackElement.hxx"
#include "TReconShowerElement.hxx"
#include "TReconClusterElement.hxx"

#include <CubeHandle.hxx>
#include <CubeReconNode.hxx>
#include <CubeLog.hxx>

#include <TUnitsTable.hxx>

#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoEltu.h>
#include <TGeoSphere.h>
#include <TGeoMatrix.h>

#include <TVectorF.h>
#include <TMatrixF.h>

#include <TGButton.h>
#include <TGListBox.h>
#include <TCollection.h>

#include <TEveManager.h>
#include <TEveGeoShape.h>
#include <TEveLine.h>
#include <TGLViewer.h>
#include <TGLCamera.h>

#include <string>
#include <sstream>

Cube::TFitChangeHandler::TFitChangeHandler() {
    fHitList = new TEveElementList("HitList","Reconstructed 3D Hits");
    fHitList->SetMainColor(kYellow);
    fHitList->SetMainAlpha(1.0);
    gEve->AddElement(fHitList);

    fFitList = new TEveElementList("FitList","Reconstructed Objects");
    fFitList->SetMainColor(kGreen);
    fFitList->SetMainAlpha(0.5);
    gEve->AddElement(fFitList);
    fShowFitsHits = true;
    fShowFitsObjects = true;
}

Cube::TFitChangeHandler::~TFitChangeHandler() {
}

void Cube::TFitChangeHandler::Apply() {
    CUBE_ERROR << "Apply" << std::endl;

    fHitList->DestroyElements();
    fFitList->DestroyElements();
    fCameraWeight = 0.0;

    if (!Cube::TEventDisplay::Get().GUI().GetShowFitsButton()->IsOn()
        && !Cube::TEventDisplay::Get().GUI().GetShowFitHitsButton()->IsOn()) {
        CUBE_ERROR << "Not On" << std::endl;
        return;
    }

    if (Cube::TEventDisplay::Get().GUI().GetShowFitHitsButton()->IsOn()) {
        fShowFitsHits = true;
    }
    else {
        fShowFitsHits = false;
    }

    if (Cube::TEventDisplay::Get().GUI().GetShowFitsButton()->IsOn()) {
        fShowFitsObjects = true;
    }
    else {
        fShowFitsObjects = false;
    }

    if (!Cube::gEvent) {
        CUBE_ERROR << "No event" << std::endl;
        return;
    }

    // Get a TList of all of the selected entries.
    TList selected;
    Cube::TEventDisplay::Get().GUI().GetResultsList()
        ->GetSelectedEntries(&selected);

    int index = 0;

    // Iterate through the list of selected entries.
    TIter next(&selected);
    TGLBEntry* lbEntry;
    while ((lbEntry = (TGLBEntry*) next())) {
        std::string objName(lbEntry->GetTitle());
        CUBE_ERROR << "Get entry " << objName << std::endl;
        Cube::Handle<Cube::ReconObjectContainer> objects
            = gEvent->GetObjectContainer(objName.c_str());
        if (!objects) {
            CUBE_ERROR << "No objects" << std::endl;
            continue;
        }
        index = ShowReconObjects(fFitList,objects,index);
    }

#ifdef SETUP_CAMERA
    if (fCameraWeight > 1
        && Cube::TEventDisplay::Get().GUI()
        .GetRecalculateViewButton()->IsOn()) {
        TGLViewer* glViewer = gEve->GetDefaultGLViewer();
        fCameraCenter *= 1.0/fCameraWeight;
        glViewer->SetDrawCameraCenter(kTRUE);
        glViewer->CurrentCamera().SetExternalCenter(kTRUE);
        glViewer->CurrentCamera().SetCenterVecWarp(fCameraCenter.X(),
                                                   fCameraCenter.Y(),
                                                   fCameraCenter.Z());
    }
#endif

}

int Cube::TFitChangeHandler::ShowReconCluster(
    TEveElementList* list,
    Cube::Handle<Cube::ReconCluster> obj,
    int index,
    bool forceUncertainty) {
    if (!obj) return index;

    Cube::Handle<Cube::ClusterState> state = obj->GetState();
    if (!state) {
        std::cout << "ClusterState missing!" << std::endl;
        return index;
    }

    // Check that the charge per hit is reasonable (this is debugging and
    // should have a "user interface".
    double charge = 0.0;
    double hits = 0.0;
    for (Cube::HitSelection::iterator h = obj->GetHitSelection()->begin();
         h != obj->GetHitSelection()->end(); ++h) {
        charge += (*h)->GetCharge();
        hits += 1.0;
    }
    charge = charge/hits;
    if (charge < 15.0) {
        return index;
    }

    // Increment the index to get a new value for the names.
    ++index;

    if (Cube::TEventDisplay::Get().GUI()
        .GetShowClusterUncertaintyButton()->IsOn()) {
        forceUncertainty = true;
    }

    Cube::TReconClusterElement *eveCluster
        = new Cube::TReconClusterElement(*obj,forceUncertainty);

    list->AddElement(eveCluster);

    if (Cube::TEventDisplay::Get().GUI()
        .GetShowClusterHitsButton()->IsOn()) {
        // Draw the hits.
        Cube::TShowHits showHits;
        Cube::Handle<Cube::HitSelection> fibers = obj->GetHitSelection();
        if (fibers) showHits(eveCluster, *fibers);
    }

    return index;
}

int Cube::TFitChangeHandler::ShowReconShower(
    TEveElementList* list,
    Cube::Handle<Cube::ReconShower> obj,
    int index) {
    if (!obj) return index;

    Cube::Handle<Cube::ShowerState> state = obj->GetState();
    if (!state) {
        std::cout << "ShowerState missing!" << std::endl;
        return index;
    }

    // Get a new index
    ++index;

    Cube::TReconShowerElement *eveShower
        = new Cube::TReconShowerElement(*obj,true);

    list->AddElement(eveShower);

    // Draw the clusters.
    if (Cube::TEventDisplay::Get().GUI()
        .GetShowConstituentClustersButton()->IsOn()) {
        for (Cube::ReconNodeContainer::iterator n = obj->GetNodes().begin();
             n != obj->GetNodes().end(); ++n) {
            index = ShowReconObject(eveShower,(*n)->GetObject(),index, false);
        }
    }

    return index;
}

int Cube::TFitChangeHandler::ShowReconTrack(
    TEveElementList* list,
    Cube::Handle<Cube::ReconTrack> obj,
    int index) {
    if (!obj) return index;
    Cube::Handle<Cube::TrackState> frontState = obj->GetState();
    if (!frontState) {
        std::cout << "TrackState is missing!" << std::endl;
        return index;
    }

    // Get a new index
    ++index;

    std::unique_ptr<TReconTrackElement> eveTrack(
        new TReconTrackElement(
            *obj, true,
            (Cube::TEventDisplay::Get().GUI().GetShowConstituentClustersButton()
             ->IsOn()),
            (Cube::TEventDisplay::Get().GUI().GetShowFitDirectionButton()
             ->IsOn())));
    if (eveTrack->Valid()) list->AddElement(eveTrack.release());

    return index;
}

int Cube::TFitChangeHandler::ShowReconVertex(
    TEveElementList* list,
    Cube::Handle<Cube::ReconVertex> obj,
    int index) {
    if (!obj) return index;

    Cube::Handle<Cube::VertexState> state = obj->GetState();
    if (!state) {
        std::cout << "VertexState missing!" << std::endl;
        return index;
    }

    std::stringstream title;

    ++index;
    TLorentzVector pos = state->GetPosition();
    TLorentzVector var = state->GetPositionVariance();
    double uncertainty = var.X();
    uncertainty = std::max(var.Y(),uncertainty);
    uncertainty = std::max(var.Z(),uncertainty);
    if (uncertainty > 0) uncertainty = std::sqrt(uncertainty);
    uncertainty = std::max(5.0*unit::mm,uncertainty);

    TEveGeoShape *vtxShape = new TEveGeoShape("vertex");
    vtxShape->SetName("vertex");
    vtxShape->SetTitle("A Vertex");
    vtxShape->SetMainColor(kRed);
    // Set the translation
    TGeoTranslation trans(pos.X(),
                          pos.Y(),
                          pos.Z());
    vtxShape->SetTransMatrix(trans);
    TGeoManager* saveGeom = gGeoManager;
    gGeoManager = vtxShape->GetGeoMangeur();
    TGeoShape* geoShape = new TGeoSphere(0.0, uncertainty);
    vtxShape->SetShape(geoShape);
    gGeoManager = saveGeom;
    list->AddElement(vtxShape);

    title << "Vertex(" << obj->GetUniqueID() << "): "
          << unit::AsString(pos.X(),std::sqrt(var.X()),"length")
          <<", "<<unit::AsString(pos.Y(),std::sqrt(var.Y()),"length")
          <<", "<<unit::AsString(pos.Z(),std::sqrt(var.Z()),"length");
    if (obj->GetHitSelection()) {
        title <<", H: "<< obj->GetHitSelection()->size();
    }

    Cube::Handle<Cube::ReconObjectContainer>
        constituents = obj->GetConstituents();
    if (constituents) {
        title << std::endl << "    N: " << constituents->size() << " -- ";
        for (Cube::ReconObjectContainer::iterator o = constituents->begin();
             o != constituents->end(); ++o) {
            Cube::Handle<Cube::ReconTrack> trk = *o;
            if (!trk) continue;
            title << " " << trk->GetUniqueID();
            TVector3 dir = trk->GetDirection();
            TVector3 p1 = pos.Vect() + 15.0*unit::cm*dir;
            TVector3 fr = trk->GetFront()->GetPosition().Vect();
            TVector3 bk = trk->GetBack()->GetPosition().Vect();
            double frDist = (pos.Vect() - fr).Mag();
            double bkDist = (pos.Vect() - bk).Mag();
            double minDist = uncertainty + 2.0*unit::cm;
            if (frDist > minDist && frDist < bkDist) {
                p1 = fr;
            }
            else if (bkDist > minDist && bkDist < frDist) {
                p1 = bk;
            }
            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetTitle("constituent");
            eveHit->SetLineWidth(1);
            eveHit->SetLineColor(kRed);
            eveHit->SetPoint(0,pos.Vect().X(),pos.Vect().Y(),pos.Vect().Z());
            eveHit->SetPoint(1,p1.X(),p1.Y(),p1.Z());
            list->AddElement(eveHit);
        }
    }
    vtxShape->SetTitle(title.str().c_str());

    std::cout << title.str() << std::endl;

    return index;
}

int Cube::TFitChangeHandler::ShowReconObject(TEveElementList* list,
                                           Cube::Handle<Cube::ReconObject> obj,
                                           int index,
                                           bool forceUncertainty) {
    if (!obj) return index;
    // Add this object to the estimated center.
    Cube::Handle<Cube::HitSelection> hits = obj->GetHitSelection();
    if (hits) {
        for (Cube::HitSelection::iterator h = hits->begin();
             h != hits->end(); ++h) {
            fCameraCenter += (*h)->GetCharge()*(*h)->GetPosition();
            fCameraWeight += (*h)->GetCharge();
        }
    }
    Cube::Handle<Cube::ReconVertex> vertex = obj;
    if (vertex) {
        index = ShowReconVertex(list, vertex, index);
        return index;
    }
    if (!fShowFitsObjects) return index;
    Cube::Handle<Cube::ReconCluster> cluster = obj;
    if (cluster) {
        index = ShowReconCluster(
            list, cluster, index, forceUncertainty);
        return index;
    }
    Cube::Handle<Cube::ReconShower> shower = obj;
    if (shower) {
        index = ShowReconShower(list, shower, index);
        return index;
    }
    Cube::Handle<Cube::ReconTrack> track = obj;
    if (track) {
        index = ShowReconTrack(list, track, index);
        return index;
    }
    return index;
}

int Cube::TFitChangeHandler::ShowReconObjects(
    TEveElementList* list,
    Cube::Handle<Cube::ReconObjectContainer> objects,
    int index) {
    if (!objects) return index;
    for (Cube::ReconObjectContainer::reverse_iterator obj = objects->rbegin();
         obj != objects->rend(); ++obj) {
        if (Cube::TEventDisplay::Get().GUI()
            .GetSkipFitVerticesButton()->IsOn()) {
            Cube::Handle<Cube::ReconVertex> vertex = *obj;
            if (vertex) continue;
        }
        if (Cube::TEventDisplay::Get().GUI().GetSkipFitTracksButton()->IsOn()) {
            Cube::Handle<Cube::ReconTrack> track = *obj;
            if (track) continue;
        }
        if (Cube::TEventDisplay::Get().GUI()
            .GetSkipFitClustersButton()->IsOn()) {
            Cube::Handle<Cube::ReconCluster> cluster = *obj;
            if (cluster) continue;
        }
        index = ShowReconObject(list,*obj, index, false);
        if (fShowFitsHits) {
            // Draw the hits.
            Cube::TShowHits showHits;
            Cube::Handle<Cube::HitSelection> fibers = (*obj)->GetHitSelection();
            if (fibers) showHits(fHitList, *fibers);
        }
    }
    return index;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
