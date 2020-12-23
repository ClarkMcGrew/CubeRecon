#include "TEventChangeManager.hxx"
#include "TVEventChangeHandler.hxx"
#include "TGUIManager.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TRecurseGeoManager.hxx"

#include "CubeAlgorithmResult.hxx"

#include <TQObject.h>
#include <TGButton.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TEveGeoShape.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>

#include <iostream>

ClassImp(Cube::TEventChangeManager)

namespace {

    /// This takes a geometry id and "clones" it into the Eve display.
    TEveGeoShape* GeometryClone(TGeoVolume* currVolume) {

        TGeoShape* currShape = currVolume->GetShape();
        TEveGeoShape *fakeShape = new TEveGeoShape(currShape->GetName());

        TGeoMatrix* currMat = gGeoManager->GetCurrentMatrix();
        TGeoMatrix* mat = currMat->MakeClone();
        fakeShape->SetTransMatrix(*mat);

        // Clone the shape so that it can be displayed.  This has to play some
        // fancy footsie to get the gGeoManager memory management right.  It
        // first saves the current manager, then gets an internal geometry
        // manager used by TEveGeoShape, and then resets the old manager once
        // the shape is created.
        TGeoManager* saveGeom = gGeoManager;
        gGeoManager = fakeShape->GetGeoMangeur();
        TGeoShape* clonedShape
            = dynamic_cast<TGeoShape*> (currShape->Clone("fakeShape"));
        fakeShape->SetShape(clonedShape);
        gGeoManager = saveGeom;

        return fakeShape;
    }

    struct PickVolumes: TRecurseGeoManager {
        std::vector<std::string> fVolumesToShow;
        TEveElementList* fEveList;
        PickVolumes(TEveElementList* eveList) : fEveList(eveList) {}
        bool Action(int) {
            if (!fEveList) return true;
            TGeoNode* node = gGeoManager->GetCurrentNode();
            std::string name(node->GetName());
            for (std::vector<std::string>::iterator v = fVolumesToShow.begin();
                 v != fVolumesToShow.end(); ++v) {
                if (name.find(*v) == std::string::npos) continue;
                std::cout << "Show volume " << node->GetName() << std::endl;
                TEveGeoShape *shape = GeometryClone(node->GetVolume());
                shape->SetMainColor(node->GetVolume()->GetLineColor());
                shape->SetMainTransparency(90);
                fEveList->AddElement(shape);
                break;
            }
            return true;
        }
    };
}

TTree* Cube::gCubeReconTree = NULL;
Cube::Event* Cube::gEvent = NULL;
Cube::AlgorithmResult* Cube::gAlgorithm = NULL;
int Cube::gCubeReconEntryNumber = 0;

Cube::TEventChangeManager::TEventChangeManager() {
    TGButton* button = Cube::TEventDisplay::Get().GUI().GetNextEventButton();
    if (button) {
        button->Connect("Clicked()",
                        "Cube::TEventChangeManager",
                        this,
                        "ChangeEvent(=1)");
    }

    button = Cube::TEventDisplay::Get().GUI().GetDrawEventButton();
    if (button) {
        button->Connect("Clicked()",
                        "Cube::TEventChangeManager",
                        this,
                        "ChangeEvent(=0)");
    }

    button = Cube::TEventDisplay::Get().GUI().GetPrevEventButton();
    if (button) {
        button->Connect("Clicked()",
                        "Cube::TEventChangeManager",
                        this,
                        "ChangeEvent(=-1)");
    }

    ClearVolumesToShow();
}

Cube::TEventChangeManager::~TEventChangeManager() {}

void Cube::TEventChangeManager::SetEventSource(TFile* source) {
    if (!source) {
        std::cout << "Invalid event source" << std::endl;
        return;
    }
    fEventSource = source;
    Cube::gCubeReconTree = (TTree*) fEventSource->Get("CubeEvents");
    if (!gCubeReconTree) {
        std::cout << "Missing the event tree" << std::endl;
        return;
    }

    Cube::gCubeReconTree->SetBranchAddress("Event",&Cube::gEvent);
    Cube::gCubeReconEntryNumber = 0;

    gCubeReconTree->GetEntry(Cube::gCubeReconEntryNumber);
    Cube::gAlgorithm = Cube::gEvent;

    fEventSource->Get("CubeReconGeometry");

    if (GetShowGeometry() && !fVolumesToShow.empty()) {
        TEveElementList* simple = new TEveElementList("simplifiedGeometry");
        PickVolumes visitor(simple);
        visitor.fVolumesToShow.clear();
        visitor.fVolumesToShow.reserve(fVolumesToShow.size());
        std::copy(fVolumesToShow.begin(), fVolumesToShow.end(),
                  std::back_inserter(visitor.fVolumesToShow));
        visitor.Apply();
        gEve->AddGlobalElement(simple);
    }
    else {
        std::cout << "Geometry not shown" << std::endl;
    }

    NewEvent();
    UpdateEvent();
}

void Cube::TEventChangeManager::AddNewEventHandler(
    Cube::TVEventChangeHandler* handler) {
    fNewEventHandlers.push_back(handler);
}

void Cube::TEventChangeManager::AddUpdateHandler(
    Cube::TVEventChangeHandler* handler) {
    fUpdateHandlers.push_back(handler);
}

void Cube::TEventChangeManager::ChangeEvent(int change) {
    std::cout << "Change Event by " << change << " entries" << std::endl;
    if (!GetEventSource()) {
        std::cout << "Event source is not available" << std::endl;
        UpdateEvent();
        return;
    }

    int oldEntry = gCubeReconEntryNumber;

    Cube::gCubeReconEntryNumber += change;
    if (Cube::gCubeReconEntryNumber >= Cube::gCubeReconTree->GetEntries()) {
        std::cout << "Cannot move past last event" << std::endl;
        Cube::gCubeReconEntryNumber = Cube::gCubeReconTree->GetEntries() - 1;
    }
    if (Cube::gCubeReconEntryNumber < 0) {
        std::cout << "Cannot move before first event" << std::endl;
        Cube::gCubeReconEntryNumber = 0;
    }
    change = Cube::gCubeReconEntryNumber - oldEntry;

    if (change != 0) NewEvent();

    UpdateEvent();
}

void Cube::TEventChangeManager::NewEvent() {
    Cube::gCubeReconTree->GetEntry(Cube::gCubeReconEntryNumber);
    Cube::gAlgorithm = Cube::gEvent;

    // Run through all of the handlers.
    for (Handlers::iterator h = fNewEventHandlers.begin();
         h != fNewEventHandlers.end(); ++h) {
        (*h)->Apply();
    }

}

void Cube::TEventChangeManager::UpdateEvent() {

    if (!Cube::gEvent) {
        std::cout << "Invalid event" << std::endl;;
        return;
    }

    std::cout << "DISPLAY EVENT NUMBER " << gEvent->GetRunId()
              << "/" << gEvent->GetEventId()
              << std::endl;

    // Run through all of the handlers.
    for (Handlers::iterator h = fUpdateHandlers.begin();
         h != fUpdateHandlers.end(); ++h) {
        (*h)->Apply();
    }

    std::cout << "EVENT NUMBER " << gEvent->GetRunId()
              << "/" << gEvent->GetEventId()
              << std::endl;

    // Make sure EVE is up to date.
    gEve->Redraw3D(kTRUE);
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
