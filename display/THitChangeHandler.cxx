#include "THitChangeHandler.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TGUIManager.hxx"
#include "TShowHits.hxx"

#include <CubeLog.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>

#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>

#include <sstream>

Cube::THitChangeHandler::THitChangeHandler() {
    fHitList = new TEveElementList("HitList","Hits in the Event");
    fHitList->SetMainColor(kCyan);
    fHitList->SetMainAlpha(1.0);
    gEve->AddElement(fHitList);
}

Cube::THitChangeHandler::~THitChangeHandler() {
}

void Cube::THitChangeHandler::Apply() {
    // CUBE_ERROR << "About to destroy hit list." << std::endl;

    fHitList->DestroyElements();

    // CUBE_ERROR << "Check to draw the hits." << std::endl;

    if (Cube::TEventDisplay::Get().GUI().GetShowHitsButton()->IsOn()) {
        // CUBE_LOG(0) << "Handle the hits" << std::endl;
        if (!Cube::gEvent) return;

        Cube::Handle<Cube::HitSelection> cubeHits
            = Cube::gEvent->GetHitSelection();

        if (!cubeHits) {
            CUBE_LOG(0) << " No hits in the event" << std::endl;
            return;
        }

        CUBE_LOG(0) << "Showing " << cubeHits->size()
                    <<" hits from "
                    << cubeHits->GetName()
                    << std::endl;

        Cube::TShowHits showCubeHits;
        showCubeHits(fHitList,*cubeHits);
    }

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
