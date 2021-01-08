#include "TEventDisplay.hxx"
#include "TGUIManager.hxx"
#include "TShowHits.hxx"

#include <TUnitsTable.hxx>

#include <CubeHitSelection.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>
#include <TEveManager.h>
#include <TEveLine.h>
#include <TEveBox.h>

#include <map>
#include <sstream>

Cube::TShowHits::TShowHits()
    : fThreshold(10.0) { }

bool Cube::TShowHits::operator () (TEveElementList* elements,
                                   const Cube::HitSelection& hits) {

    std::map<Int_t, Cube::Handle<Cube::Hit> > firstHits;

    fThreshold =
        Cube::TEventDisplay::Get().GUI().GetHitChargeThreshold()->GetNumber();


    // Find the first hit for each cube, and get the total charge in the cube.
    for (Cube::HitSelection::const_iterator h = hits.begin();
         h != hits.end(); ++h) {
        if ((*h)->GetConstituentCount() < 1) {
            if (!Cube::TEventDisplay::Get().GUI().GetShowSimpleHitsButton()
                ->IsOn()) continue;
        }
        else {
            if (!Cube::TEventDisplay::Get().GUI().GetShowCompHitsButton()
                ->IsOn()) continue;
        }

        Cube::Handle<Cube::Hit> firstHit = firstHits[(*h)->GetIdentifier()];
        if (!firstHit || ((*h)->GetTime() < firstHit->GetTime())) {
            firstHits[(*h)->GetIdentifier()] = *h;
        }
    }

    // Find the range of charges that will be drawn.  Also the average.
    double minCharge = 1E+20;
    double maxCharge = 1.0;
    double aveCharge = 0.0;
    double nCharge = 0.0;
    for (std::map<Int_t,Cube::Handle<Cube::Hit>>::iterator id
             = firstHits.begin();
         id != firstHits.end(); ++id) {
        Cube::Handle<Cube::Hit> hit = id->second;
        if (hit->GetCharge() < fThreshold) continue;
        double charge = hit->GetCharge();
        nCharge += 1.0;
        aveCharge += charge;
        minCharge = std::min(minCharge, charge);
        maxCharge = std::max(maxCharge, charge);
    }
    if (nCharge > 0.0) aveCharge /= nCharge;

    CUBE_LOG(1) << "Charge Range: " << minCharge << " to " << maxCharge
                << "   Average: " << aveCharge
                << "   With threshold of " << fThreshold << std::endl;

    int hitCount = 0;

    // Now plot all of the hits.
    TVector3 hitPos;
    TVector3 hitSize;
    for (Cube::HitSelection::const_iterator h = hits.begin();
         h != hits.end(); ++h) {
        Cube::Handle<Cube::Hit> hit = *h;
        double charge = hit->GetCharge();
        if (charge < fThreshold) continue;
        if (charge < 1.0) charge = 1.0;
        bool onlyFirstHits = false;
        if (onlyFirstHits) {
            std::map<Int_t,Cube::Handle<Cube::Hit>>::iterator fh
                = firstHits.find(hit->GetIdentifier());
            if (fh == firstHits.end()) continue;
            if (fh->second == *h) {
                std::cout << "First Hit" << std::endl;
            }
            else {
                continue;
            }
        }

        hitPos = hit->GetPosition();
        hitSize = hit->GetSize();

        double chargeScale = 1.0;
        std::ostringstream title;
        title << "Hit " << unit::AsString(hit->GetTime(),"time")
              << " " << unit::AsString(hit->GetCharge(), "pe");
        if (hit->GetConstituentCount() < 1
            && Cube::Info::Is3DST(hit->GetIdentifier())) {
            chargeScale = 3.0;
            hitPos = hitPos + hitSize;
        }
        else {
            title << " composite: " << hit->GetConstituentCount();
        }

        ++hitCount;

        int color = 0;
        if (Cube::Info::Is3DST(hit->GetIdentifier())) {
            color = TEventDisplay::Get().LogColor(chargeScale*charge,
                                                  10.0,400.0,2.0);
        }
        else {
            color = TEventDisplay::Get().LogColor(charge,
                                                  10.0,maxCharge,2.0);
        }

        double hitScale = 0.8;
        TEveBox* eveHit = new TEveBox("Hit");
        eveHit->SetTitle(title.str().c_str());
        eveHit->SetVertex(0,
                           hitPos.X()-hitScale*hitSize.X(),
                           hitPos.Y()-hitScale*hitSize.Y(),
                           hitPos.Z()-hitScale*hitSize.Z());
        eveHit->SetVertex(1,
                           hitPos.X()-hitScale*hitSize.X(),
                           hitPos.Y()+hitScale*hitSize.Y(),
                           hitPos.Z()-hitScale*hitSize.Z());
        eveHit->SetVertex(2,
                           hitPos.X()+hitScale*hitSize.X(),
                           hitPos.Y()+hitScale*hitSize.Y(),
                           hitPos.Z()-hitScale*hitSize.Z());
        eveHit->SetVertex(3,
                           hitPos.X()+hitScale*hitSize.X(),
                           hitPos.Y()-hitScale*hitSize.Y(),
                           hitPos.Z()-hitScale*hitSize.Z());
        eveHit->SetVertex(4,
                           hitPos.X()-hitScale*hitSize.X(),
                           hitPos.Y()-hitScale*hitSize.Y(),
                           hitPos.Z()+hitScale*hitSize.Z());
        eveHit->SetVertex(5,
                           hitPos.X()-hitScale*hitSize.X(),
                           hitPos.Y()+hitScale*hitSize.Y(),
                           hitPos.Z()+hitScale*hitSize.Z());
        eveHit->SetVertex(6,
                           hitPos.X()+hitScale*hitSize.X(),
                           hitPos.Y()+hitScale*hitSize.Y(),
                           hitPos.Z()+hitScale*hitSize.Z());
        eveHit->SetVertex(7,
                           hitPos.X()+hitScale*hitSize.X(),
                           hitPos.Y()-hitScale*hitSize.Y(),
                           hitPos.Z()+hitScale*hitSize.Z());
        eveHit->SetMainColor(color);

        elements->AddElement(eveHit);
    }

    CUBE_LOG(1) << "Show " << hitCount << " hits" << std::endl;

    return true;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
