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
    std::vector<double> qECal;
    double minECal = 1E+20;
    double maxECal = 1.0;
    double aveECal = 0.0;
    double nECal = 0.0;
    std::vector<double> q3DST;
    double min3DST = 1E+20;
    double max3DST = 1.0;
    double ave3DST = 0.0;
    double n3DST = 0.0;
    std::vector<double> qTPC;
    double minTPC = 1E+20;
    double maxTPC = 1.0;
    double aveTPC = 0.0;
    double nTPC = 0.0;
    for (std::map<Int_t,Cube::Handle<Cube::Hit>>::iterator id
             = firstHits.begin();
         id != firstHits.end(); ++id) {
        Cube::Handle<Cube::Hit> hit = id->second;
        double charge = hit->GetCharge();
        if (Cube::Info::IsECal(hit->GetIdentifier())) {
            nECal += 1.0;
            aveECal += charge;
            minECal = std::min(minECal, charge);
            maxECal = std::max(maxECal, charge);
            qECal.push_back(charge);
        }
        if (Cube::Info::Is3DST(hit->GetIdentifier())) {
            if (charge < fThreshold) continue;
            n3DST += 1.0;
            ave3DST += charge;
            min3DST = std::min(min3DST, charge);
            max3DST = std::max(max3DST, charge);
            q3DST.push_back(charge);
        }
        if (Cube::Info::IsTPC(hit->GetIdentifier())) {
            nTPC += 1.0;
            aveTPC += charge;
            minTPC = std::min(minTPC, charge);
            maxTPC = std::max(maxTPC, charge);
            qTPC.push_back(charge);
        }
    }
    if (nECal > 0.0) {
        aveECal /= nECal;
        std::sort(qECal.begin(),qECal.end());
        maxECal = qECal[0.95*qECal.size()];
    }
    if (n3DST > 0.0) {
        ave3DST /= n3DST;
        if (q3DST.size() > 20) {
            std::sort(q3DST.begin(),q3DST.end());
            min3DST = q3DST[0.50*q3DST.size()];
            max3DST = q3DST[0.95*q3DST.size()];
        }
    }
    if (nTPC > 0.0) {
        aveTPC /= nTPC;
        std::sort(qTPC.begin(),qTPC.end());
        maxTPC = qTPC[0.95*qTPC.size()];
    }

    CUBE_LOG(1) << "3DST Charge Range: " << min3DST << " to " << max3DST
                << "   Average: " << ave3DST
                << " with threshold of " << fThreshold << std::endl;
    CUBE_LOG(1) << "TPC Charge Range: " << minTPC << " to " << maxTPC
                << "   Average: " << aveTPC
                << std::endl;
    CUBE_LOG(1) << "ECal Charge Range: " << minECal << " to " << maxECal
                << "   Average: " << aveECal
                << std::endl;

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

        if (hit->GetConstituentCount() < 1) {
            if (!Cube::TEventDisplay::Get().GUI().GetShowSimpleHitsButton()
                ->IsOn()) continue;
        }
        else {
            if (!Cube::TEventDisplay::Get().GUI().GetShowCompHitsButton()
                ->IsOn()) continue;
        }

        hitPos = hit->GetPosition();
        hitSize = hit->GetSize();

        std::ostringstream title;
        title << "Hit " << unit::AsString(hit->GetTime(),"time")
              << " " << unit::AsString(hit->GetCharge(), "pe");
        if (hit->GetConstituentCount() < 1
            && Cube::Info::Is3DST(hit->GetIdentifier())) {
            hitPos = hitPos + hitSize;
        }
        else {
            title << " composite: " << hit->GetConstituentCount();
        }

        ++hitCount;

        int color = 0;
        if (Cube::Info::Is3DST(hit->GetIdentifier())) {
            color = TEventDisplay::Get().LogColor(charge,
                                                  min3DST,max3DST,2.0);
        }
        else if (Cube::Info::IsECal(hit->GetIdentifier())) {
            color = TEventDisplay::Get().LogColor(charge,
                                                  minECal,maxECal,2.0);
        }
        else {
            color = TEventDisplay::Get().LogColor(charge,
                                                  minTPC,maxTPC,2.0);
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
