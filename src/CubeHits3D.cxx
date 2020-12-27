#include "CubeHits3D.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"
#include "CubeUnits.hxx"
#include "CubeInfo.hxx"
#include "CubeShareCharge.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>
#include <CubeReconCluster.hxx>

#include <algorithm>
#include <set>
#include <memory>
#include <cmath>

namespace {
    struct hitCompareHitZ {
        bool operator () (const Cube::Handle<Cube::Hit>& lhs,
                          const Cube::Handle<Cube::Hit>& rhs) {
            return lhs->GetPosition().Z() < rhs->GetPosition().Z();
        }
    };

};

Cube::Hits3D::Hits3D()
    : Cube::Algorithm("Cube::Hits3D", "Build 3D Cube Hits from Fibers") {
    fShareCharge = 1;
    fConserveChargeSum = 1;
    fLightSpeed = 200.0*unit::mm/unit::ns;
}

Cube::Hits3D::~Hits3D() { }


#ifdef CUBE_TIME_CHEAT
#include <CubeG4Hit.hxx>
#include <CubeEvent.hxx>
Cube::Event *gOutputEvent = NULL;
namespace {
    std::vector<Cube::Handle<Cube::G4Hit>> CubeSegments(
        Cube::Hits3D::FiberTQ& fiberTQ) {
        std::vector<Cube::Handle<Cube::G4Hit>> result;
        if (fiberTQ.size() != 3) {
            CUBE_ERROR << "Where's the beef" << std::endl;
            throw std::runtime_error("Should be three fibers for CubeSegments");
            return result;
        }
        std::set<Cube::Handle<Cube::G4Hit>> hit0;
        for (int j = 0; j < fiberTQ[0].second->GetContributorCount(); ++j) {
            int seg = fiberTQ[0].second->GetContributor(j);
            Cube::Handle<Cube::G4Hit> g4Hit = gOutputEvent->G4Hits[seg];
            if (!g4Hit) continue;
            hit0.insert(g4Hit);
        }
        std::set<Cube::Handle<Cube::G4Hit>> hit1;
        for (int j = 0; j < fiberTQ[1].second->GetContributorCount(); ++j) {
            int seg = fiberTQ[1].second->GetContributor(j);
            Cube::Handle<Cube::G4Hit> g4Hit = gOutputEvent->G4Hits[seg];
            if (!g4Hit) continue;
            hit1.insert(g4Hit);
        }
        std::set<Cube::Handle<Cube::G4Hit>> hit2;
        for (int j = 0; j < fiberTQ[2].second->GetContributorCount(); ++j) {
            int seg = fiberTQ[2].second->GetContributor(j);
            Cube::Handle<Cube::G4Hit> g4Hit = gOutputEvent->G4Hits[seg];
            if (!g4Hit) continue;
            hit2.insert(g4Hit);
        }
        for(std::set<Cube::Handle<Cube::G4Hit>>::iterator s = hit0.begin();
            s != hit0.end(); ++s) {
            if (hit1.find(*s) == hit1.end()) continue;
            if (hit2.find(*s) == hit2.end()) continue;
            result.push_back(*s);
        }
        return result;
    }
}
#endif

// Sort the fiber times, and find the average time for hits in a narrow window
// of time.  This works so that if one fiber is "out of time" for this hit
// (usually because the first photon came from a different cube), we still get
// the average for the fibers that are in coincidence.  If none of the fibers
// are in coincidence, this returns time of the last fiber.
std::pair<double,double> Cube::Hits3D::HitTime(FiberTQ& fiberTQ) const {
    // Set the window to average hit times over.
    const double timeWindow = 2.5*unit::ns;

    // Intrinsic resolution for a hit coming from the electronics.
    const double hitRes = 0.7*unit::ns;

    // Order the fiber times.  The times are transit time corrected.
    std::sort(fiberTQ.begin(),fiberTQ.end());

    // Now take the charge weighted average and RMS of the late hits.  The
    // average is also calculated for all hits to help calculate the
    // uncertainty.  Using the late hits is an approximation of using the
    // fibers that have a localized region of energy deposition.  Because of
    // geometry, the latest will usually be from fibers where the current cube
    // is the close to the sensor.  Lot's of different methods were tried, and
    // the latest hit seems to be the best one.
    double tSum = 0.0; // The sum of times after truncation
    double ttSum = 0.0; // The sum of time squared after truncation
    double qSum = 0.0;  // The sum of charge after truncation.
    double closeFiber = 0.0; // The number of fibers "causally connected".
    double lastTime = fiberTQ.back().first;
    for (FiberTQ::iterator t = fiberTQ.begin(); t != fiberTQ.end(); ++t) {
        double qq = t->second->GetCharge();
        double tt = t->first;
        if (tt < lastTime - timeWindow) continue;
        tSum += qq*tt;
        ttSum += qq*tt*tt;
        qSum += qq;
        closeFiber += 1.0;
    }

    if (qSum <= 0.0) {
        CUBE_ERROR << "Invalid hit time" << std::endl;
        throw std::runtime_error("Invalid cube hit");
        return std::pair<double,double>(-1.0,-1.0);  // Keep the compiler happy.
    }

    // The RMS for the causal hits.
    double tHit = tSum/qSum;
    ttSum = ttSum/qSum;
    double tRMS = ttSum - tHit*tHit;
    if (tRMS < 0.0) tRMS = 0.0;
    tRMS = std::sqrt(tRMS);

#define ADD_ELECTRONICS_HIT_RESOLUTION
#ifdef ADD_ELECTRONICS_HIT_RESOLUTION
    // The variance of the electronics.
    tRMS = tRMS*tRMS + hitRes*hitRes;
    tRMS = std::sqrt(tRMS);
#endif

#ifdef CUBE_TIME_CHEAT
    #warn THE CUBE TIME IS CHEATING.  DID YOU REALLY MEAN THIS??
    std::vector<Cube::Handle<Cube::G4Hit>> segs = CubeSegments(fiberTQ);
    if (segs.empty()) return std::pair<double,double>(-1.0,-1.0);
    double tCheat = segs.front()->GetStart().T();
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator s = segs.begin();
         s != segs.end(); ++s) {
        tCheat = std::min(tCheat,(*s)->GetStart().T());
    }
    tCheat += 100.0;
    std::cout << "Use time cheat with "
              << tHit << " " << tCheat << " " << tRMS << std::endl;
    return std::pair<double,double>(tCheat,tRMS);
#endif

    if (closeFiber > 1.0) return std::pair<double,double>(tHit,tRMS);

#define ADD_FIBER_HIT_RESOLUTION
#ifdef ADD_FIBER_HIT_RESOLUTION
    // Set the timing resolution for a single fiber.  This covers the full
    // range of time that the hit might have come from when we only have the
    // fiber.  This empirically based on the width of the non-causal hit
    // residual distribution, but works out to be about twice the fiber length
    // over sqrt(12).
    const double fiberRes = 6.0*unit::ns;

    // This is a hit getting it's time from a single fiber, so increase the
    // uncertainty.
    tRMS = tRMS*tRMS + fiberRes*fiberRes;
    tRMS = std::sqrt(tRMS);
#endif

#ifdef USE_FIBER_SPREAD
    // The variance due to the possible range of times the energy hit the
    // fiber.  This uses the spread of times from the hits on the fiber, but
    // allows out of time hits to contribute to the RMS.
    double tAll = 0.0; // The sum of all times before truncation
    double ttAll = 0.0; // The sum of all time squared before truncation
    double qAll = 0.0; // The sum of all charge before truncation
    for (FiberTQ::iterator t = fiberTQ.begin(); t != fiberTQ.end(); ++t) {
        double qq = t->second->GetCharge();
        double tt = t->first;
        double ww = std::exp(0.5*(tt-lastTime)/timeWindow);
        tAll += ww*tt;
        ttAll += ww*tt*tt;
        qAll += ww;
    }

    double aHit = tAll/qAll;
    ttAll = ttAll/qAll;
    double aRMS = ttAll - aHit*aHit;
    if (aRMS < 0.0) aRMS = 0.0;
    aRMS = std::sqrt(aRMS);
    tRMS = tRMS*tRMS + aRMS*aRMS;
    tRMS = std::sqrt(tRMS);
#endif

    return std::pair<double,double>(tHit,tRMS);

}

bool Cube::Hits3D::MakeHit(Cube::HitSelection& writableHits,
                           const Cube::Handle<Cube::Hit>& hit1,
                           const Cube::Handle<Cube::Hit>& hit2,
                           const Cube::Handle<Cube::Hit>& hit3) const {

    // Figure out the identifier for the new hit.
    int cubeId = -1;
    if (hit3) {
        cubeId = Cube::Info::Combine3DST(hit1->GetIdentifier(),
                                         hit2->GetIdentifier(),
                                         hit3->GetIdentifier());
    }
    else {
        cubeId = Cube::Info::Combine3DST(hit1->GetIdentifier(),
                                         hit2->GetIdentifier(),
                                         -1);
    }

    // Figure out the position for the new hit.
    TVector3 hitPos = hit1->GetPosition();
    TVector3 hitSize = hit1->GetSize();
    if (Cube::Info::CubeNumber(hit1->GetIdentifier()) < 0) {
        hitPos.SetX(hit2->GetPosition().X());
        hitSize.SetX(hit2->GetSize().X());
    }
    else if (Cube::Info::CubeBar(hit1->GetIdentifier()) < 0) {
        hitPos.SetY(hit2->GetPosition().Y());
        hitSize.SetY(hit2->GetSize().Y());
    }
    else if (Cube::Info::CubePlane(hit1->GetIdentifier()) < 0) {
        hitPos.SetZ(hit2->GetPosition().Z());
        hitSize.SetZ(hit2->GetSize().Z());
    }

    // Set the charge to the total charge.
    double qHit = hit1->GetCharge();
    qHit += hit2->GetCharge();
    if (hit3) qHit += hit3->GetCharge();

    FiberTQ fiberTQ;

    // Get the times for the hit.
    double dHit1 = (hitPos - hit1->GetPosition()).Mag();
    double dHit2 = (hitPos - hit2->GetPosition()).Mag();

    double tHit1 = hit1->GetTime() - dHit1/fLightSpeed;
    double tHit2 = hit2->GetTime() - dHit2/fLightSpeed;

    fiberTQ.push_back(std::make_pair(tHit1,hit1));
    fiberTQ.push_back(std::make_pair(tHit2,hit2));

    double dtHit = std::abs(tHit1-tHit2);
    if (hit3) {
        double dHit3 = (hitPos - hit3->GetPosition()).Mag();
        double tHit3 = hit3->GetTime() - dHit3/fLightSpeed;
        dtHit = std::max(dtHit,std::abs(tHit1-tHit3));
        dtHit = std::max(dtHit,std::abs(tHit2-tHit3));
        fiberTQ.push_back(std::make_pair(tHit3,hit3));
    }

    // Check that the hits are "at the same time".  This is the largest
    // possible different between simultaneous cube times.  It's determined by
    // geometry.
    if (dtHit > 20*unit::ns) return false;

    // Protect against future special cases.
    if (fiberTQ.empty()) return false;

    // Create the new writable hit.
    Cube::Handle<Cube::WritableHit> hit(new Cube::WritableHit);
    hit->SetIdentifier(cubeId);
    hit->SetPosition(hitPos);
    hit->SetSize(hitSize);

    hit->SetCharge(qHit);
    hit->SetChargeUncertainty(std::sqrt(qHit));

    std::pair<double,double> timePair = HitTime(fiberTQ);
    if (timePair.second <= 0.0) {
        // Something went wrong while building the hit.
        return false;
    }
    hit->SetTime(timePair.first);
    hit->SetTimeUncertainty(timePair.second);

    // Set the position uncertainty with charge weighting so that average
    // positions are charge weighted.  This is needed so that ReconClusters
    // are correctly genererated.  The actual uncertainty for a single hit is
    // still the size of a cube.
    double dx = std::sqrt(1.0*unit::cm*unit::cm/12.0/qHit);
    TVector3 unc(dx,dx,dx);
    hit->SetUncertainty(unc);

    // Add the first hit.  There MUST be a first hit!
    hit->AddHit(hit1);
    for (int i=0; i<hit1->GetContributorCount(); ++i) {
        hit->AddContributor(hit1->GetContributor(i));
    }

    // Add the second hit.  There MUST be a second hit!
    hit->AddHit(hit2);
    for (int i=0; i<hit2->GetContributorCount(); ++i) {
        hit->AddContributor(hit2->GetContributor(i));
    }

    if (hit3) {
        // Add the third hit.  There SHOULD be a third hit, but it's not
        // required.
        hit->AddHit(hit3);
        for (int i=0; i<hit3->GetContributorCount(); ++i) {
            hit->AddContributor(hit3->GetContributor(i));
        }
    }

    writableHits.push_back(hit);

    return true;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::Hits3D::Process(const Cube::AlgorithmResult& fibers,
                      const Cube::AlgorithmResult&,
                      const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Hits3D::Process" << std::endl;

    Cube::Handle<Cube::HitSelection> fiberHits = fibers.GetHitSelection();
    if (!fiberHits) {
        CUBE_ERROR << "No input hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    typedef std::set< Cube::Handle<Cube::Hit> >  HitSet;
    HitSet usedSet;
    HitSet unusedSet;

    Cube::HitSelection xzHits;
    Cube::HitSelection yzHits;
    std::map<int, Cube::HitSelection> xyHits;
    for (Cube::HitSelection::iterator h = fiberHits->begin();
         h != fiberHits->end(); ++h) {
        int plane = Cube::Info::IdentifierProjection((*h)->GetIdentifier());
        switch (plane) {
        case Cube::Info::kXZProj:
            xzHits.push_back(*h);
            break;
        case Cube::Info::kYZProj:
            yzHits.push_back(*h);
            break;
        case Cube::Info::kXYProj:
            xyHits[(*h)->GetIdentifier()].push_back(*h);
            break;
        default:
            CUBE_ERROR << "Invalid fiber plane" << std::endl;
        }
    }

    std::sort(xzHits.begin(), xzHits.end(), hitCompareHitZ());
    std::sort(yzHits.begin(), yzHits.end(), hitCompareHitZ());

    CUBE_LOG(0) << "XZ Hits: " << xzHits.size()
                << " YZ Hits: " << yzHits.size()
                << " XY Hits: " << xyHits.size()
                << std::endl;

    Cube::HitSelection writableHits;
    for (Cube::HitSelection::iterator xz=xzHits.begin();
         xz!=xzHits.end(); ++xz) {
        Cube::HitSelection::iterator yzBegin = yzHits.begin();
        Cube::HitSelection::iterator yzEnd = yzHits.end();

        for (Cube::HitSelection::iterator yz=yzBegin;
             yz!=yzEnd; ++yz) {
            if ((*yz)->GetPosition().Z()<(*xz)->GetPosition().Z()-5) continue;
            if ((*yz)->GetPosition().Z()>(*xz)->GetPosition().Z()+5) continue;
#ifdef MAKE_CONFUSED_2D_HITS
            if (MakeHit(writableHits,*xz,*yz,Cube::Handle<Cube::Hit>())) {
                usedSet.insert(*xz);
                usedSet.insert(*yz);
            }
            continue;
#endif
            int xyFiber = Cube::Info::Identifier3DST(
                Cube::Info::CubeNumber((*xz)->GetIdentifier()),
                Cube::Info::CubeBar((*yz)->GetIdentifier()),
                -1);
            std::map<int, Cube::HitSelection>::iterator xySel
                = xyHits.find(xyFiber);
            if (xySel == xyHits.end()) continue;
            Cube::HitSelection& xy = xySel->second;
            for (Cube::HitSelection::iterator h = xy.begin();
                 h != xy.end(); ++h) {
                if (MakeHit(writableHits,*xz,*yz,*h)) {
                    usedSet.insert(*xz);
                    usedSet.insert(*yz);
                    usedSet.insert(*h);
                }
            }
        }
    }
    std::cout << "Hits Generated " << writableHits.size() << std::endl;

    // Clear the selections to make sure the handles reset.
    xzHits.clear();
    yzHits.clear();
    xyHits.clear();

    // All of the 3 fiber 3D hits have been found, now check any unassociated
    // hits.  Unassociated hits are going into the unused list.
    for (Cube::HitSelection::iterator h = fiberHits->begin();
         h != fiberHits->end(); ++h) {
        if ((*h)->GetConstituentCount() > 0) {
            continue;
        }
        if (usedSet.find(*h) != usedSet.end()) continue;
        unusedSet.insert(*h);
    }

    CUBE_LOG(0) << "Total 3D Hits: " << writableHits.size() << std::endl;

    if (writableHits.size() < 1) {
        CUBE_ERROR << "No 3D hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

#ifndef MAKE_CONFUSED_2D_HITS
    // Share the charge among the 3D hits so that the total charge in the
    // event is not overcounted.  This corrects for the attenuation in the
    // fiber.
    Cube::ShareCharge shareCharge;
    shareCharge.SetChargeConservation(fConserveChargeSum);
    if (fShareCharge == 1) {
        // Apply (almost) the same formalism as for the MaximumEntropy
        // version.  Share the charge by predicting the measurement in each
        // fiber based on the deposit in each cube.  The deposit in a cube is
        // constrained to be positive.  An approximation of the maximum
        // entropy prior is applied to break degeneracies.  The prior prefers
        // that cubes have the same charge (it's a very weak prior).
        shareCharge.OptimizeCubes(writableHits);
    }
    else if (fShareCharge == 2) {
        // Share the charge assuming that all of the fibers in a cube will
        // measure the same amount of deposit, constrained by the charge
        // present at the fiber being equal to the measured charge.  The
        // constraint on the charge at the fiber is applied using a Lagrange
        // multiplier.  The constraint at the cube is the likelihood that all
        // three measurements came from the same mean.
        shareCharge.ApplyConstraints(writableHits);
    }
#ifdef ROOT_CHANGED_API_BREAKS_THIS
    else if (fShareCharge == 3) {
        // Share the charge applying a Bayesian probability with a maximum
        // entropy prior.  The probability is based on predicting the
        // measurement in each fiber based on the deposit in each cube.  The
        // deposit in a cube is constrained to be positive.  The maximum
        // entropy prior is that all cubes should have the same charge (it's a
        // very weak prior).  This is impossibly slow for large events and has
        // precision problems.  The ApplyConstraints and OptimizeCubes
        // versions are better.
        shareCharge.MaximizeEntropy(writableHits);
    }
#endif
#endif

    // Copy the writable hits into the clustered hit selection;
    Cube::HitSelection clustered;
    for (Cube::HitSelection::iterator h = writableHits.begin();
         h != writableHits.end(); ++h) {
        Cube::Handle<Cube::WritableHit> hit = *h;
        Cube::Handle<Cube::Hit> newHit(new Cube::Hit(*hit));
        clustered.push_back(newHit);
    }

    // Build an object container with the hits clustered into a convenient
    // form. This is probably mostly used for display.
    Cube::Handle<Cube::ReconObjectContainer> finalObjects(
        new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    // Create a single cluster from the new 3D hits.
    Cube::Handle<Cube::ReconCluster> cluster3D
        = Cube::CreateCluster("hits3D", clustered.begin(),clustered.end());
    if (cluster3D) {
        finalObjects->push_back(cluster3D);
    }

    // Collect all of the hits into the used and unused hit selections.
    Cube::MakeUsed makeUsed(*fiberHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
