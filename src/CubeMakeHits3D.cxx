#include "CubeMakeHits3D.hxx"
#include "CubeTimeSlice.hxx"
#include "CubeHits3D.hxx"
#include "CubeMakeUsed.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconCluster.hxx>
#include <CubeUnits.hxx>

#include <CubeInfo.hxx>

#include <sstream>
#include <iomanip>
#include <memory>

Cube::MakeHits3D::MakeHits3D()
    : Cube::Algorithm("MakeHits3D","Build 2D hits into 3D hits") {
}

Cube::MakeHits3D::~MakeHits3D() {
}

Cube::Handle<Cube::AlgorithmResult>
Cube::MakeHits3D::Process(const Cube::AlgorithmResult& input,
                          const Cube::AlgorithmResult&,
                          const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "MakeHits3D::Process" << std::endl;

    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        CUBE_ERROR << "No hits" << std::endl;
        return result;
    }
    if (inputHits->empty()) {
        CUBE_ERROR << "Hits empty" << std::endl;
        return result;
    }

    Cube::HitSelection onlyECal;
    Cube::HitSelection only3DST;
    Cube::HitSelection onlyTPC;
    for (Cube::HitSelection::iterator h = inputHits->begin();
         h != inputHits->end(); ++h) {
        if (Cube::Info::Is3DST((*h)->GetIdentifier())) {
            only3DST.push_back((*h));
            continue;
        }
        if (Cube::Info::IsTPC((*h)->GetIdentifier())) {
            onlyTPC.push_back((*h));
            continue;
        }
        if (Cube::Info::IsECal((*h)->GetIdentifier())) {
            if (Cube::Info::ECalModule((*h)->GetIdentifier()) > 29) continue;
            onlyECal.push_back((*h));
            continue;
        }
    }

    CUBE_LOG(0) << "Read 3DST: " << only3DST.size()
                << ", TPC: " << onlyTPC.size()
                << ", ECal: " << onlyECal.size()
                << std::endl;

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer> finalObjects(
        new Cube::ReconObjectContainer("final"));

    // Save the final objects last.
    result->AddObjectContainer(finalObjects);

    // Add a hit selection for the unused hits.  This needs to happen before
    // the used hit selection is added.
    Cube::Handle<Cube::HitSelection> unusedHits(
        new Cube::HitSelection("unused"));
    result->AddHitSelection(unusedHits);

    // Create the hit selection for all of the used hits.
    Cube::Handle<Cube::HitSelection> usedHits(
        new Cube::HitSelection("used"));
    result->AddHitSelection(usedHits);

    if (!onlyTPC.empty()) {
        int problems = 0;
        Cube::Handle<Cube::ReconCluster> tpcCluster(new Cube::ReconCluster);
        for (Cube::HitSelection::iterator h = onlyTPC.begin();
             h != onlyTPC.end(); ++h) {
            if (!(*h)->HasProperty("DriftVelocity")) {
                if (!problems) CUBE_ERROR << "TPC Hit without drift velocity"
                                          << std::endl;
                ++problems;
                continue;
            }
            Cube::WritableHit hit3d(*(*h));
            TVector3 pos = hit3d.GetPosition();
            double x = pos.X();
            x -= hit3d.GetTime()*hit3d.GetProperty("DriftVelocity");
            pos.SetX(x);
            hit3d.SetPosition(pos);
            hit3d.SetTime(0.0*unit::ns);
            hit3d.AddHit(*h);
            usedHits->push_back(Cube::Handle<Cube::Hit>(new Cube::Hit(hit3d)));
        }
    }

    if (!onlyECal.empty()) {
        Cube::Handle<Cube::ReconCluster> ecalCluster(new Cube::ReconCluster);
        std::map<int,Cube::HitSelection> cells;
        std::map<Cube::Handle<Cube::Hit>,int> madeHitCount;
        for (Cube::HitSelection::iterator h = onlyECal.begin();
             h != onlyECal.end(); ++h) {
            int cell = Cube::Info::ECalModLayCel((*h)->GetIdentifier());
            cells[cell].push_back(*h);
            madeHitCount[*h] = 0;
        }
        for (std::map<int,Cube::HitSelection>::iterator c = cells.begin();
             c != cells.end(); ++c) {
            for (Cube::HitSelection::iterator h1 = c->second.begin();
                 h1 != c->second.end(); ++h1) {
                for (Cube::HitSelection::iterator h2 = h1+1;
                     h2 != c->second.end(); ++h2) {
                    if (Cube::Info::ECalEnd((*h1)->GetIdentifier())
                        == Cube::Info::ECalEnd((*h2)->GetIdentifier())) {
                        continue;
                    }
                    if (Cube::Info::ECalModule((*h1)->GetIdentifier())
                        != Cube::Info::ECalModule((*h2)->GetIdentifier())) {
                        std::cout << "BAD MODULE " << std::endl;
                    }
                    if (Cube::Info::ECalLayer((*h1)->GetIdentifier())
                        != Cube::Info::ECalLayer((*h2)->GetIdentifier())) {
                        std::cout << "BAD LAYER " << std::endl;
                    }
                    if (Cube::Info::ECalCell((*h1)->GetIdentifier())
                        != Cube::Info::ECalCell((*h2)->GetIdentifier())) {
                        std::cout << "BAD LAYER " << std::endl;
                    }
                    double tDiff = (*h2)->GetTime() - (*h1)->GetTime();
                    double vlight = 17.094*unit::cm/unit::ns;
                    if (std::abs(tDiff) > 6.0*unit::meter/vlight) {
                        continue;
                    }
                    TVector3 avg = 0.5*((*h2)->GetPosition()
                                        +(*h1)->GetPosition());
                    TVector3 dif = ((*h2)->GetPosition()
                                    - (*h1)->GetPosition()).Unit();
                    TVector3 pos = avg + (0.5*tDiff*vlight)*dif;
                    Cube::WritableHit hit3d;
                    hit3d.SetIdentifier(c->first);
                    double q = (*h1)->GetCharge() + (*h2)->GetCharge();
                    hit3d.SetCharge(q);
                    hit3d.SetChargeUncertainty(std::sqrt(q));
                    double t = 0.5*((*h1)->GetTime() + (*h2)->GetTime());
                    hit3d.SetTime(t);
                    hit3d.SetPosition(pos);
                    TVector3 size(
                        std::max((*h1)->GetSize().X(),(*h2)->GetSize().X()),
                        std::max((*h1)->GetSize().Y(),(*h2)->GetSize().Y()),
                        std::max((*h1)->GetSize().Z(),(*h2)->GetSize().Z()));
                    hit3d.SetUncertainty(2.0*0.289*size);
                    hit3d.SetSize(size);
                    hit3d.AddHit(*h1);
                    hit3d.AddHit(*h2);
                    ++madeHitCount[*h1];
                    ++madeHitCount[*h2];
                    usedHits->push_back(Cube::Handle<Cube::Hit>(
                                            new Cube::Hit(hit3d)));
                }
            }
        }
        for (Cube::HitSelection::iterator h = onlyECal.begin();
             h != onlyECal.end(); ++h) {
            if (madeHitCount[*h] < 1) usedHits->push_back(*h);
        }
    }

    // Slice the 3DST up by time.  Only fibers in the same time slice are used
    // to build cube hits.
    Cube::Handle<Cube::AlgorithmResult> timeSlice
        = Run<Cube::TimeSlice>(only3DST);
    if (timeSlice) {
        result->AddAlgorithmResult(timeSlice);

        // Process each time slice as a separate cube reconstruction.
        Cube::Handle<Cube::ReconObjectContainer> slices
            = timeSlice->GetObjectContainer("final");

        int slice = 0;
        for (Cube::ReconObjectContainer::iterator obj = slices->begin();
             obj != slices->end(); ++obj) {
            Cube::Handle<Cube::HitSelection> hits = (*obj)->GetHitSelection();
            if (!hits) continue;
            if (hits->empty()) continue;
            // This is building cube hits from fiber hits.  Notice that since
            // this is happening "out-of-band" and we need to make a local
            // copy to allow the THitSelection to be translated into a
            // Cube::AlgorithmResult.
            Cube::HitSelection local("fibers");
            std::copy(hits->begin(), hits->end(), std::back_inserter(local));
            Cube::Handle<Cube::AlgorithmResult> hits3D
                = Run<Cube::Hits3D>(local);
            if (!hits3D) continue;
            // Build the new name for the algorithm (named after the slice
            // number.
            std::ostringstream nm;
            nm << "Slice"
               << std::setfill('0')
               << std::setw(2)
               << slice++;
            hits3D->SetName(nm.str().c_str());
            // Copy the important information to the main result.
            Cube::Handle<Cube::ReconObjectContainer> objects
                = hits3D->GetObjectContainer("final");
            if (objects) {
                std::copy(objects->begin(), objects->end(),
                          std::back_inserter(*finalObjects));
            }
            Cube::Handle<Cube::HitSelection> cubes
                = hits3D->GetHitSelection("used");
            if (cubes) {
                std::copy(cubes->begin(), cubes->end(),
                          std::back_inserter(*usedHits));
            }
            Cube::Handle<Cube::HitSelection> fibers
                = hits3D->GetHitSelection("unused");
            if (fibers) {
                std::copy(fibers->begin(), fibers->end(),
                          std::back_inserter(*unusedHits));
            }
            // Save the sub result.
            result->AddAlgorithmResult(hits3D);
        }
    }

    // Make sure the hit selections don't have duplicates
    std::sort(unusedHits->begin(), unusedHits->end());
    Cube::HitSelection::iterator end
        = std::unique(unusedHits->begin(), unusedHits->end());
    unusedHits->erase(end, unusedHits->end());

    std::sort(usedHits->begin(), usedHits->end());
    end = std::unique(usedHits->begin(), usedHits->end());
    usedHits->erase(end, usedHits->end());

    // Nicely count all the hits (mostly for debugging).
    int simple3DST = 0;
    int composite3DST = 0;
    int simpleTPC = 0;
    int compositeTPC = 0;
    int simpleECal = 0;
    int compositeECal = 0;
    int otherHits = 0;
    for (Cube::HitSelection::iterator h = usedHits->begin();
         h != usedHits->end(); ++h) {
        int id = (*h)->GetIdentifier();
        if (Cube::Info::Is3DST(id)) {
            if ((*h)->GetConstituentCount() < 1) ++simple3DST;
            else ++composite3DST;
        }
        else if (Cube::Info::IsTPC(id)) {
            if ((*h)->GetConstituentCount() < 1) ++simpleTPC;
            else ++compositeTPC;
        }
        else if (Cube::Info::IsECal(id)) {
            if ((*h)->GetConstituentCount() < 1) ++simpleECal;
            else ++compositeECal;
        }
        else {
            ++otherHits;
        }
    }
    CUBE_LOG(0) << "MakeHits3D::Process"
                << " 3DST: " << composite3DST << " ("<< simple3DST << ")"
                << " TPC: " << compositeTPC << " ("<< simpleTPC << ")"
                << " ECal:" << compositeECal << " ("<< simpleECal << ")"
                << " Other: " << otherHits
                << std::endl;

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
