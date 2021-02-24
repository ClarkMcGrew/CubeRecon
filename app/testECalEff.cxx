#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>
#include <ToolRecon.hxx>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGeoManager.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

bool histInitialized = false;
TH1F* histECalSegTime = NULL;
TH1F* histECalHitTime = NULL;
TH1F* histECalHitZ = NULL;
TH1F* histECalDeltaTime = NULL;
TH1F* histECalDeltaX = NULL;
TH1F* histOtherHitTime = NULL;
TH1F* histMuonOverlap = NULL;
TH1F* histUpstreamOverlap = NULL;
TH1F* histHitsSpill = NULL;
TH1F* histHitOverlap = NULL;
TH1F* histMuonSpill = NULL;
TH1F* histUpstreamSpill = NULL;
TH1F* histTrajectoryTotal = NULL;
TH1F* histTrajectoryLarge = NULL;
TH1F* histTrajectorySmall = NULL;
TH1F* histTrajectorySize = NULL;
TH1F* histInteractionCount = NULL;

/// Assign ECal hit segments to ecal hits.  This is a test bed for later
/// analysis.
bool AnalyzeEvent(Cube::Event& event) {
    /// Flag tht this event should be saved (if not rejected for some cut).
    bool shouldSave = false;

    if (!histInitialized) {
        histInitialized = true;
        histECalSegTime = new TH1F("ecalSegTime",
                                   "ECal Segment times",
                                   1000, 0.0, 20000.0);
        histECalHitTime = new TH1F("ecalHitTime",
                                   "ECal Hit times",
                                   1000, 0.0, 20000.0);
        histECalHitZ = new TH1F("ecalHitZ",
                                "ECal Hit Z",
                                1000, 60000.0, 70000.0);
        histECalDeltaTime = new TH1F("ecalDeltaTime",
                                     "ECal Delta times",
                                     500, -100.0, 500.0);
        histECalDeltaX = new TH1F("ecalDeltaX",
                                  "ECal Delta X",
                                  120, 0.0, 600.0);
        histOtherHitTime = new TH1F("otherHitTime",
                                    "Other Hit times",
                                    1000, 0.0, 20000.0);
        histMuonOverlap = new TH1F("muonOverlap",
                                   "ECal muons in spill with overlaps",
                                   20, 0.0, 1.0);
        histUpstreamOverlap = new TH1F("upstreamOverlap",
                                   "Upstream ECal muons in spill with overlaps",
                                   20, 0.0, 1.0);
        histHitsSpill = new TH1F("hitsSpill",
                                 "ECal hits in spill",
                                 65, 0.0,1300.0);
        histHitOverlap = new TH1F("hitOverlap",
                                  "ECal hits in spill with overlaps",
                                  20, 0.0, 1.0);
        histMuonSpill = new TH1F("muonSpill",
                                 "Muons in ECal per spill",
                                 60, 0.0, 60.0);
        histUpstreamSpill = new TH1F("upstreamSpill",
                                     "Muons in upstream ECal per spill",
                                     40, 0.0, 40.0);
        histTrajectoryTotal = new TH1F("trajectoryTotal",
                                       "Total ECal tracks",
                                       60, 0.0, 6000.0);
        histTrajectoryLarge = new TH1F("trajectoryLarge",
                                       "ECal tracks with more than 2 hits",
                                       90, 0.0, 900.0);
        histTrajectorySmall = new TH1F("trajectorySmall",
                                       "ECal tracks with less than 3 hits",
                                       60, 0.0, 6000.0);
        histTrajectorySize = new TH1F("trajectorySize",
                                      "ECal muon trajectory size",
                                      50, 0.0, 100.0);
        histInteractionCount = new TH1F("interactionCount",
                                        "Iteractions per spill hitting ECal",
                                        50, 0.0, 100.0);
    }

    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    Cube::Handle<Cube::HitSelection> hits = event.GetHitSelection();
    if (!hits) return false;

    // Collect all of the ECal Hits for easy use later.
    Cube::HitSelection ecalOnly;
    for (Cube::Handle<Cube::Hit> hit : *hits) {
        if (Cube::Info::Is3DST((hit->GetIdentifier()))) {
            histOtherHitTime->Fill(hit->GetTime());
            continue;
        }
        if (!Cube::Info::IsECal(hit->GetIdentifier())) continue;
        if (hit->GetTime() < 1.0) continue;
        if (hit->GetConstituentCount() < 2) continue;
        histECalHitTime->Fill(hit->GetTime());
        histECalHitZ->Fill(hit->GetPosition().Z());
        ecalOnly.push_back(hit);
    }

    // Collect all of the ECal hit segments.
    for (Cube::Event::G4HitContainer::iterator s = event.G4Hits.begin();
         s != event.G4Hits.end(); ++s) {
        Cube::Handle<Cube::G4Hit> seg = s->second;
        if (!seg) {
            CUBE_ERROR << "Invalid track segment" << std::endl;
            continue;
        }
        TVector3 pnt = 0.5*(seg->GetStart().Vect() + seg->GetStop().Vect());
        TGeoNode* node  = gGeoManager->FindNode(pnt.X(),pnt.Y(),pnt.Z());
        if (!node) continue;
        std::string volumeName = gGeoManager->GetPath();
        if (volumeName.find("kloe_calo_volume") == std::string::npos) continue;
    }

    double total = 0;
    double overlap = 0;
    double interactionOverlap = 0;
    std::set<Cube::Handle<Cube::Hit>> hitWithOverlap;
    std::map<int,std::set<Cube::Handle<Cube::Hit>>> trajToHit;
    for (Cube::Handle<Cube::Hit> hit : ecalOnly) {
        std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
            = Cube::Tool::HitG4Hits(event,hit);
        for (Cube::Handle<Cube::G4Hit> seg : g4Hits) {
            trajToHit[seg->GetPrimaryId()].insert(hit);
            histECalSegTime->Fill(seg->GetStart().T());
            histECalDeltaTime->Fill(seg->GetStart().T()-hit->GetTime());
            histECalDeltaX->Fill(
                std::abs(seg->GetStart().X()-hit->GetPosition().X()));
        }

        double deltaX = 0;
        double deltaT = 0.0;
        double interactionDeltaT = 0.0;
        if (g4Hits.size() > 1) {
            Cube::Handle<Cube::G4Hit> seg = g4Hits.front();
            TVector3 pnt = 0.5*(seg->GetStart().Vect() + seg->GetStop().Vect());
            double x1 = pnt.X();
            double t1 = seg->GetStart().T();
            int trajId1 = Cube::Tool::PrimaryId(event,seg->GetPrimaryId());
            Cube::Event::G4TrajectoryContainer::iterator traj1
                = event.G4Trajectories.find(trajId1);
            for (Cube::Handle<Cube::G4Hit> seg2: g4Hits) {
                int trajId2 = Cube::Tool::PrimaryId(event,seg2->GetPrimaryId());
                Cube::Event::G4TrajectoryContainer::iterator traj2
                    = event.G4Trajectories.find(trajId2);
                if (traj1->second && traj2->second) {
                    double dt = traj2->second->GetInitialPosition().T()
                        - traj1->second->GetInitialPosition().T();
                    interactionDeltaT = std::max(interactionDeltaT,
                                                 std::abs(dt));
                }
                else {
                    std::cout << "Missing trajectory" << std::endl;
                    throw std::runtime_error("Missing trajectory");
                }
                pnt = 0.5*(seg2->GetStart().Vect() + seg2->GetStop().Vect());
                double x2 = pnt.X();
                double t2 = seg2->GetStart().T();
                deltaX = std::max(deltaX, std::abs(x2-x1));
                deltaT = std::max(deltaT, std::abs(t2-t1));
            }
        }
        ++total;
        if (interactionDeltaT > 0.50) {
            ++interactionOverlap;
        }
        if (std::abs(deltaX) > 50.0*unit::cm
            || std::abs(deltaT) > 5.0*unit::ns
            || interactionDeltaT > 0.50*unit::ns) {
            ++overlap;
            hitWithOverlap.insert(hit);
        }
    }
    if (total<1) total = 1;
    std::cout << "Hits " << total
              << " " << overlap
              << " " << hitWithOverlap.size()
              << " " << interactionOverlap
              << " " << overlap/total
              << std::endl;
    histHitOverlap->Fill(overlap/total);
    histHitsSpill->Fill(total);

    double trajectoryTotal = 0;
    double trajectoryLarge = 0;
    double trajectorySmall = 0;
    double muonCount = 0;
    double muonOverlaps = 0;
    double muonUpstreamCount = 0;
    double muonUpstreamOverlaps = 0;
    std::vector<int> timeVector;
    for (auto elem : trajToHit) {
        int trId = elem.first;
        if (trId < 0) {
            continue;
        }
        Cube::Event::G4TrajectoryContainer::iterator trajElem
            = event.G4Trajectories.find(trId);
        if (trajElem == event.G4Trajectories.end()) {
            std::cout << "traj id not found" << std::endl;
            continue;
        }
        int primId = Cube::Tool::PrimaryId(event,trId);
        Cube::Event::G4TrajectoryContainer::iterator primElem
            = event.G4Trajectories.find(primId);
        if (primElem != event.G4Trajectories.end() && primElem->second) {
            int integerTime = primElem->second->GetInitialPosition().T();
            timeVector.push_back(integerTime);
        }
        std::set<Cube::Handle<Cube::Hit>>& hitSet = elem.second;
        if (hitSet.size() < 1) {

            continue;
        }
        if (hitSet.size() > 2) ++trajectoryLarge;
        else if (hitSet.size() > 0) ++trajectorySmall;
        ++trajectoryTotal;
        if (std::abs(trajElem->second->GetPDGCode()) != 13) {
            continue;
        }
        histTrajectorySize->Fill(hitSet.size());
        ++muonCount;
        int trajOverlaps = 0;
        int upstreamCount = 0;
        int upstreamOverlaps = 0;
        for (Cube::Handle<Cube::Hit> hit : elem.second) {
            if (hit->GetPosition().Z() < 63900) ++upstreamCount;
            if (hitWithOverlap.find(hit) == hitWithOverlap.end()) {
                continue;
            }
            ++trajOverlaps;
            if (hit->GetPosition().Z() < 63900) ++upstreamOverlaps;
        }
        if (trajOverlaps > 0) ++muonOverlaps;
        if (upstreamCount > 0) ++muonUpstreamCount;
        if (upstreamOverlaps > 0) ++muonUpstreamOverlaps;
    }
    if (muonCount < 1) muonCount = 1;
    if (muonUpstreamCount < 1) muonUpstreamCount = 1;
    std::cout << "Trajectories " << trajectoryTotal
              << " " << muonCount
              << " " << muonOverlaps
              << " " << muonOverlaps/muonCount
              << " " << muonUpstreamOverlaps/muonUpstreamCount
              << std::endl;
    histMuonOverlap->Fill(muonOverlaps/muonCount);
    histUpstreamOverlap
        ->Fill(muonUpstreamOverlaps/muonUpstreamCount);
    histMuonSpill->Fill(muonCount);
    histUpstreamSpill->Fill(muonUpstreamCount);
    histTrajectorySmall->Fill(trajectorySmall);
    histTrajectoryLarge->Fill(trajectoryLarge);
    histTrajectoryTotal->Fill(trajectoryTotal);

    std::sort(timeVector.begin(),timeVector.end());
    std::vector<int>::iterator iend = std::unique(timeVector.begin(),
                                                  timeVector.end());
    timeVector.erase(iend,timeVector.end());
    histInteractionCount->Fill(timeVector.size());

    // Uncomment to save all of the CC0pi events.
    shouldSave = true;
    return shouldSave;
}

int main(int argc, char** argv) {
    int maxEntries = 1E+8; // Maximum to process.
    int firstEntry = 0;
    std::string outputName;

    while (true) {
        int c = getopt(argc,argv,"n:o:s:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
            break;
        }
        case 'o': {
            outputName = optarg;
            break;
        }
        case 's': {
            std::istringstream tmp(optarg);
            tmp >> firstEntry;
            break;
        }
        default: {
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
                      << "-o <number>  : Output file"
                      << std::endl
                      << "-s <number>  : Skip <number> entries"
                      << std::endl
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl;
            exit(1);
        }
        }
    }

    std::vector<std::string> inputNames;
    if (argc <= optind) throw std::runtime_error("Missing input file");
    while (optind < argc) {
        inputNames.push_back(argv[optind++]);
        std::cout << "Input Name " << inputNames.back() << std::endl;
    }

    if (outputName.empty()) {
        std::cout << "NO OUTPUT FILE!!!!" << std::endl;
    }

    // Attach to the input tree.
    std::unique_ptr<TFile> inputFile(
        new TFile(inputNames.back().c_str(),"old"));
    if (!inputFile->IsOpen()) throw std::runtime_error("Input file not open");

    TGeoManager* geom
        = dynamic_cast<TGeoManager*>(inputFile->Get("CubeReconGeometry"));
    if (!geom) {
        std::cout << "No geometry for file" << std::endl;
    }

    /// Attach to the input tree.
    std::unique_ptr<TChain> inputChain(new TChain("CubeEvents"));
    for (int i = 0; i<inputNames.size(); ++i) {
        inputChain->AddFile(inputNames[i].c_str());
    }
    Cube::Event *inputEvent = NULL;
    inputChain->SetBranchAddress("Event",&inputEvent);

    // Open the output file
    std::unique_ptr<TFile> outputFile;
    if (!outputName.empty()) {
        std::cout << "Open Output File: " << outputName << std::endl;
        outputFile.reset(new TFile(outputName.c_str(),"recreate"));
    }
    TTree *outputTree = new TTree("CubeEvents","Reconstructed Event");
    static Cube::Event *outputEvent = inputEvent;
    outputTree->Branch("Event",&outputEvent);

    // Loop through the events.
    int totalEntries = inputChain->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputChain->GetEntry(entry);
        outputEvent = inputEvent;
        std::cout << "Process event "
                  << entry
                  << "/" << inputEvent->GetRunId()
                  << "/" << inputEvent->GetEventId() << std::endl;
        bool save = AnalyzeEvent(*inputEvent);
        if (save) outputTree->Fill();
    }

    if (outputFile) {
        outputFile->Write();
        outputFile->Close();
    }

    return 0;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
