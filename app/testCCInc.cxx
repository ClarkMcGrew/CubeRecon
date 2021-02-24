#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconVertex.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolInitialize.hxx>
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
#include <TDatabasePDG.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <deque>

typedef std::deque<Cube::Handle<Cube::ReconTrack>> MatchedTracks;
namespace {
    struct CMPMatch {
        bool operator () (const MatchedTracks& lhs,
                          const MatchedTracks& rhs) {
            return (lhs.front()->GetPosition().T()
                    < rhs.front()->GetPosition().T());
        }
    };
}

bool histInitialized = false;
TH1F* histMatchDistance = NULL;
TH1F* histMatchAngle = NULL;
TH1F* histMuonTruth = NULL;
TH1F* histMuonReco = NULL;
TH1F* histBackgroundReco = NULL;

/// Assign ECal hit segments to ecal hits.  This is a test bed for later
/// analysis.
bool AnalyzeEvent(Cube::Event& event) {
    /// Flag tht this event should be saved (if not rejected for some cut).
    bool shouldSave = false;

    // The sign of muon that should be selected: neutrino is -1.0, and
    // anti-neutrino is 1.0
    const double signSelection = -1.0;
    const double fiducialCut = 5*unit::cm;

    if (!histInitialized) {
        Cube::Tool::Initialize();
        histInitialized = true;
        histMatchDistance
            = new TH1F("matchDistance",
                       "Intersection distance for 3DST and TPC tracks",
                       50, 0.0, 25.0);
        histMatchAngle
            = new TH1F("matchAngle",
                       "Intersection angle for 3DST and TPC tracks",
                       50, 0.0, 0.50);
        histMuonTruth
            = new TH1F("muonTruth",
                       "True muon momentum of contained interactions",
                       2000, 0.0, 5*unit::GeV);
        histMuonReco
            = new TH1F("muonReco",
                       "Reconstructed muon momentum for contained events",
                       2000, 0.0, 5*unit::GeV);
        histBackgroundReco
            = new TH1F("backgroundReco",
                       "Reconstructed background momentum for contained events",
                       2000, 0.0, 5*unit::GeV);
   }

    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    Cube::Handle<Cube::HitSelection> hits = event.GetHitSelection();
    if (!hits) return false;
    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return false;

    // Collect the TPC tracks.
    Cube::ReconObjectContainer tracksTPC;
    for (Cube::Handle<Cube::ReconObject> obj : *objects) {
        Cube::Handle<Cube::ReconTrack> track = obj;
        if (!track) continue;
        std::string detectors = Cube::Tool::ObjectDetectors(*obj);
        if (detectors.find(":TPC") == std::string::npos) continue;
        tracksTPC.push_back(track);
    }

    // Collect the 3DST tracks and vertices
    Cube::ReconObjectContainer tracks3DST;
    Cube::ReconObjectContainer vertex3DST;
    for (Cube::Handle<Cube::ReconObject> obj : *objects) {
        std::string detectors = Cube::Tool::ObjectDetectors(*obj);
        if (detectors.find(":3DST") == std::string::npos) continue;
        Cube::Handle<Cube::ReconTrack> track = obj;
        if (track) tracks3DST.push_back(track);
        Cube::Handle<Cube::ReconVertex> vertex = obj;
        if (vertex) vertex3DST.push_back(vertex);
    }

    std::vector<MatchedTracks> matchedTracks;
    for (Cube::Handle<Cube::ReconObject> dst : tracks3DST) {
        Cube::Handle<Cube::ReconTrack> tDST = dst;
        double bestDist = 1E+20;
        double bestCos = 1E+20;
        Cube::Handle<Cube::ReconTrack> bestTPC;
        for (Cube::Handle<Cube::ReconObject> tpc : tracksTPC) {
            Cube::Handle<Cube::ReconTrack> tTPC = tpc;
            // Get the drift velocity for this track.  The sign is the
            // important part.
            Cube::Handle<Cube::Hit> hit = tTPC->GetHitSelection()->front();
            double dv = hit->GetProperty("DriftVelocity");
            // The TPC track time should always be zero.
            double dt = tDST->GetBack()->GetPosition().T()
                - tTPC->GetPosition().T();
            TVector3 backDST = tDST->GetBack()->GetPosition().Vect();
            TVector3 dirDST =  tDST->GetBack()->GetDirection();
            TVector3 frontTPC = tTPC->GetFront()->GetPosition().Vect();
            frontTPC.SetX(frontTPC.X() + dt*dv);
            TVector3 backTPC = tTPC->GetBack()->GetPosition().Vect();
            backTPC.SetX(backTPC.X() + dt*dv);
            double f1 = (frontTPC-backDST).Mag();
            double b1 = (backTPC-backDST).Mag();
            double dCos = 2.0;
            double dMiss = 1E+20;
            if (f1 < b1) {
                dCos = std::abs(dirDST*tTPC->GetFront()->GetDirection());
                double moveDist = (frontTPC-backDST)*dirDST;
                TVector3 diff = frontTPC - (backDST + moveDist * dirDST);
                dMiss = diff.Mag();
            }
            else {
                dCos = std::abs(dirDST*tTPC->GetBack()->GetDirection());
                double moveDist = (backTPC-backDST)*dirDST;
                TVector3 diff = backTPC - (backDST + moveDist * dirDST);
                dMiss = diff.Mag();
            }
            if (dCos < 0.7) continue;
            if (bestDist < dMiss) continue;
            bestDist = dMiss;
            bestCos = dCos;
            bestTPC = tTPC;
        }
        if (bestDist < 30.0 && bestCos > 0.7) {
            MatchedTracks tmp;
            tmp.push_back(tDST);
            tmp.push_back(bestTPC);
            matchedTracks.push_back(tmp);
            histMatchDistance->Fill(bestDist);
            histMatchAngle->Fill(std::acos(bestCos));
        }
    }

    // Now check if there are 3DST tracks that need to be merged
    for (MatchedTracks& match : matchedTracks) {
        bool matchFound = false;
        do {
            Cube::Handle<Cube::ReconTrack> tDST = match.front();
            double bestDist = 1E+20;
            double bestCos = 1E+20;
            TVector3 posDST = tDST->GetFront()->GetPosition().Vect();
            TVector3 dirDST = tDST->GetFront()->GetDirection();
            int trajDST = Cube::Tool::MainTrajectory(event,*tDST);
            matchFound = false;
            for (Cube::Handle<Cube::ReconObject> dst : tracks3DST) {
                Cube::Handle<Cube::ReconTrack> cand = dst;
                if (!cand) continue;
                if (tDST == cand) continue;
                if (tDST->GetPosition().T()<=cand->GetPosition().T()) continue;
                TVector3 posCand = tDST->GetBack()->GetPosition().Vect();
                TVector3 dirCand = tDST->GetBack()->GetDirection();
                double dCos = dirCand * dirDST;
                // if (dCos < 0.8) continue;
                double dist = (posDST - posCand).Mag();
                // if (dist > 10*unit::cm) continue;
                int trajCand = Cube::Tool::MainTrajectory(event,*cand);
                if (trajCand != trajDST) continue;
                match.push_front(cand);
                matchFound = true;
                break;
            }
        } while (matchFound);
    }

    // The number of matched pairs between the 3DST and the TPC.
    std::sort(matchedTracks.begin(), matchedTracks.end(),CMPMatch());
    double lastTime = -1E+20;
    std::vector<MatchedTracks> mipCandidates;
    MatchedTracks bestMatch;
    double bestMomentum = 0;
    for (MatchedTracks match : matchedTracks) {
        int traj3DST = Cube::Tool::MainTrajectory(event,*match.front());
        int trajTPC =  Cube::Tool::MainTrajectory(event,*match.back());
        if (trajTPC < 0) {
            std::cout << "BAD TPC TRACK MATCH" << std::endl;
            continue;
        }
        int primTPC = Cube::Tool::PrimaryId(event,trajTPC);
        if (primTPC < 0) {
            std::cout << "BAD TPC PRIMARY" << std::endl;
            continue;
        }
        double deltaTime = match.front()->GetPosition().T() - lastTime;
        if (lastTime > 0 && deltaTime > 20*unit::ns) {
            if (bestMatch.size() > 0) {
                std::cout << "Best Momentum " << bestMomentum << std::endl;
                mipCandidates.push_back(bestMatch);
            }
            bestMomentum = 0.0;
            bestMatch.clear();
        }
        else if (lastTime < 0) {
            bestMomentum = 0.0;
            bestMatch.clear();
        }
        lastTime = match.front()->GetPosition().T();
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(primTPC);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO TPC PRIMARY" << std::endl;
            continue;
        }
        double tpcMom = t->second->GetInitialMomentum().P();
        if (tpcMom < bestMomentum) {
            continue;
        }

        int pdgCode =  t->second->GetPDGCode();
        const TDatabasePDG *pdg = TDatabasePDG::Instance();

        // Get the charge for the track.  This is cheating since it should
        // come from the curvature, but the current fitter doesn't do it.
        double charge = 0.0;
        switch (pdgCode) {
        case 13: charge = -1.0; break;
        case -13: charge = 1.0; break;
        case 11: charge = -1.0; break;
        case -11: charge = 1.0; break;
        default:
            if (pdgCode > 0) charge = 1.0;
            else charge = -1.0;
            break;
        }

        // The signSelection value is set at the top of this routine.
        if (signSelection*charge < 0.0) continue;

        // This is cheating, but the id is very good.
        if (pdgCode == 11 || pdgCode == -11) continue;

        // Approximate particle ID.  This is a bit better than reality.
        if (pdgCode == 2212 && tpcMom < 1.2*unit::GeV) continue;

        // Apply a fiducial volume.
        double fidDist
            = Cube::Tool::ContainedPoint(
                match.front()->GetFront()->GetPosition().Vect());
        if (fidDist < fiducialCut) continue;

        bestMatch = match;
        bestMomentum = tpcMom;
#ifdef LOUD_AND_PROUD
        std::cout << "MATCH " << match.front()->GetPosition().T()
                  << " " << traj3DST
                  << " " << trajTPC
                  << " " << primTPC
                  << " " << bestMomentum
                  << " " << pdgCode
                  << " " << charge
                  << std::endl;
#endif
    }

    // Collect information about the truth!
    std::vector<int> allPrimaries = Cube::Tool::AllPrimaries(event);
    for (int prim : allPrimaries) {
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(prim);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO PRIMARY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> traj = t->second;
        // Check the particle type.
        if (traj->GetPDGCode() != 13 && traj->GetPDGCode() != -13) continue;
        // Check the sign.
        // if (signSelection*traj->GetPDGCode() > 0) continue;
        // Check the fiducial volume.
        double fidTruth
            = Cube::Tool::ContainedPoint(
                traj->GetInitialPosition().Vect());
        if (fidTruth < fiducialCut) continue;
        // Get the momentum
        double momTruth = traj->GetInitialMomentum().P();
        if (momTruth > 4.99*unit::GeV) momTruth = 4.99*unit::GeV;
        histMuonTruth->Fill(momTruth);
        std::cout << "Truth " << momTruth << " " << fidTruth << std::endl;
    }

    // Now the mipCandidates vector has all of the 3DST and TPC candidates.
    std::cout << "MIP Candidates " << mipCandidates.size() << std::endl;
    for (MatchedTracks cand : mipCandidates) {
        bool signal = true;
        int traj3DST = Cube::Tool::MainTrajectory(event,*cand.front());
        if (cand.size() > 1) {
            int trajTPC =  Cube::Tool::MainTrajectory(event,*cand.back());
            if (traj3DST != trajTPC) signal = false;
        }
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(traj3DST);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO TRAJECTORY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> traj = t->second;
        // Check the particle type.
        if (traj->GetPDGCode() == 11 || traj->GetPDGCode() == -11) continue;
        if (traj->GetPDGCode() != 13 && traj->GetPDGCode() != -13) signal=false;
        if (signSelection*traj->GetPDGCode() > 0) signal = false;
        double fidTruth
            = Cube::Tool::ContainedPoint(
                traj->GetInitialPosition().Vect());
        if (fidTruth < 0.0*unit::cm) signal = false;
        double momReco = traj->GetInitialMomentum().P();
        if (momReco > 4.99*unit::GeV) momReco = 4.99*unit::GeV;
        if (signal) {
            histMuonReco->Fill(momReco);
            std::cout << "Sgnl " << traj->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << std::endl;
        }
        else {
            histBackgroundReco->Fill(momReco);
            std::cout << "Back " << traj->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << std::endl;
            if (std::abs(momReco - 1086) < 2) {
                std::cout << traj->GetPDGCode() << std::endl;
                std::cout << traj->GetInitialPosition().X() << std::endl;
                std::cout << traj->GetInitialPosition().Y() << std::endl;
                std::cout << traj->GetInitialPosition().Z() << std::endl;
                std::cout << traj->GetInitialPosition().T() << std::endl;
                std::cout << traj->GetInitialMomentum().X() << std::endl;
                std::cout << traj->GetInitialMomentum().Y() << std::endl;
                std::cout << traj->GetInitialMomentum().Z() << std::endl;
                std::cout << traj->GetInitialMomentum().P() << std::endl;
            }
        }
    }

    // Uncomment to save all of the CC0pi events.
    // shouldSave = true;
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
