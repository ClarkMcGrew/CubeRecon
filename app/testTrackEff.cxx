#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>
#include <CubeUnits.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

bool histInitialized = false;
TH1F* histPrimMuonMomentum = NULL;
TH1F* histPrimMuonLength = NULL;
TH1F* histPrimMuonDCos = NULL;
TH2F* histPrimMuonMomDCos = NULL;
TH2F* histPrimMuonMomLen = NULL;
TH2F* histPrimMuonDir = NULL;
TH1F* histRecoMuonMomentum = NULL;
TH1F* histRecoMuonLength = NULL;
TH1F* histRecoMuonDCos = NULL;
TH1F* histRecoMuonComplete = NULL;
TH2F* histRecoMuonLenComplete = NULL;
TH2F* histRecoMuonMomDCos = NULL;
TH2F* histRecoMuonDir = NULL;
TH1F* histEffMuonMomentum = NULL;
TH1F* histEffMuonLength = NULL;
TH1F* histEffMuonDCos = NULL;
TH2F* histEffMuonMomDCos = NULL;
TH2F* histEffMuonDir = NULL;

TH1F* histPrimProtonMomentum = NULL;
TH1F* histRecoProtonMomentum = NULL;
TH1F* histEffProtonMomentum = NULL;
TH2F* histPrimProtonDir = NULL;
TH2F* histRecoProtonDir = NULL;
TH2F* histEffProtonDir = NULL;

TH1F* histPrimPionMomentum = NULL;
TH1F* histRecoPionMomentum = NULL;
TH1F* histEffPionMomentum = NULL;
TH2F* histPrimPionDir = NULL;
TH2F* histRecoPionDir = NULL;
TH2F* histEffPionDir = NULL;

TH1F* histPositionX = NULL;
TH1F* histPositionY = NULL;
TH1F* histPositionZ = NULL;

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
bool AnalyzeEvent(Cube::Event& event) {
    const double fidDist = 25*unit::cm;

    if (!histInitialized) {
        histInitialized = true;
        histPositionX = new TH1F("positionX",
                                 "Trajectory X Position",
                                 400, -2000.0, 2000.0);
        histPositionY = new TH1F("positionY",
                                 "Trajectory Y Position",
                                 1000, -5000.0, 0.0);
        histPositionZ = new TH1F("positionZ",
                                 "Trajectory Z Position",
                                 1000, 20000.0, 50000.0);
        histPrimMuonMomentum = new TH1F("primMuonMomentum",
                                       "True Momentum of Muon",
                                       50,0.0,500.0);
        histPrimMuonLength = new TH1F("primMuonLength",
                                      "True Length of Muon",
                                      50,0.0,2000.0);
        histPrimMuonDCos = new TH1F("primMuonDCos",
                                    "Angle Relative to Closest Track (Muon)",
                                    20,-1.0,1.0);
        histPrimMuonMomLen = new TH2F("primMuonMomLen",
                                      "True Length vs Momentum of Muon",
                                      50,0.0,500.0,
                                      50,0.0,2000.0);
        histPrimMuonMomDCos = new TH2F("primMuonMomDCos",
                                       "Angle to Neighbor vs Momentum of Muon",
                                       10,0.0,500.0,
                                       10,-1.0,1.0);
        histPrimMuonDir = new TH2F("primMuonDir",
                                   "True Cosine vs Momentum for Muons",
                                   10,0.0,500.0,
                                   10,0.0,1.0);

        histRecoMuonMomentum = new TH1F("recoMuonMomentum",
                                        "True Momentum of Reconstructed Muon",
                                        50,0.0,500.0);
        histRecoMuonLength = new TH1F("recoMuonLength",
                                        "True Length of Reconstructed Muon",
                                        50,0.0,2000.0);
        histRecoMuonDCos = new TH1F("recoMuonDCos",
                                    "Angle Relative to Closest Track (Muon)",
                                    20,-1.0,1.0);
        histRecoMuonComplete = new TH1F("recoMuonComplete",
                                        "Completeness of the Track (Muon)",
                                        100,0.0,1.0);
        histRecoMuonLenComplete = new TH2F("recoMuonLenComplete",
                                           "Completeness vs Length (Muon)",
                                           10,0.0,2000.0,
                                           20,0.0,1.0);
        histRecoMuonMomDCos = new TH2F("recoMuonMomDCos",
                                       "Angle to Neighbor vs Momentum of Muon",
                                       10,0.0,500.0,
                                       10,-1.0,1.0);
        histRecoMuonDir = new TH2F("recoMuonDir",
                                   "True Cosine vs Momentum for Reco Muons",
                                   10,0.0,500.0,
                                   10,0.0,1.0);

        histEffMuonMomentum = new TH1F("effMuonMomentum",
                                        "Efficiency at True Muon Momentum",
                                        50,0.0,500.0);
        histEffMuonLength = new TH1F("effMuonLength",
                                        "Efficiency at True Muon Length",
                                        50,0.0,2000.0);
        histEffMuonDCos = new TH1F("effMuonDCos",
                                    "Efficiency Relative to Neighbors (Muon)",
                                    20,-1.0,1.0);
        histEffMuonMomDCos = new TH2F("effMuonMomDCos",
                                      "Eff. for Neighbor vs Momentum of Muon",
                                      10,0.0,500.0,
                                      10,-1.0,1.0);
        histEffMuonDir = new TH2F("effMuonDir",
                                  "Muon Efficiency for True Cosine vs Momentum",
                                   10,0.0,500.0,
                                   10,0.0,1.0);

        histPrimProtonMomentum = new TH1F("primProtonMomentum",
                                       "True Momentum of Proton",
                                       50,200.0,700.0);
        histRecoProtonMomentum = new TH1F("recoProtonMomentum",
                                        "True Momentum of Reconstructed Proton",
                                        50,200.0,700.0);
        histEffProtonMomentum = new TH1F("effProtonMomentum",
                                        "Efficiency at True Proton Momentum",
                                        50,200.0,700.0);
        histPrimProtonDir = new TH2F("primProtonDir",
                                   "True Cosine vs Momentum for Protons",
                                   10,200.0,700.0,
                                   10,0.0,1.0);
        histRecoProtonDir = new TH2F("recoProtonDir",
                                   "True Cosine vs Momentum for Reco Protons",
                                   10,200.0,700.0,
                                   10,0.0,1.0);
        histEffProtonDir = new TH2F(
            "effProtonDir",
            "Proton Efficiency for True Cosine vs Momentum",
            10,200.0,700.0,
            10,0.0,1.0);

        histPrimPionMomentum = new TH1F("primPionMomentum",
                                       "True Momentum of Pi+",
                                       50,0.0,500.0);
        histRecoPionMomentum = new TH1F("recoPionMomentum",
                                        "True Momentum of Reconstructed Pi+",
                                        50,0.0,500.0);
        histEffPionMomentum = new TH1F("effPionMomentum",
                                        "Efficiency at True Pi+ Momentum",
                                        50,0.0,500.0);
        histPrimPionDir = new TH2F("primPionDir",
                                   "True Cosine vs Momentum for Pi+",
                                   10,0.0,500.0,
                                   10,0.0,1.0);
        histRecoPionDir = new TH2F("recoPionDir",
                                   "True Cosine vs Momentum for Reco Pi+",
                                   10,0.0,500.0,
                                   10,0.0,1.0);
        histEffPionDir = new TH2F("effPionDir",
                                  "Pi+ Efficiency for True Cosine vs Momentum",
                                   10,0.0,500.0,
                                   10,0.0,1.0);
    }

    int muonFound = 0;

    // Make sure that the tracks are in the "fiducial", and not going in
    // opposite directions.
    std::vector<TVector3> directions;
    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    for (Cube::Event::G4TrajectoryContainer::iterator tr = trajectories.begin();
         tr != trajectories.end(); ++tr) {
        // only count primaries.
        if (tr->second->GetParentId() >= 0) continue;
        TVector3 pnt = tr->second->GetInitialPosition().Vect();
        histPositionX->Fill(pnt.X());
        histPositionY->Fill(pnt.Y());
        histPositionZ->Fill(pnt.Z());
        double fid= Cube::Tool::ContainedPoint(pnt);
        if (fid < fidDist) {
            std::cout << "Trajectory out of fiducial " << fid << std::endl;
            return false;
        }
        directions.push_back(tr->second->GetInitialMomentum().Vect().Unit());
    }

    for (Cube::Event::G4TrajectoryContainer::iterator tr = trajectories.begin();
         tr != trajectories.end(); ++tr) {
        // only count primaries.
        if (tr->second->GetParentId() >= 0) continue;
        double fid= Cube::Tool::ContainedPoint(
            tr->second->GetInitialPosition().Vect());
        if (fid < fidDist) {
            std::cout << "Count out of fiducial " << fid << std::endl;
            return false;
        }
        int trackId = tr->second->GetTrackId();
        int pdgCode = std::abs(tr->second->GetPDGCode());
        std::vector<Cube::Handle<Cube::G4Hit>> segments
            = Cube::Tool::TrajectoryG4Hits(event,trackId);
        if (segments.empty()) continue;
        TVector3 head = segments.front()->GetStart().Vect();
        TVector3 tail = segments.back()->GetStop().Vect();
        double len = (tail-head).Mag();
        TVector3 dir = tr->second->GetInitialMomentum().Vect().Unit();
        double mom = tr->second->GetInitialMomentum().P();
        // Make sure the track is not back to back with another.
        double dCos = 0.0;
        for (std::vector<TVector3>::iterator d1 = directions.begin();
             d1 != directions.end(); ++d1) {
            double d = (*d1) * dir;
            if (d > 0.999) continue; // Not the same track!
            if (std::abs(d) > std::abs(dCos)) dCos = d;
        }
        int pdgMuon = 13;
        if (pdgCode == pdgMuon) {
            if (dCos < 0.9 && dCos > -0.50) {
                muonFound += 10;
                histPrimMuonLength->Fill(len);
                histPrimMuonMomLen->Fill(mom,len);
                histPrimMuonMomentum->Fill(mom);
                histPrimMuonDir->Fill(mom,std::abs(dir.Z()));
            }
            histPrimMuonMomDCos->Fill(mom,dCos);
            if (len > 50) {
                histPrimMuonDCos->Fill(dCos);
            }
        }
        int pdgProton = 2212;
        if (pdgCode == pdgProton) {
            if (dCos < 0.9 && dCos > -0.50) {
                histPrimProtonMomentum->Fill(mom);
                histPrimProtonDir->Fill(mom,std::abs(dir.Z()));
            }
        }
        int pdgPion = 211;
        if (pdgCode == pdgPion) {
            if (dCos < 0.9 && dCos > -0.50) {
                histPrimPionMomentum->Fill(mom);
                histPrimPionDir->Fill(mom,std::abs(dir.Z()));
            }
        }
    }

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) {
        std::cout << "No reconstruction objects" << std::endl;
        return false;
    }

    int pdgMuon = 13;
    int pdgProton = 2212;
    int pdgPion = 211;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) continue;
        // Check that the hits are in the 3DST.
        Cube::Handle<Cube::HitSelection> hits = track->GetHitSelection();
        if (hits->empty()) continue;
        if (!Cube::Info::Is3DST(hits->front()->GetIdentifier())) continue;
        int mainTraj = Cube::Tool::MainTrajectory(event,*track);
        if (mainTraj<0) {
            continue;
        }
        // Reject tracks from secondary particles.
        int primTraj = Cube::Tool::PrimaryId(event,mainTraj);
        if (primTraj != mainTraj) {
            continue;
        }
        // Get the trajectory and check the fiducial for the starting point.
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[primTraj];
        if (!traj) throw;
        // Get the trajectory information for convenience.
        int trackId = traj->GetTrackId();
        int pdgCode = std::abs(traj->GetPDGCode());
        std::vector<Cube::Handle<Cube::G4Hit>> segments
            = Cube::Tool::TrajectoryG4Hits(event,trackId);
        if (segments.empty()) continue;
        TVector3 head = segments.front()->GetStart().Vect();
        TVector3 tail = segments.back()->GetStop().Vect();
        double len = (tail-head).Mag();
        TVector3 dir = traj->GetInitialMomentum().Vect().Unit();
        double mom = traj->GetInitialMomentum().P();
        // Get the purity and completeness of the track.
        double purity = Cube::Tool::MainPurity(event,*track);
        double completeness = Cube::Tool::MainCompleteness(event,*track);
        // Make sure the track is not back to back with another.
        double dCos = 0.0;
        for (std::vector<TVector3>::iterator d1 = directions.begin();
             d1 != directions.end(); ++d1) {
            double d = (*d1) * dir;
            if (d > 0.999) continue;
            if (std::abs(d) > std::abs(dCos)) dCos = d;
        }
        // Collect the info for the muon, pion and proton.
        if (pdgCode == pdgMuon) {
            if (dCos < 0.9 && dCos > -0.50) {
                muonFound += 1000;
                histRecoMuonMomentum->Fill(mom);
                histRecoMuonDir->Fill(mom,std::abs(dir.Z()));
                histRecoMuonLength->Fill(len);
                histRecoMuonComplete->Fill(completeness);
                histRecoMuonLenComplete->Fill(len,completeness);
            }
            histRecoMuonMomDCos->Fill(mom,dCos);
            if (len > 50) {
                histRecoMuonDCos->Fill(dCos);
            }
            pdgMuon = -1;
        }
        if (pdgCode == pdgProton) {
            if (dCos < 0.9 && dCos > -0.50) {
                histRecoProtonMomentum->Fill(mom);
                histRecoProtonDir->Fill(mom,std::abs(dir.Z()));
            }
            pdgProton = -1;
        }
        if (pdgCode == pdgPion) {
            if (dCos < 0.9 && dCos > -0.50) {
                histRecoPionMomentum->Fill(mom);
                histRecoPionDir->Fill(mom,std::abs(dir.Z()));
            }
            pdgPion = -1;
        }
    }

    if (muonFound > 0 && muonFound < 1000) {
        std::cout << "Save bad reco" << std::endl;
        return true;
    }

    return false;
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
    std::unique_ptr<TFile> inputFile(new TFile(inputNames.back().c_str(),"old"));
    if (!inputFile->IsOpen()) throw std::runtime_error("Input file not open");

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
    std::cout << firstEntry << " " << maxEntries << " "  << totalEntries << std::endl;
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputChain->GetEntry(entry);
        outputEvent = inputEvent;
        std::cout << "Process event "
                  << entry
                  << "/" << inputEvent->GetRunId()
                  << "/" << inputEvent->GetEventId()
                  << " out of " << totalEntries << std::endl;
        bool save = AnalyzeEvent(*inputEvent);
        if (save) outputTree->Fill();
    }

    histEffMuonMomentum->Divide(histRecoMuonMomentum,
                                histPrimMuonMomentum,1.0,1.0,"B");
    histEffMuonLength->Divide(histRecoMuonLength,
                              histPrimMuonLength,1.0,1.0,"B");
    histEffMuonDCos->Divide(histRecoMuonDCos,
                            histPrimMuonDCos,1.0,1.0,"B");
    histEffMuonMomDCos->Divide(histRecoMuonMomDCos,
                            histPrimMuonMomDCos,1.0,1.0,"B");
    histEffMuonDir->Divide(histRecoMuonDir,
                           histPrimMuonDir,1.0,1.0,"B");
    histEffProtonMomentum->Divide(histRecoProtonMomentum,
                                  histPrimProtonMomentum,1.0,1.0,"B");
    histEffProtonDir->Divide(histRecoProtonDir,
                             histPrimProtonDir,1.0,1.0,"B");
    histEffPionMomentum->Divide(histRecoPionMomentum,
                                histPrimPionMomentum,1.0,1.0,"B");
    histEffPionDir->Divide(histRecoPionDir,
                           histPrimPionDir,1.0,1.0,"B");

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
