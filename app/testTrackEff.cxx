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
TH1F* histRecoMuonMomentum = NULL;
TH1F* histEffMuonMomentum = NULL;
TH2F* histPrimMuonDir = NULL;
TH2F* histRecoMuonDir = NULL;
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

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    if (!histInitialized) {
        histInitialized = true;
        histPrimMuonMomentum = new TH1F("primMuonMomentum",
                                       "True Momentum of Muon",
                                       50,0.0,500.0);
        histRecoMuonMomentum = new TH1F("recoMuonMomentum",
                                        "True Momentum of Reconstructed Muon",
                                        50,0.0,500.0);
        histEffMuonMomentum = new TH1F("effMuonMomentum",
                                        "Efficiency at True Muon Momentum",
                                        50,0.0,500.0);
        histPrimMuonDir = new TH2F("primMuonDir",
                                   "True Cosine vs Momentum for Muons",
                                   10,0.0,500.0,
                                   10,0.0,1.0);
        histRecoMuonDir = new TH2F("recoMuonDir",
                                   "True Cosine vs Momentum for Reco Muons",
                                   10,0.0,500.0,
                                   10,0.0,1.0);
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

    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    for (Cube::Event::G4TrajectoryContainer::iterator tr = trajectories.begin();
         tr != trajectories.end(); ++tr) {
        // only count primaries.
        if (tr->second->GetParentId() >= 0) continue;
        double fid= Cube::Tool::ContainedPoint(
            tr->second->GetInitialPosition().Vect());
        int pdgMuon = 13;
        if (tr->second->GetPDGCode() == pdgMuon && fid > 50.0) {
            double mom = tr->second->GetInitialMomentum().P();
            histPrimMuonMomentum->Fill(mom);
            TVector3 dir = tr->second->GetInitialMomentum().Vect().Unit();
            histPrimMuonDir->Fill(mom,std::abs(dir.Z()));
        }
        int pdgProton = 2212;
        if (tr->second->GetPDGCode() == pdgProton && fid > 50.0) {
            double mom = tr->second->GetInitialMomentum().P();
            histPrimProtonMomentum->Fill(mom);
            TVector3 dir = tr->second->GetInitialMomentum().Vect().Unit();
            histPrimProtonDir->Fill(mom,std::abs(dir.Z()));
        }
        int pdgPion = 211;
        if (tr->second->GetPDGCode() == pdgPion && fid > 50.0) {
            double mom = tr->second->GetInitialMomentum().P();
            histPrimPionMomentum->Fill(mom);
            TVector3 dir = tr->second->GetInitialMomentum().Vect().Unit();
            histPrimPionDir->Fill(mom,std::abs(dir.Z()));
        }
    }

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return;

    int pdgMuon = 13;
    int pdgProton = 2212;
    int pdgPion = 211;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) continue;
        int mainTraj = Cube::Tool::MainTrajectory(event,*track);
        if (mainTraj<0) continue;
        int primTraj = Cube::Tool::PrimaryId(event,mainTraj);
        if (primTraj != mainTraj) continue;
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[primTraj];
        if (!traj) throw;
        double fid= Cube::Tool::ContainedPoint(
            traj->GetInitialPosition().Vect());
        if (traj->GetPDGCode() == pdgMuon && fid > 50.0) {
            double mom = traj->GetInitialMomentum().P();
            histRecoMuonMomentum->Fill(mom);
            TVector3 dir = traj->GetInitialMomentum().Vect().Unit();
            histRecoMuonDir->Fill(mom,std::abs(dir.Z()));
            pdgMuon = -1;
        }
        if (traj->GetPDGCode() == pdgProton && fid > 50.0) {
            double mom = traj->GetInitialMomentum().P();
            histRecoProtonMomentum->Fill(mom);
            TVector3 dir = traj->GetInitialMomentum().Vect().Unit();
            histRecoProtonDir->Fill(mom,std::abs(dir.Z()));
            pdgProton = -1;
        }
        if (traj->GetPDGCode() == pdgPion && fid > 50.0) {
            double mom = traj->GetInitialMomentum().P();
            histRecoPionMomentum->Fill(mom);
            TVector3 dir = traj->GetInitialMomentum().Vect().Unit();
            histRecoPionDir->Fill(mom,std::abs(dir.Z()));
            pdgPion = -1;
        }
    }
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

    // Loop through the events.
    int totalEntries = inputChain->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputChain->GetEntry(entry);
        std::cout << "Process event "
                  << entry
                  << "/" << inputEvent->GetRunId()
                  << "/" << inputEvent->GetEventId() << std::endl;
        AnalyzeEvent(*inputEvent);
        for (Cube::Event::G4HitContainer::iterator g
                 = inputEvent->G4Hits.begin();
             g != inputEvent->G4Hits.end(); ++g) {
        }
    }

    histEffMuonMomentum->Divide(histRecoMuonMomentum,
                                histPrimMuonMomentum,1.0,1.0,"B");
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
