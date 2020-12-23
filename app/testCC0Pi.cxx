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

#include <iostream>
#include <sstream>
#include <map>
#include <set>

bool histInitialized = false;
TH1F* histMuonRes = NULL;
TH1F* histEarlyRes = NULL;
TH1F* histTrackRes = NULL;
TH1F* histClusterRes = NULL;
TH1F* histEarlyRMS = NULL;
TH1F* histTimeAll = NULL;
TH1F* histTimeTrack = NULL;
TH1F* histTimeCluster = NULL;
TH1F* histTimeNeutron = NULL;
TH1F* histTimeOther = NULL;
TH1F* histAvgCharge = NULL;
TH1F* histEarlyAvgCharge = NULL;
TH1F* histNeutronAvgCharge = NULL;
TH1F* histOtherAvgCharge = NULL;

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
bool AnalyzeEvent(Cube::Event& event) {
    /// Flag tht this event should be saved (if not rejected for some cut).
    bool shouldSave = false;

    if (!histInitialized) {
        histInitialized = true;
        histEarlyRes = new TH1F("earlyRes",
                                "Time resolution of earliest object",
                                200, 80.0, 120.0);
        histMuonRes = new TH1F("muonRes",
                                "Time resolution of muon object",
                                200, 80.0, 120.0);
        histTrackRes = new TH1F("trackRes",
                                "Time resolution of earliest object (track)",
                                200, 80.0, 120.0);
        histClusterRes = new TH1F("clusterRes",
                                "Time resolution of earliest object (cluster)",
                                200, 80.0, 120.0);
        histEarlyRMS = new TH1F("earlyRMS",
                                "Time RMS of earliest object",
                                50, 0.0, 10.0);
        histTimeAll = new TH1F("timeAll", "Muon to earliest activity",
                               100, -10.0, 10.0);
        histTimeTrack
            = new TH1F("timeTrack", "Muon to earliest activity (tracks)",
                       100, -10.0, 10.0);
        histTimeCluster
            = new TH1F("timeCluster", "Muon to earliest activity (clusters)",
                       100, -10.0, 10.0);
        histTimeNeutron = new TH1F("timeNeutron", "Muon to neutron activity",
                                100, -10.0, 10.0);
        histTimeOther = new TH1F("timeOther", "Muon to other activity",
                                100, -10.0, 10.0);
        histAvgCharge = new TH1F("avgCharge",
                                 "Avg hit charge of earliest activity",
                                150, 0.0, 1500.0);
        histEarlyAvgCharge
            = new TH1F("earlyAvgCharge",
                       "Early avg hit charge of earliest activity",
                       150, 0.0, 1500.0);
        histNeutronAvgCharge
            = new TH1F("neutronAvgCharge",
                       "Neutron avg hit charge of earliest activity",
                       150, 0.0, 1500.0);
        histOtherAvgCharge
            = new TH1F("otherAvgCharge",
                       "Other avg hit charge of earliest activity",
                       150, 0.0, 1500.0);
    }

    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;

    std::vector<int> muons;
    std::vector<int> antiMuons;
    std::vector<int> electrons;
    std::vector<int> gammas;
    std::vector<int> pions;
    std::vector<int> neutrons;
    std::vector<int> protons;
    std::vector<int> others;

    // Get the primaries for this event.
    std::vector<int> primId = Cube::Tool::AllPrimaries(event);
    for (std::vector<int>::iterator tr = primId.begin();
         tr != primId.end(); ++tr) {
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[*tr];
        if (!traj) {
            std::cout << "This cannot happen!";
            throw std::runtime_error("Invalid primary id");
        }
        switch (traj->GetPDGCode()) {
        case 11: case -11: electrons.push_back(*tr); break;
        case 13:  muons.push_back(*tr); break;
        case -13: antiMuons.push_back(*tr); break;
        case 22: gammas.push_back(*tr); break;
        case 111: case 211: case -211: pions.push_back(*tr); break;
        case 2112: neutrons.push_back(*tr); break;
        case 2212: protons.push_back(*tr); break;
        default:
            others.push_back(*tr);
            break;
        }
    }

    if (!electrons.empty()) return false;
    if (!muons.empty()) return false;
    if (antiMuons.empty()) return false;
    if (!pions.empty()) return false;
    if (!protons.empty()) return false;
    if (neutrons.size() != 1) return false;
    if (!gammas.empty()) return false;
    if (!others.empty()) return false;

    std::cout << "e: " << electrons.size()
              << " mu-: " << muons.size()
              << " mu+: " << antiMuons.size()
              << " p: " << protons.size()
              << " n: " << neutrons.size()
              << " pi: " << pions.size()
              << " g: " << gammas.size()
              << " o: " << others.size()
              << std::endl;

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return false;

    // Make sure there's a single muon track.
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) {
            continue;
        }
        int mainTraj = Cube::Tool::MainTrajectory(event,*track);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (!traj) {
            continue;
        }
        if (std::abs(traj->GetPDGCode()) != 13) {
            continue;
        }
        ++muonLike;
        muonObject = track;
    }

    if (muonLike > 1) return false;
    if (!muonObject) {
        return false;
    }

    // Collect the earliest objects and the muon.
    double muonTime = muonObject->GetPosition().T();
    double earliestTime = 1E+8;
    Cube::Handle<Cube::ReconObject> earliestObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject,*(*o))) continue;
        Cube::Handle<Cube::ReconTrack> track = *o;
        double objTime = -1.0;
        if (track) {
            if (track == muonObject) continue;
            objTime = track->GetPosition().T();
        }
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (cluster) {
            objTime = cluster->GetPosition().T();
        }
        if (objTime < 0) continue;
        if (objTime < earliestTime) {
            earliestTime = objTime;
            earliestObject = *o;
        }

#ifdef DUMP_STACK
        int mainTraj = Cube::Tool::MainTrajectory(event,*(*o));
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        int parentTraj = mainTraj;
        std::cout << "Object stack: "
        while (parentTraj >= 0) {
            traj = trajectories[parentTraj];
            std::cout << parentTraj << "(" << traj->GetPDGCode()
                       << "," << traj->GetInitialPosition().T()
                       << "," << (traj->GetInitialMomentum().E()
                                  - traj->GetInitialMomentum().M()) << ") ";
            parentTraj = traj->GetParentId();
        }
        std::cout << std::endl;
#endif
    }
    if (!earliestObject) return false;

    std::vector<Cube::Handle<Cube::G4Hit>> earliestSegs
        = Cube::Tool::ObjectG4Hits(event,*earliestObject);
    double earliestTruth = 1E+8;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
             t = earliestSegs.begin();
         t != earliestSegs.end(); ++t) {
        if ((*t)->GetStart().T() < earliestTruth) {
            earliestTruth = (*t)->GetStart().T();
        }
    }

    std::vector<Cube::Handle<Cube::G4Hit>> muonSegs
        = Cube::Tool::ObjectG4Hits(event,*muonObject);
    double muonTruth = 1E+8;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
             t = muonSegs.begin();
         t != muonSegs.end(); ++t) {
        if ((*t)->GetStart().T() < muonTruth) {
            muonTruth = (*t)->GetStart().T();
        }
    }

    Cube::Handle<Cube::HitSelection> hits = earliestObject->GetHitSelection();
    if (!hits) return false;

    double totalCharge = 0.0;
    double timeAvg = 0.0;
    double time2Avg = 0.0;
    double latestTime = -1.0;
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        totalCharge += (*h)->GetCharge();
        timeAvg += (*h)->GetCharge()*(*h)->GetTime();
        time2Avg += (*h)->GetCharge()*(*h)->GetTime()*(*h)->GetTime();
        if ((*h)->GetTime() > latestTime) {
            latestTime = (*h)->GetTime();
        }
    }
    timeAvg = timeAvg/totalCharge;
    time2Avg = time2Avg/totalCharge;
    double timeRMS = std::sqrt(time2Avg-timeAvg*timeAvg);
    double avgCharge = totalCharge/hits->size();

    histTimeAll->Fill(earliestTime - muonTime);
    Cube::Handle<Cube::ReconTrack> track = earliestObject;
    if (track) {
        histTimeTrack->Fill(earliestTime - muonTime);
        histTrackRes->Fill(earliestTime - earliestTruth);
        if (earliestTime < muonTime) shouldSave = true;
    }
    else {
        histTimeCluster->Fill(earliestTime - muonTime);
        histClusterRes->Fill(earliestTime - earliestTruth);
    }
    histMuonRes->Fill(muonTime - muonTruth);
    histEarlyRes->Fill(earliestTime - earliestTruth);
    histEarlyRMS->Fill(timeRMS);

    histAvgCharge->Fill(avgCharge);
    if (earliestTime < muonTime) {
        histEarlyAvgCharge->Fill(avgCharge);
    }

    int earliestTraj = Cube::Tool::MainTrajectory(event,*earliestObject);
    int earliestPrim = Cube::Tool::PrimaryId(event,earliestTraj);
    if (earliestPrim == 1) {
        histTimeNeutron->Fill(earliestTime - muonTime);
        histNeutronAvgCharge->Fill(avgCharge);
    }
    else {
        histTimeOther->Fill(earliestTime - muonTime);
        histOtherAvgCharge->Fill(avgCharge);
    }

    std::cout << "TTT " << muonTime
              << " " << earliestTime
              << " " << avgCharge
              << " " << totalCharge
              << std::endl;

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
