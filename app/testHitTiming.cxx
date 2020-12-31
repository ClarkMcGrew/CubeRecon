#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeEvent.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolCubeTruth.hxx>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    static TH1F* histHitTiming = NULL;
    static TH1F* histAvgHitTiming = NULL;
    static TH1F* histMaxHitTiming = NULL;
    static TH1F* histLateHitTiming = NULL;
    static TH1F* histLateHitCount = NULL;
    static TH1F* histLateHitDiff = NULL;
    static TH1F* histWghtHitTiming = NULL;
    static TH1F* histCubeCount = NULL;

    if (!histHitTiming) {
        std::cout << "Create the histogram" << std::endl;
        histHitTiming = new TH1F("HitTiming","Time Resolution for Default",
                                 80,-20.0,20.0);
        histAvgHitTiming = new TH1F("AvgHitTiming","Time Resolution for Avg.",
                                    80,-20.0,20.0);
        histMaxHitTiming = new TH1F("MaxHitTiming",
                                    "Time Resolution for Last Hit",
                                    80,-20.0,20.0);
        histLateHitTiming = new TH1F("LateHitTiming",
                                     "Time Resolution for Avg. Near Last Hit",
                                     80,-20.0,20.0);
        histLateHitCount = new TH1F("LateHitCount",
                                    "Number of fibers near last hit",
                                    4,0.0,4.0);
        histLateHitDiff = new TH1F("LateHitDiff",
                                   "Time from last fiber",
                                   80,-20.0,0.0);
        histWghtHitTiming = new TH1F("wghtHitTiming",
                                     "Time Resolution for weighting",
                                    80,-20.0,20.0);
        histCubeCount = new TH1F("cubeCount",
                                 "Number of cubes with energy on a fiber",
                                 10,0.0, 10.0);
    }

    Cube::Handle<Cube::HitSelection> eventHits
        = event.GetHitSelection();
    if (!eventHits) return;

    // Build a map of 2D to 3D hits.
    std::map<Cube::Handle<Cube::Hit>,std::set<Cube::Handle<Cube::Hit>>> D2ToD3;
    for (Cube::HitSelection::iterator h = eventHits->begin();
         h != eventHits->end(); ++h) {
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            D2ToD3[ch].insert(*h);
        }
    }

    for (Cube::HitSelection::iterator h = eventHits->begin();
         h != eventHits->end(); ++h) {
        if (!(*h)) continue;
        double energy = Cube::Tool::CubeDeposit(event,*h);
        double xtalk =  Cube::Tool::CubeCrossTalk(event,*h);
        double trueTime =  Cube::Tool::CubeTime(event,*h);
        histHitTiming->Fill((*h)->GetTime()-trueTime-100.0);
        double avgT = 0.0;
        double avgW = 0.0;
        double maxT = -1E+20;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - dd/200.0;
            avgT += tt*ch->GetCharge();
            avgW += ch->GetCharge();
            maxT = std::max(tt,maxT);
        }
        avgT = avgT/avgW;
        histMaxHitTiming->Fill(maxT-trueTime-100.0);
        histAvgHitTiming->Fill(avgT-trueTime-100.0);

        /////////////////////////////////////////////////////////////
        // Calculate the hit time by looking at the latest hit after
        // correcting for distance.  Since the fiber time usually comes from
        // the cube hit on the fiber, this tends to make second and third
        // cubes on the fiber have an early estimated time.  The hit times
        // have a small early bias.
        double lateT = 0.0;
        double lateW = 0.0;
        if ((*h)->GetConstituentCount() != 3) {
            std::cout <<" OH! " << (*h)->GetConstituentCount() << std::endl;
        }
        int lateCount = 0;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            if (ch->GetContributorCount() < 1) {
                std::cout << " Charge " << (ch)->GetCharge() << std::endl;
                continue;
            }
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - dd/200.0;
            if (tt < maxT - 1E-10) histLateHitDiff->Fill(tt-maxT);
            if (tt < maxT - 2.5) continue;
            ++lateCount;
            lateT += tt*ch->GetCharge();
            lateW += ch->GetCharge();
        }
        lateT = lateT/lateW;
        histLateHitTiming->Fill(lateT-trueTime-100.0);
        histLateHitCount->Fill(lateCount+0.5);

        for (auto hh:  D2ToD3) {
            histCubeCount->Fill(hh.second.size()+0.5);
        }

        ////////////////////////////////////////////////////////////
        // Calculate hit time assuming that all deposites along the fiber are
        // close to "simultaneous".  This means that the time for the fiber is
        // corrected to the distance of the closest crossing fiber.  This
        // tends to make the hit time to be bias a bit late.

        // Find the latest time for the entire cube.
        maxT = 0.0;
        double lateFiberT = 1E+8;
        double lateFiberD = 1E+8;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            double minD = 1E+8;
            double maxD = 0;
            for (auto hh:  D2ToD3[ch]) {
                double dd = (hh->GetPosition() - ch->GetPosition()).Mag();
                if (dd < minD) minD = dd;
                if (dd > maxD) maxD = dd;
            }
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - minD/200.0;
            if (maxT < tt) maxT = tt;
            for (int s = 0; s < ch->GetContributorCount(); ++s) {
                int seg = ch->GetContributor(s);
                Cube::Handle<Cube::G4Hit> g4Hit = event.G4Hits[seg];
                if (!g4Hit) continue;
                double tf = g4Hit->GetStart().T();
                if (tf < lateFiberT) {
                    lateFiberT = tf;
                    lateFiberD = (g4Hit->GetStart().Vect()
                                  - ch->GetPosition()).Mag();
                }
            }
        }
        // Now take the average of times close to the maxT.
        lateT = 0.0;
        lateW = 0.0;
        double lateQ = 0.0;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            double minD = 1E+8;
            double maxD = 0;
            for (auto hh:  D2ToD3[ch]) {
                double dd = (hh->GetPosition() - ch->GetPosition()).Mag();
                if (dd < minD) minD = dd;
                if (dd > maxD) maxD = dd;
            }
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - minD/200.0;
            if (tt < maxT - 2.5) continue;
            lateT += tt*ch->GetCharge();
            lateW += ch->GetCharge();
        }
        lateT = lateT/lateW;
        histWghtHitTiming->Fill(lateT-trueTime-100.0);
    }

}

int main(int argc, char** argv) {
    int maxEntries = 1E+8; // Maximum to process.
    int firstEntry = 0;

    while (true) {
        int c = getopt(argc,argv,"n:s:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
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
                      << "-s <number>  : Skip <number> entries"
                      << std::endl
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl;
            exit(1);
        }
        }
    }

    if (argc <= optind) throw std::runtime_error("Missing input file");
    std::string inputName(argv[optind++]);
    std::cout << "Input Name " << inputName << std::endl;

    std::string outputName;
    if (argc > optind) {
        outputName = argv[optind++];
    }
    else {
        std::cout << "NO OUTPUT FILE!!!!" << std::endl;
    }

    // Attach to the input tree.
    std::unique_ptr<TFile> inputFile(new TFile(inputName.c_str(),"old"));
    if (!inputFile->IsOpen()) throw std::runtime_error("Input file not open");

    /// Attach to the input tree.
    TTree* inputTree = (TTree*) inputFile->Get("CubeEvents");
    if (!inputTree) throw std::runtime_error("Missing the event tree");
    Cube::Event *inputEvent = NULL;
    inputTree->SetBranchAddress("Event",&inputEvent);

    // Open the output file
    std::unique_ptr<TFile> outputFile;
    if (!outputName.empty()) {
        std::cout << "Open Output File: " << outputName << std::endl;
        outputFile.reset(new TFile(outputName.c_str(),"recreate"));
    }

    // Loop through the events.
    int totalEntries = inputTree->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputTree->GetEntry(entry);
        std::cout << "Process event " << inputEvent->GetRunId()
                    << "/" << inputEvent->GetEventId() << std::endl;
        AnalyzeEvent(*inputEvent);
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
