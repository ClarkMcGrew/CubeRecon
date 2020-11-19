#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconVertex.hxx>
#include <CubeReconTrack.hxx>
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
#include <TPad.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

bool histInitialized = false;
TH1F* histVertexResolution = NULL;
TH1F* histVertexXResidual = NULL;
TH1F* histVertexYResidual = NULL;
TH1F* histVertexZResidual = NULL;
TH1F* histVertexTResidual = NULL;

TH1F* histTrackResolution = NULL;
TH1F* histTrackXResidual = NULL;
TH1F* histTrackYResidual = NULL;
TH1F* histTrackZResidual = NULL;
TH1F* histTrackTResidual = NULL;
TH1F* histTrackPerp = NULL;
TH1F* histTrackLong = NULL;

void PlotHistogram(TH1* hist) {
    std::string name(hist->GetName());
    std::cout << name + ".C" << std::endl;
    std::cout << name + ".png" << std::endl;
    std::cout << name + ".pdf" << std::endl;
    hist->Draw();
    gPad->Update();
    gPad->Print((name + ".png").c_str());
    gPad->Print((name + ".pdf").c_str());
    gPad->Print((name + ".C").c_str());
}

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    if (!histInitialized) {
        histInitialized = true;
        histVertexResolution = new TH1F("vertexResolution",
                                        "Vertex Resolution (multi-track)",
                                        50,0.0,50.0);
        histVertexResolution->SetXTitle("Residual (mm)");

        histVertexXResidual = new TH1F(
            "vertexXResidual",
            "Vertex Residual (multi-track) on X-Axis",
            100,-50.0,50.0);
        histVertexXResidual->SetXTitle("Residual (mm)");

        histVertexYResidual = new TH1F(
            "vertexYResidual",
            "Vertex Residual (multi-track) on Y-Axis",
            100,-50.0,50.0);
        histVertexYResidual->SetXTitle("Residual (mm)");

        histVertexZResidual = new TH1F(
            "vertexZResidual",
            "Vertex Residual (multi-track) on Z-Axis",
            100,-50.0,50.0);
        histVertexZResidual->SetXTitle("Residual (mm)");

        histVertexTResidual = new TH1F(
            "vertexTResidual",
            "Vertex Residual (multi-track) in Time",
            100,-10.0,10.0);
        histVertexTResidual->SetXTitle("Residual (ns)");

        histTrackResolution = new TH1F("trackResolution",
                                       "Vertex Resolution (single-track)",
                                        50,0.0,50.0);
        histTrackResolution->SetXTitle("Residual (mm)");

        histTrackXResidual = new TH1F(
            "trackXResidual",
            "Vertex Residual (single-track) on X-Axis",
            100,-50.0,50.0);
        histTrackXResidual->SetXTitle("Residual (mm)");

        histTrackYResidual = new TH1F(
            "trackYResidual",
            "Vertex Residual (single-track) on Y-Axis",
            100,-50.0,50.0);
        histTrackYResidual->SetXTitle("Residual (mm)");

        histTrackZResidual = new TH1F(
            "trackZResidual",
            "Vertex Residual (single-track) on Z-Axis",
            100,-50.0,50.0);
        histTrackZResidual->SetXTitle("Residual (mm)");

        histTrackTResidual = new TH1F(
            "trackTResidual",
            "Vertex Residual (single-track) in Time",
            100,-10.0,10.0);
        histTrackTResidual->SetXTitle("Residual (ns)");

        histTrackPerp = new TH1F("trackPerpResolution",
                                 "Vertex (single-track) Impact Parameter",
                                 50,0.0,50.0);
        histTrackPerp->SetXTitle("Residual (mm)");

        histTrackLong = new TH1F("trackLongResolution",
                                 "Vertex (single-track) Longitudinal Residual",
                                 100,-50.0,50.0);
        histTrackLong->SetXTitle("Residual (mm)");
    }

    // Find the primary vertex.
    TLorentzVector primaryVertex;
    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    for (Cube::Event::G4TrajectoryContainer::iterator tr = trajectories.begin();
         tr != trajectories.end(); ++tr) {
        // only count primaries.
        if (tr->second->GetParentId() >= 0) continue;
        primaryVertex = tr->second->GetInitialPosition();
        break;
    }
    primaryVertex.SetT(primaryVertex.T() + 100*unit::ns);

    // Only look at contained primary vertices.
    double fid= Cube::Tool::ContainedPoint(primaryVertex.Vect());
    if (fid < 0.0) return;

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return;

    Cube::Handle<Cube::ReconVertex> bestVertex;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconVertex> vertex = *o;
        if (!vertex) continue;
        double newDist =  (vertex->GetPosition().Vect()
                           -primaryVertex.Vect()).Mag();
        if (newDist > 10.0*unit::cm) continue;
        if (!bestVertex) {
            bestVertex = vertex;
            continue;
        }
        double bestDist = (bestVertex->GetPosition().Vect()
                           -primaryVertex.Vect()).Mag();
        if (newDist > bestDist) continue;
        bestVertex = vertex;
    }

    Cube::Handle<Cube::ReconTrack> bestTrack;
    double bestTrackDist = 1E+20;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) continue;
        double frontDist =  (track->GetFront()->GetPosition().Vect()
                             -primaryVertex.Vect()).Mag();
        double backDist =  (track->GetBack()->GetPosition().Vect()
                            -primaryVertex.Vect()).Mag();
        double newDist = std::min(frontDist,backDist);
        if (newDist > 10.0*unit::cm) continue;
        if (!bestTrack) {
            bestTrack = track;
            continue;
        }
        if (newDist > bestTrackDist) continue;
        bestTrack = track;
    }

    if (!bestVertex && bestTrack) {
        double frontDist =  (bestTrack->GetFront()->GetPosition().Vect()
                             -primaryVertex.Vect()).Mag();
        double backDist =  (bestTrack->GetBack()->GetPosition().Vect()
                            -primaryVertex.Vect()).Mag();
        TLorentzVector offset;
        double perp = 1E+6;
        double para = 1E+6;
        if (frontDist < backDist) {
            offset = bestTrack->GetFront()->GetPosition() - primaryVertex;
            para = offset.Vect()*bestTrack->GetFront()->GetDirection();
            perp = (offset.Vect()
                    - para*bestTrack->GetFront()->GetDirection()).Mag();
            // Correct for the sense of the direction
            para = -para;
        }
        else {
            offset = bestTrack->GetBack()->GetPosition() - primaryVertex;
            para = offset.Vect()*bestTrack->GetBack()->GetDirection();
            perp = (offset.Vect()
                    - para*bestTrack->GetBack()->GetDirection()).Mag();
        }
        histTrackResolution->Fill(offset.Vect().Mag());
        histTrackXResidual->Fill(offset.X());
        histTrackYResidual->Fill(offset.Y());
        histTrackZResidual->Fill(offset.Z());
        histTrackTResidual->Fill(offset.T());
        histTrackPerp->Fill(perp);
        histTrackLong->Fill(para);

        return;
    }

    if (!bestVertex) return;

    TLorentzVector offset = bestVertex->GetPosition() - primaryVertex;
    if (offset.Vect().Mag() > 10.0*unit::cm) return;

    histVertexResolution->Fill(offset.Vect().Mag());
    histVertexXResidual->Fill(offset.X());
    histVertexYResidual->Fill(offset.Y());
    histVertexZResidual->Fill(offset.Z());
    histVertexTResidual->Fill(offset.T());

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

    PlotHistogram(histVertexResolution);
    PlotHistogram(histVertexXResidual);
    PlotHistogram(histVertexYResidual);
    PlotHistogram(histVertexZResidual);
    PlotHistogram(histVertexTResidual);
    PlotHistogram(histTrackResolution);
    PlotHistogram(histTrackXResidual);
    PlotHistogram(histTrackYResidual);
    PlotHistogram(histTrackZResidual);
    PlotHistogram(histTrackTResidual);
    PlotHistogram(histTrackPerp);
    PlotHistogram(histTrackLong);

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
