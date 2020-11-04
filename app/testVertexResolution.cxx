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

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    if (!histInitialized) {
        histInitialized = true;
        histVertexResolution = new TH1F("vertexResolution","Vertex Resolution",
                                        50,0.0,50.0);
        histVertexXResidual = new TH1F("vertexXResidual",
                                       "Vertex residual on X axis",
                                       100,-50.0,50.0);
        histVertexYResidual = new TH1F("vertexYResidual",
                                       "Vertex residual on Y axis",
                                       100,-50.0,50.0);
        histVertexZResidual = new TH1F("vertexZResidual",
                                       "Vertex residual on Z axis",
                                       100,-50.0,50.0);
        histVertexTResidual = new TH1F("vertexTResidual",
                                       "Vertex residual in time",
                                       100,-10.0,10.0);

        histTrackResolution = new TH1F("trackResolution","Track Resolution",
                                        50,0.0,50.0);
        histTrackXResidual = new TH1F("trackXResidual",
                                       "Track residual on X axis",
                                       100,-50.0,50.0);
        histTrackYResidual = new TH1F("trackYResidual",
                                       "Track residual on Y axis",
                                       100,-50.0,50.0);
        histTrackZResidual = new TH1F("trackZResidual",
                                       "Track residual on Z axis",
                                       100,-50.0,50.0);
        histTrackTResidual = new TH1F("trackTResidual",
                                       "Track residual in time",
                                       100,-10.0,10.0);
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
        if (frontDist < backDist) {
            offset = bestTrack->GetFront()->GetPosition() - primaryVertex;
        }
        else {
            offset = bestTrack->GetBack()->GetPosition() - primaryVertex;
        }
        histTrackResolution->Fill(offset.Vect().Mag());
        histTrackXResidual->Fill(offset.X());
        histTrackYResidual->Fill(offset.Y());
        histTrackZResidual->Fill(offset.Z());
        histTrackTResidual->Fill(offset.T());
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
