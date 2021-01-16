#include "CubeEvent.hxx"
#include "CubeMakeHits3D.hxx"
#include "CubeRecon.hxx"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVector.h>
#include <TGeoManager.h>

#include <iostream>
#include <sstream>
#include <exception>
#include <memory>

Cube::Event *gOutputEvent = NULL;

int main(int argc, char **argv) {
    std::cout << "CubeRecon: Hello World" << std::endl;
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

    if (argc <= optind) {
        throw std::runtime_error("Missing input file");
    }
    std::string inputName(argv[optind++]);
    std::cout << "Input Name " << inputName << std::endl;

    if (argc <= optind) {
        throw std::runtime_error("Missing output file");
    }
    std::string outputName(argv[optind++]);
    std::cout << "Output Name " << outputName << std::endl;

    std::unique_ptr<TFile> inputFile(new TFile(inputName.c_str(),"old"));
    if (!inputFile->IsOpen()) {
        throw std::runtime_error("File not open");
    }
    std::cout << "Input File " << inputFile->GetName() << std::endl;

    /// Attach to the input tree.
    TTree* inputTree = (TTree*) inputFile->Get("CubeEvents");
    if (!inputTree) {
        std::cout << "Missing the event tree" << std::endl;
        return 1;
    }
    Cube::Event *inputEvent = NULL;
    inputTree->SetBranchAddress("Event",&inputEvent);

    // Attach to the output tree.
    std::unique_ptr<TFile> outputFile(new TFile(outputName.c_str(),"recreate"));
    TTree *outputTree = new TTree("CubeEvents","Reconstructed Event");
    gOutputEvent = inputEvent;
    outputTree->Branch("Event",&gOutputEvent);

    TGeoManager* geom
        = dynamic_cast<TGeoManager*>(inputFile->Get("CubeReconGeometry"));
    if (!geom) {
        std::cout << "No geometry for file" << std::endl;
    }

    if (geom && outputFile) {
        geom->SetName("CubeReconGeometry");
        geom->Write();
    }

    // Loop through the events.
    int totalEntries = inputTree->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputTree->GetEntry(entry);
        gOutputEvent = inputEvent;
        gOutputEvent->MakeCurrentEvent();

        CUBE_LOG(0) << "#####################################################"
                    << std::endl;
        CUBE_LOG(0) << "Process event "
                    << entry << "/" << gOutputEvent->GetRunId()
                    << "/" << gOutputEvent->GetEventId() << std::endl;
        CUBE_LOG(0) << "#####################################################"
                    << std::endl;

        // Check if MakeHits3D has been run.  If it is missing, then run it.
        // This will leave the 3D hits as the main hit selection for the
        // event.
        Cube::Handle<Cube::AlgorithmResult> makeHits3D
            = gOutputEvent->GetAlgorithmResult("MakeHits3D");
        if (!makeHits3D) {
            std::unique_ptr<Cube::MakeHits3D>
                algoMakeHits3D(new Cube::MakeHits3D);
            makeHits3D = algoMakeHits3D->Process(*gOutputEvent);
            gOutputEvent->AddAlgorithmResult(makeHits3D);
            gOutputEvent->AddHitSelection(makeHits3D->GetHitSelection());
        }

        // Get the main hits for the event.
        Cube::Handle<Cube::HitSelection> hits3D
            = gOutputEvent->GetHitSelection();

        // Run the main reconstruction.
        if (hits3D) {
            std::unique_ptr<Cube::Recon>
                algoRecon(new Cube::Recon);
            Cube::Handle<Cube::AlgorithmResult>
                recon = algoRecon->Process(*hits3D);
            gOutputEvent->AddAlgorithmResult(recon);
            Cube::Handle<Cube::ReconObjectContainer> finalObjects
                = recon->GetObjectContainer("final");
            if (finalObjects) gOutputEvent->AddObjectContainer(finalObjects);
        }


        CUBE_LOG(0) << "Finished event " << gOutputEvent->GetRunId()
                    << "/" << gOutputEvent->GetEventId() << std::endl;
        outputTree->Fill();
    }

    outputFile->Write();

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
