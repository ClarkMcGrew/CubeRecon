#include "CubeEvent.hxx"
#include "CubeMakeHits3D.hxx"
#include "CubeRecon.hxx"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVector.h>

#include <iostream>
#include <sstream>
#include <exception>
#include <memory>

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
    static Cube::Event *outputEvent = inputEvent;
    outputTree->Branch("Event",&outputEvent);

    // Loop through the events.
    int totalEntries = inputTree->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputTree->GetEntry(entry);
        outputEvent = inputEvent;
        outputEvent->MakeCurrentEvent();

        CUBE_LOG(0) << "Process event " << outputEvent->GetRunId()
                    << "/" << outputEvent->GetEventId() << std::endl;

        // Check if MakeHits3D has been run.  If it is missing, then run it.
        // This will leave the 3D hits as the main hit selection for the
        // event.
        Cube::Handle<Cube::AlgorithmResult> makeHits3D
            = outputEvent->GetAlgorithmResult("MakeHits3D");
        if (!makeHits3D) {
            std::unique_ptr<Cube::MakeHits3D>
                algoMakeHits3D(new Cube::MakeHits3D);
            makeHits3D = algoMakeHits3D->Process(*outputEvent);
            outputEvent->AddAlgorithmResult(makeHits3D);
            outputEvent->AddHitSelection(makeHits3D->GetHitSelection());
        }

        // Get the main hits for the event.
        Cube::Handle<Cube::HitSelection> hits3D
            = outputEvent->GetHitSelection();

        // Run the main reconstruction.
        if (hits3D) {
            std::unique_ptr<Cube::Recon>
                algoRecon(new Cube::Recon);
            Cube::Handle<Cube::AlgorithmResult>
                recon = algoRecon->Process(*hits3D);
            outputEvent->AddAlgorithmResult(recon);
            Cube::Handle<Cube::ReconObjectContainer> finalObjects
                = recon->GetObjectContainer("final");
            if (finalObjects) outputEvent->AddObjectContainer(finalObjects);
        }


        CUBE_LOG(0) << "Finished event " << outputEvent->GetRunId()
                    << "/" << outputEvent->GetEventId() << std::endl;
        outputTree->Fill();
    }

    outputFile->Write();

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
