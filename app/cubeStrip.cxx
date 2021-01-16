#include "CubeEvent.hxx"
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

int main(int argc, char **argv) {
    std::cout << "CubeStrip: Hello World" << std::endl;
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
            std::cout << "Strip reconstruction information from an event. This"
                      << std::endl
                      << "   returns the file to the state it was in after"
                      << std::endl
                      << "   cubeERepTranslate."
                      << std::endl;
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
                      << "-s <number>  : Skip <number> entries"
                      << std::endl
                      << "   "
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl
                      << "   "
                      << "-h           : Do not strip the 3D hits"
                      << " [default: strip the hits]"
                      << std::endl;
            exit(1);
        }
        }
    }

    if (argc <= optind) {
        std::cout << "Use the -h option for help" << std::endl;
        throw std::runtime_error("Missing input file");
    }
    std::string inputName(argv[optind++]);
    std::cout << "Input Name " << inputName << std::endl;

    if (argc <= optind) {
        std::cout << "Use the -h option for help" << std::endl;
        throw std::runtime_error("Missing output file");
    }
    std::string outputName(argv[optind++]);
    std::cout << "Output Name " << outputName << std::endl;

    std::unique_ptr<TFile> inputFile(new TFile(inputName.c_str(),"old"));
    if (!inputFile->IsOpen()) {
        std::cout << "Use the -h option for help" << std::endl;
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
        outputEvent = inputEvent;
        outputEvent->MakeCurrentEvent();

        CUBE_LOG(0) << "Process event " << outputEvent->GetRunId()
                    << "/" << outputEvent->GetEventId() << std::endl;

        outputEvent->GetResultsContainer().clear();
        outputEvent->GetObjectContainers().clear();
        if (outputEvent->GetHitSelections().size() > 1) {
            outputEvent->GetHitSelections().erase(
                outputEvent->GetHitSelections().begin()+1,
                outputEvent->GetHitSelections().end());
        }
        for (Cube::Event::HitSelections::iterator s
                 = outputEvent->GetHitSelections().begin();
             s != outputEvent->GetHitSelections().end(); ++s) {
            std::cout << "Hit Selection " << (*s)->GetName() << std::endl;
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
