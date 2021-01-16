#include "CubeEvent.hxx"
#include "CubeHit.hxx"
#include "CubeHitSelection.hxx"
#include "CubeAlgorithmResult.hxx"

#include "CubeInfo.hxx"
#include "CubeERepSim.hxx"
#include "ERepSimInput.hxx"

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
    std::cout << "CubeRecon: Hello World" << std::endl;
    int maxEntries = 1E+8; // Maximum to process.

    while (true) {
        int c = getopt(argc,argv,"n:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
            break;
        }
        default: {
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
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

    ERepSim::Input::Get().Attach(inputFile.get());
    int totalEntries = ERepSim::Input::Get().DataTree->GetEntries();

    Cube::Event event;
    std::unique_ptr<TFile> outputFile(new TFile(outputName.c_str(),"recreate"));
    TTree *dataTree = new TTree("CubeEvents","Reconstructed Event");
    static Cube::Event *pEvent = &event;
    dataTree->Branch("Event",&pEvent);

    TGeoManager* geom
        = dynamic_cast<TGeoManager*>(inputFile->Get("EDepSimGeometry"));
    if (!geom) {
        geom = dynamic_cast<TGeoManager*>(inputFile->Get("CubeReconGeometry"));
    }

    if (geom && outputFile) {
        geom->SetName("CubeReconGeometry");
        geom->Write();
    }

    totalEntries = std::min(totalEntries,maxEntries);
    for (int entry = 0; entry < totalEntries; ++entry) {
        ERepSim::Input::Get().GetEntry(entry);

        Cube::ConvertERepSim(event);

        dataTree->Fill();
    }

    outputFile->Write();

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
