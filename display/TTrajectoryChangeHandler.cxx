#include <CubeUnits.hxx>

#include "TTrajectoryChangeHandler.hxx"
#include "TEventDisplay.hxx"
#include "TGUIManager.hxx"

#include "TEventManager.hxx"

#include <TGeoManager.h>
#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>
#include <TDatabasePDG.h>

#include <sstream>
#include <iostream>

Cube::TTrajectoryChangeHandler::TTrajectoryChangeHandler() {
    fTrajectoryList = new TEveElementList("g4Trajectories",
                                          "Geant4 Trajectories");
    fTrajectoryList->SetMainColor(kRed);
    fTrajectoryList->SetMainAlpha(1.0);
    gEve->AddElement(fTrajectoryList);
}

Cube::TTrajectoryChangeHandler::~TTrajectoryChangeHandler() {
}

void Cube::TTrajectoryChangeHandler::Apply() {
    const TDatabasePDG *pdg = TDatabasePDG::Instance();

    fTrajectoryList->DestroyElements();

    if (!Cube::TEventDisplay::Get().GUI().GetShowTrajectoriesButton()->IsOn()) {
        std::cout <<"Trajectories are off" << std::endl;
        return;
    }

    if (!Cube::gEvent) return;
    std::cout <<"Handle the trajectories" << std::endl;

    for (Cube::Event::G4TrajectoryContainer::iterator t
             = Cube::gEvent->G4Trajectories.begin();
         t != Cube::gEvent->G4Trajectories.end();
         ++t) {
        Cube::Handle<Cube::G4Trajectory> traj = t->second;

#ifndef INCLUDE_ROCK_TRAJECTORIES
        {
            TVector3 pos = traj->GetInitialPosition().Vect();
            TVector3 center(0.0, -3000.0*unit::mm, 63.0*unit::meter);
            if ((pos-center).Mag() > 5*unit::meter) continue;
        }
#endif
        double eKin = traj->GetInitialMomentum().E() -
            traj->GetInitialMomentum().M();
        if (eKin < 5.0) eKin = 5.0;

        std::string pdgName = "undefined";
        int charged = 0;

        const TParticlePDG* pdgParticle = pdg->GetParticle(traj->GetPDGCode());
        if (pdgParticle) {
            pdgName = pdgParticle->GetName();
            if (pdgParticle->Charge() > 0.0) charged = +1;
            if (pdgParticle->Charge() < 0.0) charged = -1;
        }
        else {
            std::cout << "######################################" << std::endl;
            std::cout << "# No Trajectory Code: " << traj->GetPDGCode()
                      << std::endl;
            std::cout << "######################################" << std::endl;
        }

        double beta
            = traj->GetInitialMomentum().P()/traj->GetInitialMomentum().E();

        double wght = 1.0/beta/beta;
        int iWght = 2*wght;
        if (iWght < 1) iWght = 1;

        std::ostringstream label;
        label << traj->GetTrackId() << " " << pdgName
              << " from " << traj->GetParentId()
              << " (" << eKin << " MeV)";

        if (traj->GetParentId() < 0) {
            std::cout << "Primary " << traj->GetTrackId()
                      << " " << pdgName
                      << " " << eKin << " MeV"
                      << " " << beta
                      << " " << wght << " " << iWght
                      << " " << charged << std::endl;
        }

        TEveLine *track = new TEveLine();

        track->SetName("trajectory");
        track->SetTitle(label.str().c_str());

        if (charged<0) {
            track->SetLineColor(kAzure-4);
            track->SetLineWidth(iWght);
        }
        else if (charged>0) {
            track->SetLineColor(kCyan);
            track->SetLineWidth(iWght);
        }
        else {
            track->SetLineColor(kGreen);
            track->SetLineWidth(1);
        }

        track->SetLineStyle(3);

        TVector3 ep = traj->GetInitialPosition().Vect()
            + (eKin*traj->GetInitialMomentum().Vect().Unit());
        track->SetPoint(0,
                        traj->GetInitialPosition().X(),
                        traj->GetInitialPosition().Y(),
                        traj->GetInitialPosition().Z());
        track->SetPoint(1,ep.X(),ep.Y(),ep.Z());

        fTrajectoryList->AddElement(track);
    }
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
