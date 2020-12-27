#include "ToolPrimaryId.hxx"

#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>

// Find the primary trajectory id (this is the trajectory that was started
// by a GEANT4 primary particle).
int Cube::Tool::PrimaryId(Cube::Event& event, int trajId) {
    Cube::Handle<Cube::G4Trajectory> traj = event.G4Trajectories[trajId];
    if (!traj) return -1;
    if (traj->GetParentId() < 0) return traj->GetTrackId();
    return PrimaryId(event,traj->GetParentId());
}

std::vector<int> Cube::Tool::AllPrimaries(Cube::Event& event) {
    std::vector<int> output;
    for (Cube::Event::G4TrajectoryContainer::iterator t
             = event.G4Trajectories.begin();
         t != event.G4Trajectories.end(); ++t) {
        if (!t->second) {
            // This happens when a non-existant trajectory id is accessed.
            // That creates a new entry in the map that doesn't have a
            // trajectory attached.
            std::cout << "Invalid trajectory entry for id " << t->first
                      << std::endl;
            continue;
        }
        // If the parent id is negative, this is a "real" primary.
        if (t->second->GetParentId() < 0) {
            output.push_back(t->second->GetTrackId());
            continue;
        }
        // Check if this is a photon from pizero decay
        //  NOT IMPLEMENTED
        // Check if this is a michel electron from muon decay
        //  NOT IMPLEMENTED
        // Make sure that the first particle created by a neutron is recorded
        //  NOT IMPLEMENTED... SHOULD THIS BE IMPLEMENTED?
    }
    return output;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
