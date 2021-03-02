#include "ToolPrimaryId.hxx"

#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeUnits.hxx>
#include <CubeG4Trajectory.hxx>

// Find the primary trajectory id (this is the trajectory that was started
// by a GEANT4 primary particle).
int Cube::Tool::PrimaryId(const Cube::Event& event, int trajId) {
    Cube::Event::G4TrajectoryContainer::const_iterator
        t = event.G4Trajectories.find(trajId);
    if (t == event.G4Trajectories.end()) return -1;
    Cube::Handle<Cube::G4Trajectory> traj = t->second;
    if (!traj) return -1;
    if (traj->GetParentId() < 0) return traj->GetTrackId();
    return PrimaryId(event,traj->GetParentId());
}

std::vector<int> Cube::Tool::AllPrimaries(const Cube::Event& event) {
    std::vector<int> output;
    for (Cube::Event::G4TrajectoryContainer::const_iterator t
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

std::vector<int> Cube::Tool::PrimaryVertex(const Cube::Event& event,
                                           int trajId) {
    std::vector<int> output;

    // Find the location of the primary associated with the trajectory id.
    int primId = Cube::Tool::PrimaryId(event,trajId);
    Cube::Event::G4TrajectoryContainer::const_iterator
        t = event.G4Trajectories.find(primId);
    if (t == event.G4Trajectories.end()) return output;
    if (!t->second) return output;
    TLorentzVector vertex = t->second->GetInitialPosition();

    // Add all of the primaries coming from vertex to an output vector.
    for (t = event.G4Trajectories.begin();
         t != event.G4Trajectories.end(); ++t) {
        if (!t->second) {
            // This happens when a non-existant trajectory id is accessed.
            // That creates a new entry in the map that doesn't have a
            // trajectory attached.
            continue;
        }

        //  Skip the trajectory if it is a secondary particle.
        if (0 <= t->second->GetParentId()) continue;

        // The actual vertex information is lost to /dev/null, so check if
        // this came from the same vertex by looking at the separation between
        // the initial positions.
        TLorentzVector sep = t->second->GetInitialPosition() - vertex;

        // Separated times so not the same interaction
        if (std::abs(sep.T()) > 1.0*unit::ns) continue;

        // Separated locations so not the same interaction.
        if (sep.Vect().Mag() > 1.0*unit::mm) continue;

        // We have a winner.
        output.push_back(t->second->GetTrackId());
    }

    return output;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
