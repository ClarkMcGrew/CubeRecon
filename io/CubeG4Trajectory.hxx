#ifndef CubeG4Trajectory_hxx_seen
#define CubeG4Trajectory_hxx_seen

#include <TLorentzVector.h>
#include <TObject.h>

namespace Cube {
    class G4Trajectory;
}

/// A class to save a G4 trajectory into a root output file without linking to
/// geant.  A trajectory is the truth information about the path of a particle
/// through the G4 simulation. It saves the parent trajectory that generated
/// this particle, the initial momentum of the particle.
class Cube::G4Trajectory : public TObject {
public:
    G4Trajectory(void)
        : fTrackId(-1), fParentId(-1),
          fPDGCode(0) {}

    virtual ~G4Trajectory();

    /// The TrackId of this trajectory.
    int GetTrackId() const {return fTrackId;}
    void SetTrackId(int i) {fTrackId = i;}

    /// The unique Id of the parent trajectory (The TrackId of the parent).
    int GetParentId() const {return fParentId;}
    void SetParentId(int i) {fParentId = i;}

    /// The PDG encoding of the particle.
    int GetPDGCode() const {return fPDGCode;}
    void SetPDGCode(int i) {fPDGCode = i;}

    /// The initial position of the particle
    const TLorentzVector& GetInitialPosition() const {return fInitialPosition;}
    void SetInitialPosition(const TLorentzVector& v) {fInitialPosition = v;}

    /// The initial momentum of the particle
    const TLorentzVector& GetInitialMomentum() const {return fInitialMomentum;}
    void SetInitialMomentum(const TLorentzVector& v) {fInitialMomentum = v;}

private:

    /// The TrackId of this trajectory.
    Int_t fTrackId;

    /// The unique Id of the parent trajectory (The TrackId of the parent).
    Int_t fParentId;

    /// The PDG encoding of the particle.
    Int_t fPDGCode;

    /// The initial position of the particle
    TLorentzVector fInitialPosition;

    /// The initial momentum of the particle
    TLorentzVector fInitialMomentum;

    ClassDef(G4Trajectory,1)
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
