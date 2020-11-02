#ifndef CubeBuildPairwiseVertices_hxx_seen
#define CubeBuildPairwiseVertices_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconVertex.hxx>

#include <TVector3.h>

namespace Cube {
    class BuildPairwiseVertices;
};

/// Build vertices and connect tracks to the vertices.  This uses the tracks
/// in the input algorithm final container.  The tracks are not copied to the
/// final container.  A vertex is defined as a position where two or more
/// tracks are expected to have connected.  Tracks can be separated from the
/// vertex (e.g. the projected "tracks" from a pizero), but a vertex cannot
/// occur in the middle of a track (i.e. only the end points of a track are
/// considered).  The allowed overlap between a tracks can be set.  This is
/// implemented with a "Kalman-like" incremental fit for the vertex position.
class Cube::BuildPairwiseVertices: public Cube::Algorithm {
public:

    BuildPairwiseVertices(const char* name="BuildPairwiseVertices");
    virtual ~BuildPairwiseVertices();

    /// Apply the algorithm.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// Set the acceptable likelihood (between 0.0 and 1.0) to accept a
    /// vertex. The likelihood is calculated based on the Chi2/DOF of a
    /// vertex, and a typical value might be 0.05.
    void SetLikelihoodCut(double v) {fLikelihoodCut = v;}

    /// Set the allowed overlap for a track.  The vertex must overlap the
    /// track by less than this distance for one end of a track (A typical
    /// value is 10 mm which means the track can "overshoot" the vertex by 1
    /// cube).
    void SetOverlap(double v) {fMaxOverlap = v;}

    /// Set the allowed distance from the track to the vertex.  The vertex
    /// must be within this distance of one end of a track (A typical value is
    /// 100 mm).
    void SetDistance(double v) {fMaxDistance = v;}

    /// Set the maximum allowed approach distance to be considered for a
    /// vertex.
    void SetMaxApproach(double v) {fMaxApproach = v;}

    /// Set the minimum track length to be considered for a vertex.  Tracks
    /// shorter than this should be added to identified vertices (by a
    /// different algorithm), but not used to select a vertex.
    void SetMinTrackLength(double v) {fMinTrackLength = v;}

private:

    // Stop combining vertices when the likelihood goes below this.  This is a
    // number between 0.0 and 1.0.  The likelihood is calculated with
    // TMath::Prob.
    double fLikelihoodCut;

    // The maximum allowed overshoot for a track with a vertex.  This is only
    // applied when building the initial vertices.  It is not enforced when
    // the vertices are combined.  A typical value would be 1 cm so that the
    // track can overlap the vertex position by one cube.
    double fMaxOverlap;

    // The maximum allowed distance between the track and the vertex.  This
    // prevents creating vertices from tracks that are a long way from the
    // intersection point.
    double fMaxDistance;

    // The maximum allowed approach distance.  Tracks must pass within this
    // distance of each other to be formed into a vertex.  This is not
    // enforced when vertices are combined.
    double fMaxApproach;

    // The minimum track length for adding to a vertex.
    double fMinTrackLength;

    // Find a rough estimate of the track length.  This is the distance
    // between the front and back positions.
    double TrackLength(Cube::Handle<Cube::ReconTrack> track);

    // For two lines (A,Ad) and (B,Bd), find the distance of closest approach.
    // The points, A and B, with the directions, Ad and Bd, define the lines.
    double ClosestApproach(const TVector3& A, const TVector3& Ad,
                           const TVector3& B, const TVector3& Bd);

    // Find the closest approach between two tracks.
    double ClosestApproach(Cube::Handle<Cube::ReconTrack> track1,
                           Cube::Handle<Cube::ReconTrack> track2);

    // Find the travel distance from A along Ad to the point of closest
    // approach.
    double TravelDistance(const TVector3& A, const TVector3& Ad,
                          const TVector3& B, const TVector3& Bd);

    // Find the best vertex based on the closest approach of the line (A,Ad)
    // and (B,Bd).  This needs the position variance at for each line at the
    // point of closest approach.
    TLorentzVector PairVertex(
        const TLorentzVector& A, const TVector3& Ad, double Avar,
        const TLorentzVector& B, const TVector3& Bd, double Bvar);

    // Build a vertex based on two tracks.  This looks at both the front and
    // back ends of each track.
    Cube::Handle<Cube::ReconVertex>
    MakePairVertex(Cube::Handle<Cube::ReconTrack> track1,
                   Cube::Handle<Cube::ReconTrack> track2);

    // Combine two vertices into a third.
    Cube::Handle<Cube::ReconVertex>
    CombineVertices(Cube::Handle<Cube::ReconVertex> vertex1,
                    Cube::Handle<Cube::ReconVertex> vertex2);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
