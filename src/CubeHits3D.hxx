#ifndef CubeHits3D_hxx_seen
#define CubeHits3D_hxx_seen
#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>

namespace Cube {
    class Hits3D;
};

/// Take a TAlgorithmResult containing a hit selection of 2D (MPPC) hits and
/// build the 2D hits into 3D hits (class TReconHit) that are returned as new
/// hit selection.  The intention is that this will be a lower level algorithm
/// which is used to build 3D hits that are then transformed into clusters.
/// For each input selection 2D hits, all possible combinations of 3D hits are
/// formed which are saved in a single output TReconCluster object.
class Cube::Hits3D
    : public Cube::Algorithm {
public:
    Hits3D();
    virtual ~Hits3D();

    /// Take a TAlgorithmResult containing THitSelection objects with 2D fiber
    /// hits, and group them into a THitselection containing 3D XYZT hits.
    /// Since a THitSelection can be converted into a TAlgorithmResult this
    /// can also be called with a THandle to a THitSelection.  The output
    /// TAlgorithmResult will contain:
    ///
    ///   * unused -- A THitSelection of any 2D hits that were not used in 3D
    ///               hits.  It should be empty.
    ///
    ///   * used -- A THitSelection with all of the 3D TReconHit objects from
    ///              this algorithm.  This is the last hit selection added, so
    ///              it's the returned as the default result of
    ///              "GetHitSelection()".
    ///
    ///   * final -- A TReconObjectContainer with a cluster constructed from
    ///               all of the 3D hits.
    Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// Set the type of charge sharing to apply:
    ///
    /// 1: OptimizeCubes: Apply (almost) the same formalism as for the
    ///       MaximumEntropy version.  Share the charge by predicting the
    ///       measurement in each fiber based on the deposit in each cube.
    ///       The deposit in a cube is constrained to be positive.  The
    ///       maximum entropy prior is that all cubes should have the same
    ///       charge (it's a very weak prior).  There is a prior applied based
    ///       on an approximate entropy calculation so that this runs more
    ///       quickly.
    ///
    /// 2: ApplyConstraints: Share the charge assuming that all of the fibers
    ///       in a cube will measure the same amount of deposit, constrained
    ///       by the charge present at the fiber being equal to the measured
    ///       charge.  The constraint on the charge at the fiber is applied
    ///       using a Lagrange multiplier.  The constraint at the cube is the
    ///       likelihood that all three measurements came from the same mean.
    ///
    ///
    /// 3: MaximumEntropy: Share the charge applying a Bayesian probability
    ///       with a maximum entropy prior.  The probability is based on
    ///       predicting the measurement in each fiber based on the deposit in
    ///       each cube.  The deposit in a cube is constrained to be positive.
    ///       The maximum entropy prior is that all cubes should have the same
    ///       charge (it's a very weak prior).  This is impossibly slow for
    ///       large events and has precision problems.  The ApplyConstraints
    ///       and OptimizeCubes versions are better.
    ///
    /// Any other value won't apply charge sharing.
    void SetShareCharge(int i) {fShareCharge = i;}

    typedef std::vector<std::pair<double, Cube::Handle<Cube::Hit>>> FiberTQ;
private:

    /// Take three fiber hits (one can be empty) and make a 3D
    /// TWritableReconHit handle objects.  The hits need to be from different
    /// projections and the fibers MUST cross (may not be explicitly checked).
    /// This will return false if there is a problem constructing the hit.
    /// The new hit is pushed to the writableHits hit selection.
    ///
    /// NOTE: This has an odd convention for the first constituent of the
    /// TReconHit.  The first hit is a "fake" that is constructed to hold the
    /// geometry id of the cube.  The other constituents correspond to the
    /// fibers.  This is done to workaround a design problem in the oaEvent
    /// TReconHit.
    bool MakeHit(Cube::HitSelection& writableHits,
                 const Cube::Handle<Cube::Hit>& hit1,
                 const Cube::Handle<Cube::Hit>& hit2,
                 const Cube::Handle<Cube::Hit>& hit3) const;

    /// Take the hit times and charges and estimate the best time.  The return
    /// value is the time (first), and rms (second).
    std::pair<double,double> HitTime(FiberTQ& fiberTQ) const;

    /// The type of charge sharing to apply. 0) no sharing, 1) OptimizeCubes,
    /// 2) ShareContraints, or 3) MaximumEntropy.  The default is 1.
    int fShareCharge;

    /// Flag whether the final shared charge should be adjusted to match the
    /// total measured charge in the MPPCs.  This tends to increase the total
    /// charge since during the sharing, some MPPCs will not be matched into a
    /// hit.
    int fConserveChargeSum;

    /// The velocity of the light in the fiber.
    double fLightSpeed;

};
#endif
