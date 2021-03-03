#ifndef CubeReconBase_hxx_seen
#define CubeReconBase_hxx_seen

#include <TObject.h>

#include "CubeLog.hxx"
#include "CubeHandle.hxx"
#include "CubeHitSelection.hxx"
#include "CubeReconState.hxx"

namespace Cube {
    class ReconObject;
    class ReconObjectContainer;
    class ReconNodeContainer;
}

/// Define a base class for all of the reconstruction objects.  This contains
/// the shared methods required to implement the basic reconstruction object
/// functionality.  Each ReconObject object has it's TObject::fUniqueID value
/// set.  These values can be used to associate global reconstruction objects
/// to sub-detector reconstruction objects.
class Cube::ReconObject : public TNamed {
public:
    /// The bits defining the state of the reconstruction that created this
    /// object.  The state bits come in two flavors: Status bits are used to
    /// summarize the completion status of the algorithm.  Detector bits are
    /// used to record which detectors provided information to the algorithm.
    /// The status bits are accessed using SetStatus() and GetStatus().  The
    /// detector bits are accessed using AddDetector() and UsesDetector().
    /// The full status field can be accessed using the GetStatus() method.
    typedef enum {
        /// The algorithm ran to completion.  This bit is set using
        /// SetStatus() and checked using CheckStatus()
        kRan            = BIT(0),

        /// The algorithm succeeded and the result is reliable.  This bit is
        /// set using SetStatus and checked using CheckStatus().
        kSuccess        = BIT(1),

        /// A chi2 fitting method was used and the quality corresponds to a
        /// chi2.  This bit is set using SetStatus and checked using
        /// CheckStatus().
        kChi2Fit        = BIT(2),

        /// A maximum likelihood fitting methods was used.  This bit is set
        /// using SetStatus and checked using CheckStatus().
        kLikelihoodFit  = BIT(3),

        /// A analytic Bayesian filter fitting algorithm was used (i.e. a form
        /// of a Kalman filter).  This bit is set using SetStatus and checked
        /// using CheckStatus().
        kKalmanFit      = BIT(4),

        /// A stocastic Bayesian filter fitting algorithm was used (i.e. a form
        /// of a particle filter).  This bit is set using SetStatus and checked
        /// using CheckStatus().
        kStocasticFit      = BIT(5),

        /// A mask of all the fitter bits.
        kFitMask       = (kChi2Fit|kLikelihoodFit|kKalmanFit|kStocasticFit),

        /// A mask of all of the status bits.
        kStatusMask  = (kRan|kSuccess|kFitMask),

        /////////////////////////////////////////////////////////////////
        // The remaining bits define which detectors contribute information
        // to the reconstruction object.
        /////////////////////////////////////////////////////////////////

        /// A mask to check if the result includes data from the TPC.  This
        /// mask is checked using UsesDetector().
        kDST         = BIT(12),

        kTPCBack     = BIT(13),
        kTPCTop      = BIT(14),
        kTPCBot      = BIT(15),
        kTPC         = (kTPCBack|kTPCBot|kTPCTop),

        kECalBarrel  = BIT(16),
        kECalLeft    = BIT(17),
        kECalRight   = BIT(18),
        kECal        = (kECalBarrel|kECalLeft|kECalRight),

        /// A mask for all of the detector bits.
        kDetectorMask   = (kDST|kTPC|kECal),

    } StateBits;

    /// A status value with one or more bits sets
    typedef unsigned long Status;

public:
    virtual ~ReconObject();

    /// Get the name of the algorithm that created this reconstruction object.
    std::string GetAlgorithmName() const {return fAlgorithm;}

    /// Set the name of the algorithm.
    void SetAlgorithmName(const char* name) {fAlgorithm = name;}

    /// Check the status of the reconstruction that created this
    /// reconstruction result.  This method checks if any of the bit
    bool CheckStatus(Cube::ReconObject::Status status) const;

    /// Set the status of the algorithm that generated this object.
    void SetStatus(Cube::ReconObject::Status status);

    /// Clear the status of the algorithm that generated this object.
    void ClearStatus(Cube::ReconObject::Status status);

    /// Get the status field for this object.
    Cube::ReconObject::Status GetStatus() const;

    /// Get the detector bits that are set for this object.
    Cube::ReconObject::Status GetDetectors() const;

    /// Return true if the object uses a particular detector
    bool UsesDetector(Cube::ReconObject::Status detector) const;

    /// Add one or more detector to the object.
    void AddDetector(Cube::ReconObject::Status detector);

    /// Remove a detector from the object.
    void RemoveDetector(Cube::ReconObject::Status detector);

    /// Return the goodness of fit for the reconstruction.
    double GetQuality() const {return fQuality;}

    /// Set the goodness of the fit from the reconstruction.
    void SetQuality(double quality) {fQuality = quality;}

    /// Get the number of degrees of freedom in the reconstruction.
    double GetNDOF() const {return fNDOF;}

    /// Set the number of degrees of freedom.
    void SetNDOF(double n) {fNDOF = n;}

    /// Get the state associated with this reconstruction object.  The state
    /// is created by the derived class.  This should generally be used in
    /// user code with a specific Cube::Handle type.  For instance, if you
    /// want to access a Cube::ReconCluster state, you should write code
    /// like:
    ///
    /// \code
    /// Cube::Handle<Cube::TClusterState> state = recObj->GetReconState();
    /// if (state) {
    ///    // the state
    /// }
    /// \endcode
    ///
    /// Which hides the dynamic cast and makes sure that the local variable
    /// has the right type.
    ///
    /// The state holds (most of) the reconstruction object parameter values.
    /// And is used to get and set the values and covariances.  As an example,
    /// to set the direction of the ReconTrack object:
    /// \code
    /// Cube::Handle<Cube::TTrackState> trackState
    ///               = reconTrack->GetReconState();
    /// TVector3 theDir(0,0,1);
    /// trackState->SetDirection(theDir); // This works
    /// trackState->SetDirection(0,0,1);  // This also works.
    /// // Assume dirCov is a TMatrixD holding the covariance.
    /// for (int i=0; i< TMDirectionState::GetSize(); ++i) {
    ///     for (int j=0; j< TMDirectionState::GetSize(); ++j) {
    ///         trackState->SetCovariance(i+trackState->GetDirectionIndex(),
    ///                                   j+trackState->GetDirectionIndex(),
    ///                                   dirCov(i,j));
    ///     }
    /// }
    /// \endcode
    Cube::Handle<Cube::ReconState> GetReconState() const {
        if (!fState) {
            CUBE_ERROR << "ReconObject with NULL State: "
                       << "State must be created in the derived class."
                       << std::endl;
            throw std::runtime_error("ReconObject Error");
        }
        return Cube::Handle<Cube::ReconState>(fState,false);
    }

    /// Get the parameter values at each stage of the reconstruction.
    /// - For the track and particle objects (see the \ref reconTrack and \ref
    ///   reconParticle sections), the parameter values will be the associated
    ///   with each intermediate position.
    /// - For an incremental fit, the parameter values are associated with each
    ///   stage of the fit.
    /// - For other objects, the parameter values will be associated with the
    ///   reconstruction objects that have been added to the fit.
    const Cube::ReconNodeContainer& GetNodes() const {
        if (!fNodes) {
            CUBE_ERROR << "ReconObject without a ReconNodeContainer:"
                       << "Must be created in derived class constructor"
                       << std::endl;
            throw std::runtime_error("ReconObject Error");
        }
        return *fNodes;
    }

    /// Provide a non-constant version of the node container so that nodes can
    /// be added.
    Cube::ReconNodeContainer& GetNodes() {
        if (!fNodes) {
            CUBE_ERROR << "ReconObject without a ReconNodeContainer:"
                       << "Must be created in derived class constructor"
                       << std::endl;
            throw std::runtime_error("ReconObject Error");
        }
        return *fNodes;
    }

    /// Provide the hits associated with this recon object.  The hits are
    /// saved as an entry in the data vector (with name "hits"), and may not
    /// be present for all objects.
    Cube::Handle<Cube::HitSelection> GetHitSelection() const {return fHits;}

    /// Add a hit selection to the recon object.  The recon object takes
    /// ownership of the hit (as per the pointer ownership policy), so you
    /// must not add a hit selection currently owned by another object.  The
    /// name of the hit selection will be changed to "hits".
    void SetHitSelection(Cube::Handle<Cube::HitSelection> hits) {fHits = hits;}

    /// Return the median X, Y, Z, and T.  Be careful since there may
    /// not be an actual hit at the median position (X, Y, and Z can be
    /// from different hits).
    TLorentzVector GetMedian() const;

    /// Provide the reconstruction objects used in the algorithm that created
    /// this object.  The consituients are used to help understand how the
    /// object was fit, and are different from the nodes which represent
    /// actual intermediate steps of the fitting process.
    Cube::Handle<Cube::ReconObjectContainer> GetConstituents() const
        {return fConstituents;}

    /// Add a new constituent to the ReconObject.  This stores the
    /// constituents in the order in which they are added.  The constituents
    /// are saved in a ReconObjectContainer, and are available using the
    /// GetConstituents method.
    void AddConstituent(Cube::Handle<Cube::ReconObject> obj);

    /// Add several constituents in one go by adding all objects in a
    /// ReconObjectContainer
    void AddConstituents(Cube::Handle<Cube::ReconObjectContainer> objs);

    /// @{Override methods in the base TObject class.
    virtual void ls(Option_t* opt = "") const;
    virtual void Print(Option_t* opt = "") const {ls(opt);}
    /// @}

    /// Turn the object status into a string.
    std::string ConvertStatus() const;

    /// Turn the contributing detectors into a string.
    std::string ConvertDetector() const;

protected:
    /// Default constructor.
    ReconObject();

    /// Construct a named reconstruction object.  This constructs a new empty
    /// ReconObject object and is defined as explicit so that it will not
    /// function as a type conversion between a const char* and the ReconObject
    /// class.  This can only be used in one of the derived classes.
    explicit ReconObject(const char* name,
                        const char* title = "Reconstruction Object");

    /// Copy the object (this is not a deep copy).
    ReconObject(const Cube::ReconObject& object);

    /// Implement the listing of the base class fields.  This is used by ls()
    void ls_base(Option_t* opt) const;

    /// The all the hits used in this object.
    Cube::Handle<Cube::HitSelection> fHits;

    /// Any sub objects that contribute to this one.
    Cube::Handle<Cube::ReconObjectContainer> fConstituents;

    /// The quality of reconstruction;
    float fQuality;

    /// The state for this object.
    Cube::ReconState* fState;

    /// The nodes for the object.
    Cube::ReconNodeContainer* fNodes;

    /// The status of the reconstruction
    Cube::ReconObject::Status fStatus;

    /// The number of degrees of freedom.
    float fNDOF;

    /// The name of the reconstruction algorithm.
    std::string fAlgorithm;

    ClassDef(ReconObject,1);
};

/// A container class for reconstruction objects inheriting from a std::vector
/// of Handle<ReconObject>.
class Cube::ReconObjectContainer
    : public TNamed,
      public std::vector< Cube::Handle<Cube::ReconObject> > {
public:
    ReconObjectContainer();
    ReconObjectContainer(const char* name,
                          const char* title = "Recon Object Container");
    virtual ~ReconObjectContainer();

    virtual void push_back(Cube::Handle<Cube::ReconObject> data);

    /// @{Override methods in the base TObject class.
    virtual void ls(Option_t* opt = "") const;
    /// @}

private:
    ClassDef(ReconObjectContainer,1);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
