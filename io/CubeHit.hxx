#ifndef Hit_hxx_seen
#define Hit_hxx_seen

#include "CubeHandle.hxx"

#include <TObject.h>
#include <TVector3.h>

#include <vector>
#include <string>
#include <map>

namespace Cube {
    class Hit;
    class WritableHit;
}

/// A calibrated hit detector element where the hit information has been
/// reconstructed from the raw hits.  A Hit object represents a hit that
/// has been generated from one or more hits where the position and other
/// properties are calculated from the constituient hits.  A common usage is to
/// represent a hit that has been constructed from two or more hits on a wires
/// that cross at a single X/Y positon.  This is a lighter object that
/// TReconCluster and is primarily used for building up 2D hits into 3D
/// objects.  Complex combinations of hits should be placed into a
/// TReconCluster.
///
/// The Hit class can't be directly instantiated.  It is created
/// using the WritableHit class and is accessed as a Hit class.
class Cube::Hit : public TObject {
public:
    /// Define the status bits used by the Hit object.  These can't collide
    /// with any status bits defined in TObject (the parent class for Hit),
    /// and none of the Hit children can define a status bit that collides
    /// with these definitions.  Bits 14 to 23 are available for use.
    enum EStatusBits {
        kInvalidTime        = BIT(18),   // if time is not valid.
        kInvalidCharge      = BIT(19),   // if charge is not valid.
    };

    Hit();
    Hit(const WritableHit& val);
    virtual ~Hit();

    /// Return a unique identifier for the sensor that generated this hit.
    /// This should be sufficient for the reconstruction to find the sensor(s)
    /// in the geometry.  The reconstruction will use this to determine if two
    /// hits occurred on the same sensor.
    virtual int GetIdentifier(void) const;

    /// Return the calibrated "charge" for the hit.
    virtual double GetCharge(void) const;

    /// Return the uncertainty in the charge for the hit.
    virtual double GetChargeUncertainty() const;

    /// Return true if the calibrated charge is valid.
    virtual bool HasValidCharge(void) const {return !TestBit(kInvalidCharge);}

    /// Return the calibrated "time" for the hit representing the mean (or
    /// central) time of the hit.
    virtual double GetTime(void) const;

    /// Return the uncertainty on the hit time.
    virtual double GetTimeUncertainty(void) const;

    /// Return true if the calibrated time is valid.
    virtual bool HasValidTime(void) const {return !TestBit(kInvalidCharge);}

    /// The position of this hit.
    virtual const TVector3& GetPosition(void) const;

    /// Return the uncertainty of the hit position (approximated as diagonal).
    virtual const TVector3& GetUncertainty(void) const;

    /// Return the physical size of the hit in the detector.  This is used to
    /// find out if two hits are in contact.  When the hit doesn't have a well
    /// defined size, this may be the RMS.
    virtual const TVector3& GetSize(void) const;

    /// Return a constituent hit.  If the index is out of range, this will
    /// throw an EHitOutOfRange exception.  By default this will throw an
    /// EHitOutOfRange, but it may be over-ridden in a derived class.
    virtual Cube::Handle<Cube::Hit> GetConstituent(int i=0) const;

    /// Return the number of "sub" hits that part of this hit.
    virtual int GetConstituentCount() const;

    /// Return an integer to identify where the truth information for this
    /// hit came from.  The exact definition will be nailed down as I'm
    /// working, but it is probably going to be the index of the hit segment.
    virtual int GetContributor(int i=0) const;

    /// Return the number of contributor entries.
    virtual int GetContributorCount() const;

    /// Check if a property exists.
    virtual bool HasProperty(std::string name) const;

    /// Get a property value.  Common propertie values that may exist are the
    /// attenuation constants ("Atten1", "Atten2", "Ratio12", "Reflectivity"),
    /// and other details that might be in a database.  This will throw a
    /// runtime_error if the property doesn't exist.
    virtual double GetProperty(std::string name) const;

    /// Print the hit information.
    virtual void ls(Option_t *opt = "") const;

protected:

    /// Set the validity of the calibrated time.
    virtual void SetTimeValidity(bool valid) {SetBit(kInvalidTime, !valid);}

    /// Set the validity of the calibrated charge.
    virtual void SetChargeValidity(bool valid) {SetBit(kInvalidCharge, !valid);}

private:

    /// Fill all of the geometry related fields from the single hits and
    /// apply corrections.
    void Initialize();

protected:

    /// A hit identifier.  This is used to describe the type of the hit.
    Int_t fIdentifier;

    /// The measured "charge" for this hit.
    Float_t fCharge;

    Float_t fChargeUncertainty;

    /// The measured "time" for this hit.
    Float_t fTime;

    /// The reconstructed time uncertainty
    Float_t fTimeUncertainty;

    /// The reconstructed position of the hit in global coordinates.
    TVector3 fPosition;

    /// The uncertainty of the hit position in global coordinates.  The
    /// uncertainty of the hit is in the global coordinates, and is by
    /// approximation diagonal (to save space).   For a more complete
    /// representation of the covariance, use a ReconCluster.
    TVector3 fUncertainty;

    /// The physical size of the hit in the detector.  This is used to find
    /// out if two hits are in contact.  When the hit doesn't have a well
    /// defined size, this may be the RMS.
    TVector3 fSize;

    /// Any ts that make up this reconstructed hit.
    std::vector< Cube::Handle < Cube::Hit > > fConstituents;

    /// Integer identifiers for truth information about where this hit came
    /// from.  This is usually the index of the hit segments,
    std::vector< int > fContributors;

    /// A map of properties associated with this it.
    std::map<std::string, double> fProperties;

    ClassDef(Hit,1);
};

/// Provide a writable interface to a Hit that can be used to fill the
/// object.  All of the Hit objects added to a single Hit object must
/// be in the same geometry object (e.g. all in the save bar).  The geometry
/// information (i.e. the geometry id, spread, &c) is taken from the first
/// hit.  The other hits are not used by this class, but are checked to make
/// sure that they are in the same geometry volume.  If you need to combine
/// hits from different geometry volumes, use a TReconCluster.
class Cube::WritableHit : public Hit {
public:
    WritableHit();
    WritableHit(const Hit& h);
    WritableHit(const WritableHit& h);
    virtual ~WritableHit();

    /// Set the identifier.
    void SetIdentifier(Int_t identifier);

    /// Set the charge for the hit.
    void SetCharge(double q);

    void SetChargeUncertainty(double unc);

    /// Set the time for the hit.
    void SetTime(double t);

    /// Set the time uncertainty for the hit.
    void SetTimeUncertainty(double tunc);

    /// Set the position for the hit.
    void SetPosition(const TVector3& pos);

    /// Set the position uncertainty for the hit in global coordinates.  The
    /// position covariance is approximated as diagonal.  For a more complete
    /// representation of the covariance, use a ReconCluster.
    void SetUncertainty(const TVector3& unc);

    /// Set the physical size of the hit.
    void SetSize(const TVector3& siz);

    /// Add one more hit to an existing WritableHit.
    void AddHit(Cube::Handle<Cube::Hit> hit);

    /// Add a contributor to the writable hit.  This is an integer identifier
    /// for where the truth information came from.  It's usually the hit
    /// segment identifier.
    void AddContributor(int contrib);

    /// Set a property value
    void SetProperty(std::string name, double value);

    ClassDef(WritableHit,1);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
