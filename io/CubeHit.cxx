#include "CubeHit.hxx"
#include "TUnitsTable.hxx"

#include <TGeoManager.h>
#include <TROOT.h>

#include <exception>
#include <iostream>

ClassImp(Cube::Hit);
ClassImp(Cube::WritableHit);

Cube::Hit::Hit()
    : fIdentifier(0), fCharge(-9999), fChargeUncertainty(-9999),
      fTime(-9999), fTimeUncertainty(-9999) {}

Cube::Hit::Hit(const Cube::WritableHit& h)
    : fIdentifier(h.fIdentifier), fCharge(h.fCharge),
      fChargeUncertainty(h.fChargeUncertainty),
      fTime(h.fTime), fTimeUncertainty(h.fTimeUncertainty),
      fPosition(h.fPosition), fUncertainty(h.fUncertainty), fSize(h.fSize),
      fConstituents(h.fConstituents), fContributors(h.fContributors),
      fProperties(h.fProperties) {}

Cube::Hit::~Hit() { }

int Cube::Hit::GetIdentifier(void) const {return fIdentifier;}

double Cube::Hit::GetCharge(void) const {return fCharge;}

double Cube::Hit::GetChargeUncertainty(void) const {
    return fChargeUncertainty;
}

double Cube::Hit::GetTime(void) const {return fTime;}

double Cube::Hit::GetTimeUncertainty(void) const {return fTimeUncertainty;}

const TVector3& Cube::Hit::GetPosition(void) const {
    return fPosition;
}

const TVector3& Cube::Hit::GetUncertainty(void) const {
    return fUncertainty;
}

const TVector3& Cube::Hit::GetSize(void) const {
    return fSize;
}

Cube::Handle<Cube::Hit> Cube::Hit::GetConstituent(int i) const {
    if (i<0 || fConstituents.size()<= (unsigned) i) {
        throw std::runtime_error("Constituent hit index out of range");
    }
    return fConstituents[i];
}

int Cube::Hit::GetConstituentCount() const {
    return fConstituents.size();
}

int Cube::Hit::GetContributor(int i) const {
    if (i<0 || fContributors.size()<= (unsigned) i) {
        throw std::runtime_error("Contributor hit index out of range");
    }
    return fContributors[i];
}

int Cube::Hit::GetContributorCount() const {
    return fContributors.size();
}

bool Cube::Hit::HasProperty(std::string name) const {
    std::map<std::string,double>::const_iterator p = fProperties.find(name);
    return (p != fProperties.end());
}

double Cube::Hit::GetProperty(std::string name) const {
    std::map<std::string,double>::const_iterator p = fProperties.find(name);
    if (p != fProperties.end()) return p->second;
    throw std::runtime_error("Cube::Hit Property does not exist");
    return 0.0;
}

// WritableHits.
Cube::WritableHit::WritableHit() {
    // Explicitly initialize here for clarity.
    fIdentifier = 0;
    fCharge = -9999;
    fChargeUncertainty = -9999;
    fTime = -9999;
    fTimeUncertainty = -9999;
}

Cube::WritableHit::WritableHit(const Cube::Hit& h)
    : Cube::Hit(h) {}

Cube::WritableHit::WritableHit(const Cube::WritableHit& h)
    : Cube::Hit(h) {}

Cube::WritableHit::~WritableHit() {}

//////////////////////////////////////////////////
// Setter methods for Cube::WritableHit
//////////////////////////////////////////////////
void Cube::WritableHit::AddHit(Cube::Handle<Cube::Hit> hit) {
    fConstituents.push_back(hit);
}

void Cube::WritableHit::AddContributor(int c) {fContributors.push_back(c);}

void Cube::WritableHit::SetIdentifier(int id) {fIdentifier = id;}

void Cube::WritableHit::SetCharge(double q) {fCharge = q;}

void Cube::WritableHit::SetChargeUncertainty(double q) {
    fChargeUncertainty = q;
}

void Cube::WritableHit::SetTime(double t) {fTime = t;}

void Cube::WritableHit::SetTimeUncertainty(double tunc) {
    fTimeUncertainty = tunc;
}

void Cube::WritableHit::SetPosition(const TVector3& pos) {
    fPosition = pos;
}

void Cube::WritableHit::SetUncertainty(const TVector3& unc){
    fUncertainty = unc;
}

void Cube::WritableHit::SetSize(const TVector3& siz){
    fSize = siz;
}

void Cube::WritableHit::SetProperty(std::string name, double value) {
    fProperties[name] = value;
}

void Cube::Hit::ls(Option_t *opt) const {
    int prec = std::cout.precision();
    std::cout << "Hit S:" << GetIdentifier()
              << " T: " << unit::AsString(GetTime(),
                                          GetTimeUncertainty(),"time")
              << " Q: " << unit::AsString(GetCharge(),
                                           GetChargeUncertainty(),
                                           "charge")
              << std::endl;
    TROOT::DecreaseDirLevel();

    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "P: " << unit::AsString(GetPosition(),"length")
              << std::endl;
    TROOT::DecreaseDirLevel();
    std::string option(opt);
    if (option.find("hits") != std::string::npos) {
        TROOT::IncreaseDirLevel();
        TROOT::IndentLevel();
        std::cout
            << "Size          "
            << std::setprecision(0)
            << unit::AsString(GetSize(),"length")
            << std::endl;
        TROOT::IndentLevel();
        std::cout << "Uncertainty: "
                  << std::setprecision(2)
                  << "("
                  << unit::AsString(GetUncertainty().X(),"length")
                  << ", "
                  << unit::AsString(GetUncertainty().Y(),"length")
                  << ", "
                  << unit::AsString(GetUncertainty().Z(),"length")
                  << ", "
                  << unit::AsString(GetTimeUncertainty(),"time")
                  << ")"
                  << std::endl;
        TROOT::DecreaseDirLevel();
    }
    std::cout.precision(prec);

    TROOT::IncreaseDirLevel();
    for (std::vector< Cube::Handle< Cube::Hit > >::const_iterator h
             = fConstituents.begin();
         h != fConstituents.end();
         ++h) {
        h->ls(opt);
    }
    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
