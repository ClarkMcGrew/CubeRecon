#include <iostream>
#include <iomanip>

#include <TROOT.h>

#include "CubeReconObject.hxx"
#include "CubeReconNode.hxx"
#include "CubeReconState.hxx"

ClassImp(Cube::ReconObject);

namespace {
    // This keeps track of the previous ReconObject unique identifier.
    static UInt_t gReconObjectId = 0;
}

Cube::ReconObject::ReconObject()
    : TNamed("unnamed","Reconstruction Object"),
      fQuality(0), fState(NULL), fNodes(NULL), fStatus(0), fNDOF(0) {
}

Cube::ReconObject::ReconObject(const char* name, const char* title)
    : TNamed(name,title),
      fQuality(0), fState(NULL), fNodes(NULL), fStatus(0), fNDOF(0) {
    SetUniqueID(++gReconObjectId);
}

Cube::ReconObject::ReconObject(const Cube::ReconObject& object)
    : TNamed(object), fState(NULL), fNodes(NULL) {
    SetUniqueID(++gReconObjectId);

    fQuality = object.GetQuality();
    fStatus = object.GetStatus();
    fNDOF = object.GetNDOF();

    if (object.GetHitSelection()) {
        Cube::Handle<Cube::HitSelection> hits(new HitSelection());
        std::copy(object.GetHitSelection()->begin(),
                  object.GetHitSelection()->end(),
                  std::back_inserter(*hits));
        SetHitSelection(hits);
    }

    if (object.GetConstituents()) {
        AddConstituents(object.GetConstituents());
    }

    // Don't copy the nodes and the state since they depend of the object type
}

Cube::ReconObject::~ReconObject() {
    if (fNodes) delete fNodes;
    if (fState) delete fState;
}

Cube::ReconObject::Status Cube::ReconObject::GetStatus() const {
    // Notice that this returns both the status bits and the deetector bits.
    return fStatus;
}

void Cube::ReconObject::SetStatus(Cube::ReconObject::Status status) {
    fStatus |= (status & kStatusMask);
}

void Cube::ReconObject::ClearStatus(Cube::ReconObject::Status status) {
    fStatus &= ~ (status & kStatusMask);
}

bool Cube::ReconObject::CheckStatus(Cube::ReconObject::Status status) const {
    return 0 != (fStatus & (status & kStatusMask));
}

Cube::ReconObject::Status Cube::ReconObject::GetDetectors() const {
    return (fStatus & kDetectorMask);
}

void Cube::ReconObject::AddDetector(Cube::ReconObject::Status status) {
    fStatus |= (status & kDetectorMask);
}

bool Cube::ReconObject::UsesDetector(Cube::ReconObject::Status status) const{
    return 0 != (fStatus & (status & kDetectorMask));
}

void Cube::ReconObject::RemoveDetector(Cube::ReconObject::Status status) {
    fStatus &= ~(status & kDetectorMask);
}

void Cube::ReconObject::AddConstituent(Cube::Handle<Cube::ReconObject> obj) {
    if (!fConstituents) {
        fConstituents
            = Cube::Handle<Cube::ReconObjectContainer>(
                new ReconObjectContainer("constituents",
                                         "Constituents of this object"));
    }
    fConstituents->push_back(obj);
}

void Cube::ReconObject::AddConstituents(
    Cube::Handle<Cube::ReconObjectContainer> objs) {
    Cube::ReconObjectContainer::const_iterator it;
    for (it = objs->begin(); it!= objs->end(); ++it) {
        AddConstituent((*it));
    }
}

std::string Cube::ReconObject::ConvertStatus() const {
    std::string s("(");
    bool notFirst = false;
    if (CheckStatus(kRan)) {
        s+= "ran";
        notFirst = true;
    }
    if (CheckStatus(kSuccess)) {
        if (notFirst) s+=":";
        s+= "success";
        notFirst = true;
    }
    if (CheckStatus(kChi2Fit)) {
        if (notFirst) s+=":";
        s+= "chi2";
        notFirst = true;
    }
    if (CheckStatus(kLikelihoodFit)) {
        if (notFirst) s+=":";
        s+= "likelihood";
        notFirst = true;
    }
    if (CheckStatus(kKalmanFit)) {
        if (notFirst) s+=":";
        s+= "kalman";
        notFirst = true;
    }
    if (CheckStatus(kStocasticFit)) {
        if (notFirst) s+=":";
        s+= "stocastic";
        notFirst = true;
    }
    s+= ")";
    return s;
}

std::string Cube::ReconObject::ConvertDetector() const {
    std::string s("(");
    bool notFirst = false;
    if (UsesDetector(kDST)) {
        if (notFirst) s+="-";
        s+= "DST";
        notFirst = true;
    }
    s+= ")";
    return s;
}

void Cube::ReconObject::ls(Option_t *opt) const {
    ls_base(opt);
    // The class specific stuff.

    TROOT::IncreaseDirLevel();
    std::string option(opt);
    if (fState) {
        TROOT::IncreaseDirLevel();
        fState->ls(opt);
        TROOT::DecreaseDirLevel();
    }
    if (fNodes && (option.find("dump") != std::string::npos
                   || option.find("recon") != std::string::npos)) {
        TROOT::IncreaseDirLevel();
        fNodes->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    TROOT::DecreaseDirLevel();
}

void Cube::ReconObject::ls_base(Option_t *opt) const {
    // Print the standard header
    std::string option(opt);
    TROOT::IndentLevel();
    std::cout << ClassName() << "(" << this << ")::"
              << GetName();
    if (option.find("title")!= std::string::npos) {
        std::cout << "/" << GetTitle();
    }
    if (GetUniqueID() > 0) std::cout << " uid: " << GetUniqueID();
    if (option.find("size") != std::string::npos && Class()) {
        std::cout << " size: " << Class()->Size() << " b";
    }
    std::cout << std::endl;
    // Do the class specific stuff.
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "Status: " << ConvertStatus();
    std::string tmp;
    tmp = GetAlgorithmName();
    if (tmp == "") tmp = "not-set";
    if (tmp != "") std::cout << " Algorithm: \"" << tmp << "\"";
    tmp = ConvertDetector();
    if (tmp == "()") tmp = "(not-set)";
    std::cout << "  In: " << tmp;
    std::cout << std::endl;

    TROOT::IndentLevel();
    std::ios::fmtflags save = std::cout.flags();
    std::cout << "Quality: " << std::setprecision(4) << GetQuality()
              << " for " << GetNDOF() << " d.o.f." << std::endl;
    std::cout.flags(save);
    TROOT::DecreaseDirLevel();
}

ClassImp(Cube::ReconObjectContainer);

Cube::ReconObjectContainer::ReconObjectContainer()
    : TNamed("unnamed","Recon Object Container") {}

Cube::ReconObjectContainer::ReconObjectContainer(const char* name,
                                                 const char* title)
    : TNamed(name,title) {}

Cube::ReconObjectContainer::~ReconObjectContainer() {}

void Cube::ReconObjectContainer::push_back(
    Cube::Handle<Cube::ReconObject> data) {
    std::string name = data->GetName();
    if (name == "unnamed") data->SetName(data->ClassName());
    std::vector< Cube::Handle<Cube::ReconObject> >::push_back(data);
}

void Cube::ReconObjectContainer::ls(Option_t* opt) const {
    // Print the standard header
    std::string option(opt);
    TROOT::IndentLevel();
    std::cout << ClassName() << "(" << this << ")::"
              << GetName();
    if (option.find("title")!= std::string::npos) {
        std::cout << "/" << GetTitle();
    }
    if (GetUniqueID() > 0) std::cout << " uid: " << GetUniqueID();
    if (option.find("size") != std::string::npos && Class()) {
        std::cout << " size: " << Class()->Size() << " b";
    }
    std::cout << std::endl;
    // Do the class specific stuff.
    TROOT::IncreaseDirLevel();
    for (const_iterator t = begin();
         t != end();
         ++t) {
        t->ls(opt);
    }
    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
