#include "CubeAlgorithmResult.hxx"

#include "CubeAlgorithm.hxx"
#include "CubeLog.hxx"

#include <TROOT.h>

ClassImp(Cube::AlgorithmResult);

const Cube::AlgorithmResult Cube::AlgorithmResult::Empty;

Cube::AlgorithmResult::AlgorithmResult(): fStatusSummary("") { }

Cube::AlgorithmResult::AlgorithmResult(const char* name, const char* title)
    : TNamed(name, title), fStatusSummary("") { }

Cube::AlgorithmResult::AlgorithmResult(const Cube::Algorithm& algo)
    : TNamed(algo.GetName(), "An Algorithm Result"), fStatusSummary("") { }

Cube::AlgorithmResult::AlgorithmResult(const Cube::HitSelection& hits)
    : TNamed(hits.GetName(), "An Algorithm Result"), fStatusSummary("") {
    Cube::Handle<Cube::HitSelection> internalHits(
        new Cube::HitSelection(hits.GetName()));
    for (Cube::HitSelection::const_iterator h = hits.begin();
         h != hits.end();
         ++h) {
        internalHits->push_back(*h);
    }
    AddHitSelection(internalHits);
}

Cube::AlgorithmResult::~AlgorithmResult() { }

bool Cube::AlgorithmResult::IsEmpty() const {
    if (!fResultsContainer.empty()) return false;
    if (!fObjectContainers.empty()) return false;
    if (!fHitSelections.empty()) return false;
    return fStatusSummary.empty();
}

void Cube::AlgorithmResult::AddStatus(const char* s) {
    std::string status(s);
    fStatusSummary = "(" + status + ") " + fStatusSummary;
}

void Cube::AlgorithmResult::AddStatus(const std::string& status) {
    fStatusSummary = "(" + status + ") " + fStatusSummary;
}

void Cube::AlgorithmResult::SetStatus(const char* status) {
    fStatusSummary = status;
}

void Cube::AlgorithmResult::SetStatus(const std::string& status) {
    fStatusSummary = status;
}

std::string Cube::AlgorithmResult::GetStatus() const {
    return fStatusSummary;
}

void Cube::AlgorithmResult::AddAlgorithmResult(
    Cube::Handle<Cube::AlgorithmResult> result) {
    fResultsContainer.push_back(result);
}

Cube::Handle<Cube::AlgorithmResult>
Cube::AlgorithmResult::GetAlgorithmResult(const char* name) const {
    if (fResultsContainer.empty()) {
        return Cube::Handle<Cube::AlgorithmResult>();
    }
    if (!name) return fResultsContainer.back();
    std::string searchName(name);
    if (searchName.empty()) return fResultsContainer.back();
    for (AlgorithmResults::const_reverse_iterator s
             = fResultsContainer.rbegin();
         s != fResultsContainer.rend(); ++s) {
        std::string nm((*s)->GetName());
        if (nm == searchName) return *s;
    }
    return Cube::Handle<Cube::AlgorithmResult>();
}


void Cube::AlgorithmResult::AddObjectContainer(
    Cube::Handle<Cube::ReconObjectContainer> objects) {
    fObjectContainers.push_back(objects);
}

Cube::Handle<Cube::ReconObjectContainer>
Cube::AlgorithmResult::GetObjectContainer(
    const char* name) const {
    if (fObjectContainers.empty()) {
        return Cube::Handle<Cube::ReconObjectContainer>();
    }
    if (!name) return fObjectContainers.back();
    std::string searchName(name);
    for (ReconObjects::const_reverse_iterator s
             = fObjectContainers.rbegin();
         s != fObjectContainers.rend(); ++s) {
        std::string nm((*s)->GetName());
        if (nm == searchName) return *s;
    }
    return Cube::Handle<Cube::ReconObjectContainer>();
}


void Cube::AlgorithmResult::AddHitSelection(
    Cube::Handle<Cube::HitSelection> hits) {
    fHitSelections.push_back(hits);
}

Cube::Handle<Cube::HitSelection> Cube::AlgorithmResult::GetHitSelection(
    const char* name) const {
    if (fHitSelections.empty()) return Cube::Handle<Cube::HitSelection>();
    if (!name) return fHitSelections.back();
    std::string nm(name);
    for (HitSelections::const_reverse_iterator s = fHitSelections.rbegin();
         s != fHitSelections.rend(); ++s) {
        std::string sel((*s)->GetName());
        if (nm == sel) return *s;
    }
    return Cube::Handle<Cube::HitSelection>();
}

/// Print the hit information.
void Cube::AlgorithmResult::ls(Option_t *opt) const {
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
    // Print the class specific stuff.
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "Status: ";
    if (fStatusSummary != "") std::cout << fStatusSummary;
    else std::cout << "None";
    std::cout << std::endl;
    TROOT::IndentLevel();
    std::cout << "Hit Selections: " << fHitSelections.size() << std::endl;
    TROOT::IncreaseDirLevel();
    for (auto h : fHitSelections) {
        if (option.find("quiet") != std::string::npos) {
            TROOT::IndentLevel();
            if (h) {
                std::cout << h->ClassName() << ":: " << h->GetName()
                          << " with " << h->size() << " hits"
                          << std::endl;
            }
        }
        else h.ls(opt);
    }
    TROOT::DecreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "Recon Object Containers: " << fObjectContainers.size()
              << std::endl;
    TROOT::IncreaseDirLevel();
    for (auto h : fObjectContainers) {
        if (option.find("quiet") != std::string::npos) {
            TROOT::IndentLevel();
            if (h) {
                std::cout << h->ClassName() << ":: " << h->GetName()
                          << " with " << h->size() << " entries"
                          << std::endl;
            }
        }
        else h.ls(opt);
    }
    TROOT::DecreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "Algorithm Results: " << fResultsContainer.size()
              << std::endl;
    TROOT::IncreaseDirLevel();
    for (auto h : fResultsContainer) {
        h.ls(opt);
    }
    TROOT::DecreaseDirLevel();
    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
