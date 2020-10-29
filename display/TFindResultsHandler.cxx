#include "TFindResultsHandler.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TGUIManager.hxx"

#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>

#include <TGeoManager.h>
#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>

#include <TPRegexp.h>

#include <sstream>

Cube::TFindResultsHandler::TFindResultsHandler() { }

Cube::TFindResultsHandler::~TFindResultsHandler() { }

void Cube::TFindResultsHandler::FillEntries(
    std::vector<std::string>& entries,
    std::string baseName,
    Cube::Handle<Cube::AlgorithmResult> current) {
    baseName += std::string(current->GetName()) + "/";
    for (Cube::AlgorithmResult::ReconObjects::reverse_iterator o
             = current->GetObjectContainers().rbegin();
         o != current->GetObjectContainers().rend(); ++o) {
        entries.push_back(baseName + (*o)->GetName());
    }
    // Add the daughter results to the stack
    for (Cube::AlgorithmResult::AlgorithmResults::reverse_iterator r
             = current->GetResultsContainer().rbegin();
         r !=  current->GetResultsContainer().rend(); ++r) {
        FillEntries(entries,baseName,*r);
    }
}

void Cube::TFindResultsHandler::Apply() {
    if (!Cube::gEvent) return;

    static int lastEvent = -1;
    if (lastEvent == Cube::gEvent->GetEventId()) return;
    lastEvent = Cube::gEvent->GetEventId();

    TGTextEntry* defResult
        = Cube::TEventDisplay::Get().GUI().GetDefaultResult();
    TGListBox* resultsList = Cube::TEventDisplay::Get().GUI().GetResultsList();

    std::string defaultResult(defResult->GetText());
    TPRegexp regularExp(defResult->GetText());

    resultsList->RemoveAll();
    int id = 0;
    // Forage the results...
    Cube::Handle<Cube::AlgorithmResult> topResult(Cube::gEvent,false);
    std::vector<std::string> existingEntries;
    std::string baseName;
    FillEntries(existingEntries,baseName,topResult);
    for(std::vector<std::string>::iterator e = existingEntries.begin();
        e != existingEntries.end(); ++e) {
        resultsList->AddEntry((*e).c_str(),++id);
        // Check to see if this result should be selected
        if (defaultResult.size() == 0) continue;
        if (!regularExp.Match((*e).c_str())) continue;
        resultsList->Select(id);
    }
    resultsList->Layout();
    resultsList->MapSubwindows();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
