#include "TGUIManager.hxx"
#include "TEventDisplay.hxx"

#include <TGFrame.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>

#include <TEveManager.h>
#include <TEveBrowser.h>

#include <TSystem.h>

Cube::TGUIManager::TGUIManager() {
    MakeResultsTab();
    MakeControlTab();
}

void Cube::TGUIManager::MakeControlTab() {
    TEveBrowser* browser = gEve->GetBrowser();

    // Define the generic layout.  The last four parameters are the padding
    // around the widgets.
    TGLayoutHints* layoutHints
        = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,
                            2, 2, 2, 2);

    // Embed a new frame in the event browser to have control buttons.
    browser->StartEmbedding(TRootBrowser::kLeft);
    TGMainFrame* mainFrame = new TGMainFrame(gClient->GetRoot(), 1000, 600);
    mainFrame->SetWindowName("Event Control");
    mainFrame->SetCleanup(kDeepCleanup);

    TGVerticalFrame* hf = new TGVerticalFrame(mainFrame);
    TGTextButton* textButton;

    // The general action buttons.
    textButton = new TGTextButton(hf, "Previous Event");
    textButton->SetToolTipText("Go to previous event.");
    textButton->SetTextJustify(36);
    textButton->SetMargins(0,0,0,0);
    textButton->SetWrapLength(-1);
    hf->AddFrame(textButton, layoutHints);
    fPrevEventButton = textButton;

    textButton = new TGTextButton(hf, "Redraw Event");
    textButton->SetToolTipText("Refresh the current view.");
    textButton->SetTextJustify(36);
    textButton->SetMargins(0,0,0,0);
    textButton->SetWrapLength(-1);
    hf->AddFrame(textButton, layoutHints);
    fDrawEventButton = textButton;

    textButton = new TGTextButton(hf, "Next Event");
    textButton->SetToolTipText("Go to previous event.");
    textButton->SetTextJustify(36);
    textButton->SetMargins(0,0,0,0);
    textButton->SetWrapLength(-1);
    hf->AddFrame(textButton, layoutHints);
    fNextEventButton = textButton;

    TGGroupFrame *groupFrame = new TGGroupFrame(hf, "Event Number");
    fInputEvent = new TGNumberEntry(
        groupFrame, 0, 9,999,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        0, 99999);
    groupFrame->AddFrame(fInputEvent,layoutHints);
    hf->AddFrame(groupFrame, layoutHints);

    /////////////////////
    // Buttons to control how truth objects are drawn.
    /////////////////////
    TGCheckButton *checkButton;

    checkButton = new TGCheckButton(hf,"Show Trajectories");
    checkButton->SetToolTipText(
        "Show the GEANT4 trajectories and trajectory points.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowTrajectoriesButton = checkButton;

    checkButton = new TGCheckButton(hf,"Show G4 Hits");
    checkButton->SetToolTipText(
        "Show the GEANT4 hits.  This shows the energy deposition in the "
        "defined sensitive detectors.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    checkButton->SetOn();
    hf->AddFrame(checkButton, layoutHints);
    fShowG4HitsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Recalculate view");
    checkButton->SetToolTipText(
        "Recalculate center of rotation based.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    checkButton->SetOn();
    hf->AddFrame(checkButton, layoutHints);
    fRecalculateViewButton = checkButton;

    //////////////////////////
    // Buttons to control the drawing of hits.
    //////////////////////////

    checkButton = new TGCheckButton(hf,"Show hits");
    checkButton->SetToolTipText(
        "Show top level hits, modified by later selections.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowHitsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Select simple hits");
    checkButton->SetToolTipText(
        "Select simple hits to be drawn.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowSimpleHitsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Select composite hits");
    checkButton->SetToolTipText(
        "Select composite hits to be drawn.");
    checkButton->SetOn();
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowCompHitsButton = checkButton;

    groupFrame = new TGGroupFrame(hf, "Hit charge threshold");
    fHitChargeThreshold = new TGNumberEntry(
        groupFrame, 10.0, 9, 5,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMin,
        0, 1E+20);
    groupFrame->AddFrame(fHitChargeThreshold,layoutHints);
    hf->AddFrame(groupFrame, layoutHints);

    //////////////////////////////////////////////////////
    // Control how fits are drawn
    //////////////////////////////////////////////////////
    checkButton = new TGCheckButton(hf,"Show reconstruction objects");
    checkButton->SetToolTipText(
        "Show the reconstruction objects in the selected container.");
    checkButton->SetOn();
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowFitsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Skip reconstructed vertices.");
    checkButton->SetToolTipText(
        "Skip any vertices in the reconstructed objects.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fSkipFitVerticesButton = checkButton;

    checkButton = new TGCheckButton(hf,"Skip reconstructed clusters.");
    checkButton->SetToolTipText(
        "Skip any clusters in the reconstructed objects.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fSkipFitClustersButton = checkButton;

    checkButton = new TGCheckButton(hf,"Skip long reconstructed tracks.");
    checkButton->SetToolTipText(
        "Skip any long tracks in the reconstructed objects.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fSkipFitTracksButton = checkButton;

    groupFrame = new TGGroupFrame(hf, "Nodes in Long Track");
    fSkipFitTrackNodes = new TGNumberEntry(
        groupFrame, 0, 9,999,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        0, 100);
    groupFrame->AddFrame(fSkipFitTrackNodes,layoutHints);
    hf->AddFrame(groupFrame, layoutHints);

    checkButton = new TGCheckButton(hf,"Show reconstruction object hits");
    checkButton->SetToolTipText(
        "Show the hits in the reconstruction objects.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowFitHitsButton = checkButton;

    //////////////////////////////////////////////////////////////////
    // Control how clusters are drawn
    //////////////////////////////////////////////////////////////////
    checkButton = new TGCheckButton(hf,"Connect first and last cluster hit.");
    checkButton->SetToolTipText(
        "Draw a line connecting the first and last hit in the clusters.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowFitEndsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Show connections in clusters.");
    checkButton->SetToolTipText(
        "Show the order of hits in the clusters.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowClusterConnectionsButton = checkButton;

    checkButton = new TGCheckButton(hf,"Show cluster position uncertainty.");
    checkButton->SetToolTipText(
        "Show cluster centroid uncertainty, not the cluster size.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowClusterUncertaintyButton = checkButton;

    checkButton = new TGCheckButton(hf,"Show cluster hits.");
    checkButton->SetToolTipText(
        "Show cluster hits (almost the same as showing all hits.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowClusterHitsButton = checkButton;

    //////////////////////////////////////////////////////////////////
    // Control how tracks and showers are drawn
    //////////////////////////////////////////////////////////////////
    checkButton = new TGCheckButton(hf,"Show object directions.");
    checkButton->SetToolTipText(
        "Show the directions of tracks and nodes.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowFitDirectionButton = checkButton;

    checkButton = new TGCheckButton(hf,"Show object constiuents.");
    checkButton->SetToolTipText(
        "Show the node objects for tracks and showers.");
    checkButton->SetTextJustify(36);
    checkButton->SetMargins(0,0,0,0);
    checkButton->SetWrapLength(-1);
    hf->AddFrame(checkButton, layoutHints);
    fShowConstituentClustersButton = checkButton;

    /////////////////////////////////////////////////////////////
    // Do the final layout and mapping.
    /////////////////////////////////////////////////////////////
    mainFrame->AddFrame(hf, layoutHints);
    mainFrame->MapSubwindows();
    mainFrame->Resize();
    mainFrame->MapWindow();
    browser->StopEmbedding();
    browser->SetTabTitle("Control", 0);
}

void Cube::TGUIManager::MakeResultsTab() {
    TEveBrowser* browser = gEve->GetBrowser();

    // Embed a new frame in the event browser to select fit objects.
    browser->StartEmbedding(TRootBrowser::kLeft);
    TGMainFrame* mainFrame = new TGMainFrame(gClient->GetRoot(), 1000, 600);
    mainFrame->SetWindowName("Fit Selection");
    mainFrame->SetCleanup(kDeepCleanup);

    TGVerticalFrame* hf = new TGVerticalFrame(mainFrame);

    // Create the listbox for the results.
    fResultsList = new TGListBox(hf);
    fResultsList->SetMultipleSelections(true);
    TGLayoutHints* layoutList = new TGLayoutHints(kLHintsLeft
                                                  | kLHintsTop
                                                  | kLHintsExpandX
                                                  | kLHintsExpandY);
    hf->AddFrame(fResultsList,layoutList);

    // Create a text entry button to choose a default result.
    TGLabel* txt = new TGLabel(hf,"Default Result");
    TGLayoutHints* layoutHints = new TGLayoutHints(kLHintsLeft
                                                   | kLHintsTop
                                                   | kLHintsExpandX,
                                                   2,2,2,2);
    hf->AddFrame(txt,layoutHints);

    // Create a text entry button to choose a default result.  This also sets
    // the default default result.  It should be possible to override this
    // from the command line!
    fDefaultResult = new TGTextEntry(hf);
    fDefaultResult->SetText("CubeEvent/final");
    fDefaultResult->SetToolTipText(
        "Enter a regular expression to select algorithm results\n"
        "to be shown.  The usual regexp syntax applys:\n"
        "    .   -- Match any character\n"
        "    .*  -- Match any sequence of characters\n"
        "    |   -- \'or\'\n"
        "See grep for more details.  A result is show if the\n"
        "regexpn matchs any sub-string in the result name." );

    hf->AddFrame(fDefaultResult,layoutHints);

    // Do the final layout and mapping.
    TGLayoutHints* layoutFrame
        = new TGLayoutHints(kLHintsLeft
                            | kLHintsTop
                            | kLHintsExpandX
                            | kLHintsExpandY,
                            2, 2, 2, 2);
    mainFrame->AddFrame(hf, layoutFrame);
    mainFrame->MapSubwindows();
    mainFrame->Resize();
    mainFrame->MapWindow();
    browser->StopEmbedding();
    browser->SetTabTitle("Recon", 0);
}

Cube::TGUIManager::~TGUIManager() { }

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
