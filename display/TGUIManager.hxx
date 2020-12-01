#ifndef TGUIManager_hxx_seen
#define TGUIManager_hxx_seen

#include <TGButton.h>
#include <TGListBox.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>

namespace Cube {
    class TGUIManager;
}

/// This creates the GUI interface for the event display, and then provides
/// handles so that other functions can connect to the interface.  A pointer
/// to each interface widget can be accessed by a method named by it's
/// function.  For instance, the "Read Next Event" button is named
/// "GetNextEventButton()", and returns a TGButton pointer.  This is managed
/// by the Cube::TEventDisplay class which returns it via the GetGUI() method.
///
/// \note This does not connect any actions to the buttons.  It only
/// constructs the GUI and then allows access to the GUI widget.
///
/// \note The default state of widgets (e.g. the button widgets) is set when
/// the widget is created.  This is how the default drawing attributes can be
/// set.
class Cube::TGUIManager {
public:
    /// Actually construct all the GUI.
    TGUIManager();
    ~TGUIManager();

    //////////////////////
    // CONTROL TAB WIDGETS
    //////////////////////

    /// Get the next event button widget.
    TGButton* GetNextEventButton() {return fNextEventButton;}

    /// Get the redraw current event button widget.
    TGButton* GetDrawEventButton() {return fDrawEventButton;}

    /// Get the previous event button widget.
    TGButton* GetPrevEventButton() {return fPrevEventButton;}

    // Get specific event
    TGNumberEntry* GetEventField() {return fInputEvent;}

    /////////////////////////////////////////////////////
    // Buttons to control showing truth information.
    /////////////////////////////////////////////////////

    /// Get the check button selecting if trajectories should be shown.
    TGButton* GetShowTrajectoriesButton() {return fShowTrajectoriesButton;}

    /// Get the check button selecting if G4 hits should be shown.
    TGButton* GetShowG4HitsButton() {return fShowG4HitsButton;}

    /////////////////////////////////////////////////////
    // Buttons to control showing the hits.
    /////////////////////////////////////////////////////

    /// If on, draw the simple hits (i.e. things like fibers)
    TGButton* GetShowSimpleHitsButton() {return fShowSimpleHitsButton;}

    /// If on, draw the composite hits (i.e. things like cubes of voxels).
    TGButton* GetShowCompHitsButton() {return fShowCompHitsButton;}

    /// If on, draw the main hits.
    TGButton* GetShowHitsButton() {return fShowHitsButton;}

    /// Get the check button selecting if all the hits in the event should be
    /// shown.  This draws all the hit selections in ~/hits/.
    TGButton* GetShowOtherHitsButton() {return fShowOtherHitsButton;}

    /// Get the check button selecting if all of the fiber hits in the event
    /// should be shown.  If ~/hits/sfg_3d is drawn, this controls if the hits
    /// for fibers are drawn (but never draw the lines for sfg_3d).
    TGButton* GetShowFibersButton() {return fShowFibersButton;}

    /// Get the check button selecting if reconstructed objects are shown.
    TGButton* GetShowFitsButton() {return fShowFitsButton;}

    /// Get the check button selecting if reconstructed clusters are shown.
    TGButton* GetSkipFitVerticesButton() {return fSkipFitVerticesButton;}

    /// Get the check button selecting if reconstructed tracks are shown.
    TGButton* GetSkipFitTracksButton() {return fSkipFitTracksButton;}

    /// Get the check button selecting if reconstructed clusters are shown.
    TGButton* GetSkipFitClustersButton() {return fSkipFitClustersButton;}

    /// Get the check button selecting if hits associated with shown
    /// reconstruction objects should be drawn.
    TGButton* GetShowFitHitsButton() {return fShowFitHitsButton;}

    /// Get the check button selecting if a line shouble be drawn between the
    /// ends of a reconstruction object.
    TGButton* GetShowFitEndsButton() {return fShowFitEndsButton;}

    /// Get the check button selecting if directions at each node should be
    /// drawn for objects that have directions.
    TGButton* GetShowFitDirectionButton() {return fShowFitDirectionButton;}

    /// Get the button to show constituent clusters of an object or not.
    TGButton* GetShowConstituentClustersButton()
        {return fShowConstituentClustersButton;}

    /// Get the button to select if hits are drawn for a cluster.  This is
    /// almost redundent with GetShowFitHitsButton, but it only draws the
    /// hits for the clusters.  This is good when you want to look at the
    /// constituent clusters and see which hits are there.
    TGButton* GetShowClusterHitsButton() {return fShowClusterHitsButton;}

    /// Get the check button selecting if a line shouble be drawn between the
    /// each of the hits of a cluster.
    TGButton* GetShowClusterConnectionsButton()
        {return fShowClusterConnectionsButton;}

    /// Get the button to select whether the cluster uncertainty, or cluster
    /// size is draw.
    TGButton* GetShowClusterUncertaintyButton()
        {return fShowClusterUncertaintyButton;}

    /// Get the check button selecting if G4 hits should be shown.
    TGButton* GetShowOnlySFGButton() {return fShowOnlySFGButton;}

    /// Get the button to draw the time and charge summary histograms.
    TGButton* GetDrawTimeChargeButton() {return fDrawTimeChargeButton;}

    //////////////////////
    // RESULTS TAB WIDGETS
    //////////////////////

    /// Get the list box with the select results to show.
    TGListBox* GetResultsList() {return fResultsList;}

    /// The get text entry widget for the default result to show.
    TGTextEntry* GetDefaultResult() {return fDefaultResult;}

private:

    /// Make a tab in the browser for control buttons.
    void MakeControlTab();

    ////////////////////////////////////////
    // Widgets in the control tab.
    ////////////////////////////////////////
    TGButton* fNextEventButton;
    TGButton* fDrawEventButton;
    TGButton* fPrevEventButton;
    TGNumberEntry* fInputEvent;

    TGButton* fShowTrajectoriesButton;
    TGButton* fShowG4HitsButton;

    TGButton* fShowSimpleHitsButton;
    TGButton* fShowCompHitsButton;
    TGButton* fShowHitsButton;

    TGButton* fShowOtherHitsButton;
    TGButton* fShowFibersButton;
    TGButton* fShowFitsButton;
    TGButton* fSkipFitVerticesButton;
    TGButton* fSkipFitTracksButton;
    TGButton* fSkipFitClustersButton;
    TGButton* fShowFitHitsButton;
    TGButton* fShowFitEndsButton;
    TGButton* fShowFitDirectionButton;
    TGButton* fShowConstituentClustersButton;
    TGButton* fShowClusterHitsButton;
    TGButton* fShowClusterConnectionsButton;
    TGButton* fShowClusterUncertaintyButton;
    TGButton* fShowOnlySFGButton;
    TGButton* fRecalculateViewButton;
    TGButton* fDrawHitButton;
    TGButton* fDrawTimeChargeButton;

    /// Make a tab in the browser to select algorithms shown.
    void MakeResultsTab();

    ////////////////////////////////////////
    // Widgets in the results tab.
    ////////////////////////////////////////

    /// A widget containing a list of possible TAlgorithmResult object names
    /// to be displayed (when selected).
    TGListBox* fResultsList;

    /// A regular expression to select the default result(s) to be selected.
    TGTextEntry* fDefaultResult;

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
