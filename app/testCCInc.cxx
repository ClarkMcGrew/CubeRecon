#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconVertex.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolInitialize.hxx>
#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>
#include <ToolRecon.hxx>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGeoManager.h>
#include <TDatabasePDG.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <cmath>

typedef std::deque<Cube::Handle<Cube::ReconTrack>> MatchedTracks;
namespace {
    struct CMPMatch {
        bool operator () (const MatchedTracks& lhs,
                          const MatchedTracks& rhs) {
            return (lhs.front()->GetPosition().T()
                    < rhs.front()->GetPosition().T());
        }
    };

    struct ORDERMatch {
        bool operator () (const Cube::Handle<Cube::ReconTrack>& lhs,
                          const Cube::Handle<Cube::ReconTrack>& rhs) {
            std::string lhsDetectors = Cube::Tool::ObjectDetectors(*lhs);
            std::string rhsDetectors = Cube::Tool::ObjectDetectors(*rhs);
            if (rhsDetectors.find(":TPC") != std::string::npos
                && lhsDetectors.find(":TPC") != std::string::npos) {
                throw std::runtime_error("Can't have two TPC segments");
            }
            if (lhsDetectors.find(":TPC") != std::string::npos) return false;
            if (rhsDetectors.find(":TPC") != std::string::npos) return true;
            return lhs->GetPosition().T() < rhs->GetPosition().T();
        }
    };
}


bool histInitialized = false;
TH1F* histMatchDistance = NULL;
TH1F* histMatchAngle = NULL;
TH1F* histMuonTruth = NULL;
TH1F* histMuonTruthX = NULL;
TH1F* histMuonTruthY = NULL;
TH1F* histMuonTruthZ = NULL;
TH1F* histMissedMomentum = NULL;
TH1F* histMissedCosine = NULL;
TH1F* histMissedPhi = NULL;
TH1F* histMissedZ = NULL;
TH1F* histMuonReco = NULL;
TH1F* histBackgroundReco = NULL;

/// Assign ECal hit segments to ecal hits.  This is a test bed for later
/// analysis.
bool AnalyzeEvent(Cube::Event& event) {
    /// Flag tht this event should be saved (if not rejected for some cut).
    bool shouldSave = false;

    // The sign of muon that should be selected: neutrino is -1.0, and
    // anti-neutrino is 1.0
    const double signSelection = -1.0;
    const double fiducialCut = 10*unit::cm;

    if (!histInitialized) {
        Cube::Tool::Initialize();
        const int momentumBins = 20;
        const int momentumMax = 250.0*unit::MeV*momentumBins;
        histInitialized = true;
        histMatchDistance
            = new TH1F("matchDistance",
                       "Intersection distance for 3DST and TPC tracks",
                       50, 0.0, 25.0);
        histMatchAngle
            = new TH1F("matchAngle",
                       "Intersection angle for 3DST and TPC tracks",
                       50, 0.0, 0.50);
        histMuonTruth
            = new TH1F("muonTruth",
                       "True muon momentum of contained interactions",
                       momentumBins, 0.0, momentumMax);
        histMuonTruthX
            = new TH1F("muonTruthX",
                       "True X vertex",
                       600, -3000.0, 3000.0);
        histMuonTruthY
            = new TH1F("muonTruthY",
                       "True Y vertex",
                       600, -6000.0, 0.0);
        histMuonTruthZ
            = new TH1F("muonTruthZ",
                       "True Z vertex",
                       600, 61000.0, 67000.0);

        histMissedMomentum
            = new TH1F("missedMomentum",
                       "True muon momentum of missed interactions",
                       momentumBins, 0.0, momentumMax);

        histMissedCosine
            = new TH1F("missedCosine",
                       "True muon direction cosine of missed interactions",
                       50, -1.0, 1.0);

        histMissedPhi
            = new TH1F("missedPhi",
                       "True muon direction phi angle of missed interactions",
                       72, -180.0, 180.0);

        histMissedZ
            = new TH1F("missedZ",
                       "True muon vertex of missed interactions",
                       50, 61000.0, 67000.0);

        histMuonReco
            = new TH1F("muonReco",
                       "Reconstructed muon momentum for contained events",
                       momentumBins, 0.0, momentumMax);
        histBackgroundReco
            = new TH1F("backgroundReco",
                       "Reconstructed background momentum for contained events",
                       momentumBins, 0.0, momentumMax);
   }

    Cube::Event::G4TrajectoryContainer& trajectories = event.G4Trajectories;
    Cube::Handle<Cube::HitSelection> hits = event.GetHitSelection();
    if (!hits) return false;
    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return false;

    // Get the PDG database.
    const TDatabasePDG *pdgData = TDatabasePDG::Instance();

    // Collect information about the interactions in the 3DST.
    std::vector<int> allPrimaries = Cube::Tool::AllPrimaries(event);
    std::vector<int> signalPrimaries;
    std::vector<Cube::Handle<Cube::G4Trajectory>> signalTrajectories;
    for (int prim : allPrimaries) {
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(prim);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO PRIMARY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> traj = t->second;
        // Check the particle type.
        int pdgCode = - traj->GetPDGCode() * signSelection;
        if (pdgCode != 13) continue;

        // Check the fiducial volume.
        double fidTruth
            = Cube::Tool::Contained3DSTPoint(
                traj->GetInitialPosition().Vect());
        if (fidTruth < fiducialCut) continue;

        // This is a CC interaction inside the 3DST fiducial volume.
        signalPrimaries.push_back(prim);
        signalTrajectories.push_back(traj);

        // Get the momentum
        double momTruth = traj->GetInitialMomentum().P();
        histMuonTruth->Fill(momTruth);
        histMuonTruthX->Fill(traj->GetInitialPosition().X());
        histMuonTruthY->Fill(traj->GetInitialPosition().Y());
        histMuonTruthZ->Fill(traj->GetInitialPosition().Z());
        std::cout << "Truth "
                  << " " << pdgCode
                  << " " << traj->GetPDGCode()
                  << " " << momTruth
                  << " " << fidTruth
                  << " " << traj->GetInitialPosition().T()
                  << std::endl;
    }

    // Cache all of the main trajectory results
    std::map<Cube::Handle<Cube::ReconTrack>,int> mainTrajectoryCache;
    for (Cube::Handle<Cube::ReconObject> obj : *objects) {
        Cube::Handle<Cube::ReconTrack> track = obj;
        if (!track) continue;
        int mainTrajectory = Cube::Tool::MainTrajectory(event,*obj);
        if (mainTrajectory<0) continue;
        mainTrajectoryCache[track] = mainTrajectory+1;
    }

    // Collect the TPC tracks.
    Cube::ReconObjectContainer tracksTPC;
    for (Cube::Handle<Cube::ReconObject> obj : *objects) {
        Cube::Handle<Cube::ReconTrack> track = obj;
        if (!track) continue;
        std::string detectors = Cube::Tool::ObjectDetectors(*obj);
        if (detectors.find(":TPC") == std::string::npos) continue;
        if (detectors.find(":3DST") != std::string::npos) {
            std::cout << "3DST in TPC Track" << std::endl;
            continue;
        }
        tracksTPC.push_back(track);
    }

    // Collect the 3DST clusters tracks and vertices
    Cube::ReconObjectContainer clusters3DST;
    Cube::ReconObjectContainer tracks3DST;
    Cube::ReconObjectContainer vertex3DST;
    for (Cube::Handle<Cube::ReconObject> obj : *objects) {
        std::string detectors = Cube::Tool::ObjectDetectors(*obj);
        if (detectors.find(":3DST") == std::string::npos) continue;
        if (detectors.find(":TPC") != std::string::npos) {
            std::cout << "TPC in 3DST track" << std::endl;
            continue;
        }
        Cube::Handle<Cube::ReconCluster> cluster = obj;
        if (cluster) {
            clusters3DST.push_back(cluster);
            continue;
        }
        Cube::Handle<Cube::ReconTrack> track = obj;
        if (track) {
            tracks3DST.push_back(track);
            continue;
        }
        Cube::Handle<Cube::ReconVertex> vertex = obj;
        if (vertex) {
            vertex3DST.push_back(vertex);
            continue;
        }
    }

    std::vector<Cube::Handle<Cube::ReconObject>> matchCand3DST;
    for (Cube::Handle<Cube::ReconObject> dst : tracks3DST) {
        Cube::Handle<Cube::ReconTrack> track3DST = dst;
        int trajId3DST = mainTrajectoryCache[track3DST] - 1;
        if (trajId3DST < 0) continue;
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(trajId3DST);
        if (t == event.G4Trajectories.end()) {
            std::cout << "3DST TRACK NO TRAJECTORY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> traj3DST = t->second;
        // 3DST electron id is excellent.
        if (traj3DST->GetPDGCode() == 11) continue;
        if (traj3DST->GetPDGCode() == -11) continue;
        // 3DST low momentum proton id is excellent.
        if (traj3DST->GetPDGCode() == 2212
            && traj3DST->GetInitialMomentum().P() < 900*unit::MeV) continue;
        // Check that the 3DST track is long enough.
        double len = Cube::Tool::TrackLength(*track3DST);
        if (len < 0.9*fiducialCut) continue;
        // Check that the 3DST track exits the 3DST
        double fidDist = Cube::Tool::Contained3DSTObject(*dst);
        if (fidDist > 25.0*unit::mm) continue;
        matchCand3DST.push_back(dst);
    }

    std::vector<MatchedTracks> matchedTracks;
    for (Cube::Handle<Cube::ReconObject> tpc : tracksTPC) {
        Cube::Handle<Cube::ReconTrack> trackTPC = tpc;
        int trajIdTPC = mainTrajectoryCache[trackTPC]-1;
        if (trajIdTPC < 0) continue;
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(trajIdTPC);
        if (t == event.G4Trajectories.end()) {
            std::cout << "TPC TRACK NO TRAJECTORY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> trajTPC = t->second;
#define ASSUME_TPC_ELECTRON_ID
#ifdef ASSUME_TPC_ELECTRON_ID
        // Reject electrons.  This is cheating, but the TPC id is very good.
        if (trajTPC->GetPDGCode() == 11) continue;
        if (trajTPC->GetPDGCode() == -11) continue;
#endif
#define ASSUME_TPC_PROTON_ID
#ifdef ASSUME_TPC_PROTON_ID
        // TPC Proton id is excellent
        if (trajTPC->GetPDGCode() == 2212
            && trajTPC->GetInitialMomentum().P() < 1.1*unit::GeV) continue;
#endif
        // Get the drift velocity for this track.  The sign is the
        // important part.
        Cube::Handle<Cube::Hit> hit = trackTPC->GetHitSelection()->front();
        double driftVelocity = hit->GetProperty("DriftVelocity");
        if (std::abs(driftVelocity) < 1E-8) {
            throw std::runtime_error("inconceivable");
        }
        double bestDist = 1E+20;
        double bestCos = 1E+20;
        Cube::Handle<Cube::ReconTrack> best3DST;
        for (Cube::Handle<Cube::ReconObject> dst : matchCand3DST) {
            Cube::Handle<Cube::ReconTrack> track3DST = dst;
            // The TPC track time should always be zero.
            double dt = track3DST->GetBack()->GetPosition().T()
                - trackTPC->GetPosition().T();
            TVector3 backDST = track3DST->GetBack()->GetPosition().Vect();
            TVector3 dirDST =  track3DST->GetBack()->GetDirection().Unit();
            // Correct the TPC track position to the time of this 3DST track.
            // The TPC track might be oriented in either directions, so do
            // both the front and back.
            TVector3 frontTPC = trackTPC->GetFront()->GetPosition().Vect();
            TVector3 frontDirTPC = trackTPC->GetFront()->GetDirection().Unit();
            TVector3 backTPC = trackTPC->GetBack()->GetPosition().Vect();
            TVector3 backDirTPC = trackTPC->GetBack()->GetDirection().Unit();
            frontTPC.SetX(frontTPC.X() + dt*driftVelocity);
            backTPC.SetX(backTPC.X() + dt*driftVelocity);
            // Figure out which end of the TPC track matchs to the 3DST track
            double frontDist = (frontTPC-backDST).Mag();
            double backDist = (backTPC-backDST).Mag();
            // Make a rather horrid linear approximation to see if the tracks
            // intersect.  This won't work well for side exiting 3DST tracks
            // that go quite aways before hitting the TPC, but is OK for most
            // of the phase space where the tracks "immediately" enter the
            // tpc.
            double dCos;
            TVector3 dMiss;
            if (frontDist < backDist) {
                dCos = std::abs(dirDST*frontDirTPC);
                double moveDist = (frontTPC-backDST)*dirDST;
                dMiss = frontTPC - (backDST + moveDist * dirDST);
            }
            else {
                dCos = std::abs(dirDST*backDirTPC);
                double moveDist = (backTPC-backDST)*dirDST;
                dMiss = backTPC - (backDST + moveDist * dirDST);
            }
            // Tracks need to be in about the same direction!
            if (dCos < 0.7) continue;
            if (bestDist < dMiss.Mag()) continue;
#define PID_THREE_DST
#ifdef PID_THREE_DST
            // Might be a match, so check the PID.
            int trajId3DST = mainTrajectoryCache[track3DST]-1;
            if (trajId3DST < 0) continue;
            Cube::Event::G4TrajectoryContainer::iterator
                t = event.G4Trajectories.find(trajId3DST);
            if (t == event.G4Trajectories.end()) {
                std::cout << "3DST TRACK NO TRAJECTORY!!!!!" << std::endl;
                continue;
            }
            Cube::Handle<Cube::G4Trajectory> traj3DST = t->second;
            // 3DST electron id is excellent.
            if (traj3DST->GetPDGCode() == 11) continue;
            if (traj3DST->GetPDGCode() == -11) continue;
            // 3DST low momentum proton id is excellent.
            if (traj3DST->GetPDGCode() == 2212
                && traj3DST->GetInitialMomentum().P() < 900*unit::MeV) continue;
#endif
            bestDist = dMiss.Mag();
            best3DST = track3DST;
            bestCos = dCos;
        }
        histMatchDistance->Fill(bestDist);
        histMatchAngle->Fill(std::acos(bestCos));
        if (bestDist < 15*unit::mm) {
            MatchedTracks tmp;
            tmp.push_back(best3DST);
            tmp.push_back(trackTPC);
            matchedTracks.push_back(tmp);
        }
    }

    // Now check if there are 3DST tracks that need to be merged.  This cheats
    // and assumes we will be able to perfectly match track segments.  The
    // 3DST segments are ordered by the times
    for (MatchedTracks& match : matchedTracks) {
        Cube::Handle<Cube::ReconTrack> track3DST = match.front();
        double bestDist = 1E+20;
        double bestCos = 1E+20;
        TVector3 posDST = track3DST->GetFront()->GetPosition().Vect();
        TVector3 dirDST = track3DST->GetFront()->GetDirection();
        int trajId3DST = mainTrajectoryCache[track3DST]-1;
        double bestDistance = 1E+20;
        for (Cube::Handle<Cube::ReconObject> dst : tracks3DST) {
            Cube::Handle<Cube::ReconTrack> cand = dst;
            // Only look at tracks
            if (!cand) continue;
            // Make sure the candidate isn't the same as the track.
            if (track3DST == cand) continue;
            int trajIdCand = mainTrajectoryCache[cand]-1;
            if (trajIdCand != trajId3DST) continue;
            match.push_front(cand);
        }
        std::sort(match.begin(),match.end(),ORDERMatch());
    }

    // At this point, the elements of matchedTracks will contain a deque of
    // tracks that are ordered along the particle path.  The front() track is
    // where the particle is reconstructed as starting.  The back() track is
    // the TPC track for the momentum.

    // Sort the matches tracks by time.
    std::sort(matchedTracks.begin(), matchedTracks.end(),CMPMatch());

    // Loop over the matched tracks, identify the time slices, and then pick
    // the highest energy mip track in the slice as the mip candidates.  At
    // the end of this the mipCandidates vector will contain all of the
    // (anti)neutrino candidate tracks.
    MatchedTracks bestCandidate;
    double lastTime = -1E+20;
    double bestMomentum = 0;
    std::vector<MatchedTracks> mipCandidates;
    for (MatchedTracks match : matchedTracks) {
        // Check if the new match is starting a candidate interaction vertex.
        // It's a new "interaction" if the match is more than 20 ns from the
        // current best match. (Remember, the matchs are sorted by time).  The
        // 20 ns cut is based on the size of the 3DST.
        double deltaTime = match.front()->GetPosition().T() - lastTime;
        if (lastTime < 0) {
            // lastTime is negative the first time through the lop.
            bestMomentum = 0.0;
            bestCandidate.clear();
        }
        else if (deltaTime > 20*unit::ns) {
            // If we have seen a candidate then save it.
            if (bestCandidate.size() > 0) {
                mipCandidates.push_back(bestCandidate);
            }
            bestMomentum = 0.0;
            bestCandidate.clear();
        }
        // The new track gives the new last time.
        lastTime = match.front()->GetPosition().T();

        // Get the truth information...
        int trajId3DST = mainTrajectoryCache[match.front()]-1;
        int trajIdTPC =  mainTrajectoryCache[match.back()]-1;
        if (trajIdTPC < 0) {
            std::cout << "BAD TPC TRACK MATCH" << std::endl;
            continue;
        }
        int primIdTPC = Cube::Tool::PrimaryId(event,trajIdTPC);
        if (primIdTPC < 0) {
            std::cout << "BAD TPC PRIMARY" << std::endl;
            continue;
        }

        // Get the TPC trajectory.
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(trajIdTPC);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO TPC TRAJECTORY" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> trajTPC = t->second;

        // Get the TPC primary.
        t = event.G4Trajectories.find(primIdTPC);
        if (t == event.G4Trajectories.end()) {
            std::cout << "NO TPC PRIMARY" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> primTPC = t->second;


        int pdgCode =  trajTPC->GetPDGCode();

        // Check the sign of the TPC Track.  This is cheating, but the sign is
        // VERY well determined from curvature.
        const TParticlePDG *pdgParticle = pdgData->GetParticle(pdgCode);
        double tpcCharge = 0.0;
        if (pdgParticle) tpcCharge = pdgParticle->Charge();

        ////////////////////////////////////////////////////////////
        // Check the sign of the track.  Reject wrong sign tracks.  The
        // signSelection value is set at the top of this routine.
        ////////////////////////////////////////////////////////////
        if (signSelection*tpcCharge < 0.0 || std::abs(tpcCharge) < 1.0) {
            continue;
        }

        ///////////////////////////////////////////////////////////.
        // Apply a fiducial volume cut
        ///////////////////////////////////////////////////////////.
        double fidDist
            = Cube::Tool::Contained3DSTPoint(
                match.front()->GetFront()->GetPosition().Vect());
        if (fidDist < fiducialCut) {
            continue;
        }

#define VETO_UPSTREAM_ACTIVITY
#ifdef VETO_UPSTREAM_ACTIVITY
        //////////////////////////////////////////////////////////////////
        // Veto tracks that have coincident activity on the upstream face of
        // the 3DST.
        //////////////////////////////////////////////////////////////////
        double eventTime = match.front()->GetFront()->GetPosition().T();
        double upstreamDist = 1E+20;
        const double externalTimeCut = 10*unit::ns;
        for (Cube::Handle<Cube::ReconObject> obj : clusters3DST) {
            Cube::Handle<Cube::ReconCluster> cluster = obj;
            double dT = std::abs(cluster->GetPosition().T() - eventTime);
            if (dT > externalTimeCut) continue;
            for (Cube::Handle<Cube::Hit> hit : *obj->GetHitSelection()) {
                upstreamDist = std::min(upstreamDist,
                                        Cube::Tool::Contained3DSTUpstream(
                                            hit->GetPosition()));
            }
        }
        for (Cube::Handle<Cube::ReconTrack> track : tracks3DST) {
            double dT = std::abs(track->GetPosition().T() - eventTime);
            if (dT > externalTimeCut) continue;
            for (Cube::Handle<Cube::Hit> hit : *track->GetHitSelection()) {
                upstreamDist = std::min(upstreamDist,
                                        Cube::Tool::Contained3DSTUpstream(
                                            hit->GetPosition()));
            }
        }
        if (upstreamDist < 5*unit::cm) {
            continue;
        }
#endif

        // Get the "reconstructed" TPC momentum.
        double tpcMom = trajTPC->GetInitialMomentum().P();

        //////////////////////////////////////////////////////////////////
        // Take the highest momentum candidate in the time slice.
        //////////////////////////////////////////////////////////////////
        if (tpcMom <= bestMomentum) {
            continue;
        }

        // If we get here, then the current set of matched tracks will be the
        // best candidate for the current time slice.
#define LOUD_AND_PROUD
#ifdef LOUD_AND_PROUD
        std::cout << "CANDIDATE " << match.front()->GetPosition().T()
                  << " T " << trajTPC->GetInitialPosition().T()
                  << " id " << trajId3DST
                  << " " << trajIdTPC
                  << " " << primIdTPC
                  << " P " << tpcMom << " -> " << bestMomentum
                  << " C " << trajTPC->GetPDGCode()
                  << " " << tpcCharge
                  << std::endl;
#endif
        bestCandidate = match;
        bestMomentum = tpcMom;
    }

    // Save the last candidate (if it exists).  Out of caution, this resets
    // bestCandidates and bestMomentum, but they shouldn't be used after this.
    if (bestCandidate.size() > 0) {
        mipCandidates.push_back(bestCandidate);
        bestCandidate.clear();
        bestMomentum = 0.0;
    }

    // Figure out which true interactions that were missed.  The
    // interactionMissed vector containes all of the muons where the
    // interaction was entirely missed.  The signalMissed vector containes the
    // muons where the interaction was found, but the wrong particle was
    // selected as the MIP candidate.[q
    std::vector<Cube::Handle<Cube::G4Trajectory>> signalMissed;
    std::vector<Cube::Handle<Cube::G4Trajectory>> interactionMissed;
    for (Cube::Handle<Cube::G4Trajectory> signal : signalTrajectories) {
        double  minDeltaT = 1E+20;
        bool matchedTrajectory = false;
        for (MatchedTracks cand : mipCandidates) {
            int trajId3DST = mainTrajectoryCache[cand.front()]-1;
            int primId3DST = Cube::Tool::PrimaryId(event,trajId3DST);
            Cube::Event::G4TrajectoryContainer::iterator
                t = event.G4Trajectories.find(primId3DST);
            if (t == event.G4Trajectories.end()) {
                std::cout << "MIP CANDIDATE WITHOUT PRIMARY!!" << std::endl;
                continue;
            }
            Cube::Handle<Cube::G4Trajectory> prim3DST = t->second;
            if (prim3DST == signal) {
                matchedTrajectory = true;
                minDeltaT = 0.0;
                break;
            }
            double dt = std::abs(signal->GetInitialPosition().T()
                                 - prim3DST->GetInitialPosition().T());
            minDeltaT = std::min(minDeltaT,dt);
        }
        if (matchedTrajectory) continue;
        if (minDeltaT > 1.0*unit::ns) {
            interactionMissed.push_back(signal);
            continue;
        }
        signalMissed.push_back(signal);
    }

    ////////////////////////////////////////////////////////////////////
    //
    // All of the information has been collected, now plot it.
    //
    ////////////////////////////////////////////////////////////////////


    // Now the mipCandidates vector has all of the selected CC candidates.
    // std::cout << "MIP Candidates " << mipCandidates.size() << std::endl;
    for (MatchedTracks cand : mipCandidates) {
        bool wrongParticle = false;
        bool wrongSign = false;
        bool externalInteraction = false;

        int trajId3DST = mainTrajectoryCache[cand.front()]-1;
        int trajIdTPC = mainTrajectoryCache[cand.back()]-1;
        Cube::Event::G4TrajectoryContainer::iterator
            t = event.G4Trajectories.find(trajIdTPC);
        if (t == event.G4Trajectories.end()) {
            std::cout << "MIP CANDIDATE WITH NO TRAJECTORY!!!!!" << std::endl;
            continue;
        }
        Cube::Handle<Cube::G4Trajectory> trajTPC = t->second;

        // Check the particle type.
        if (trajTPC->GetPDGCode() != 13 && trajTPC->GetPDGCode() != -13) {
            wrongParticle = true;
        }

        // Check the sign
        const TParticlePDG *pdgParticle
            = pdgData->GetParticle(trajTPC->GetPDGCode());
        double charge = 0.0;
        if (pdgParticle) charge = pdgParticle->Charge();
        if (signSelection*charge < 0.0 || std::abs(charge) < 1.0) {
            wrongSign = true;
        }

        // Check if the event entered the detector.
        double fidTruth
            = Cube::Tool::Contained3DSTPoint(
                trajTPC->GetInitialPosition().Vect());
        if (fidTruth < 0.0*unit::cm) externalInteraction = true;

        double momReco = trajTPC->GetInitialMomentum().P();

        bool signal = !externalInteraction && !wrongSign && !wrongParticle;
        if (signal) {
            histMuonReco->Fill(momReco);
            std::cout << "Sgnl " << trajTPC->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << " " << trajTPC->GetInitialPosition().T()
                      << std::endl;
            continue;
        }
        histBackgroundReco->Fill(momReco);

        if (externalInteraction) {
            std::cout << "Back External " << trajTPC->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << " " << trajTPC->GetInitialPosition().T()
                      << std::endl;
            continue;
        }

        if (wrongParticle) {
            std::cout << "Back Particle " << trajTPC->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << " " << trajTPC->GetInitialPosition().T()
                      << std::endl;
            continue;
        }

        if (wrongSign) {
            std::cout << "Back Sign " << trajTPC->GetPDGCode()
                      << " " << momReco
                      << " " << fidTruth
                      << " " << trajTPC->GetInitialPosition().T()
                      << std::endl;
            continue;
        }

    }

    for (Cube::Handle<Cube::G4Trajectory> traj : interactionMissed) {
        histMissedMomentum->Fill(traj->GetInitialMomentum().P());
        double cosZ = traj->GetInitialMomentum().Vect().Unit().Z();
        histMissedCosine->Fill(cosZ);
        histMissedZ->Fill(traj->GetInitialPosition().Z());
        // Get the angle relative to vertical.
        double phi = std::atan2(traj->GetInitialMomentum().X(),
                                traj->GetInitialMomentum().Y());
        while (phi < -180.0*unit::degree) phi += 360.0*unit::degree;
        while (phi > 180.0*unit::degree) phi -= 360.0*unit::degree;
        histMissedPhi->Fill(phi/unit::degree);
        std::cout << "Missed " << traj->GetPDGCode()
                  << " " << traj->GetInitialMomentum().P()
                  << " " << traj->GetInitialPosition().T()
                  << std::endl;
    }

    for (Cube::Handle<Cube::G4Trajectory> traj : signalMissed) {
        std::cout << "Wrong " << traj->GetPDGCode()
                  << " " << traj->GetInitialMomentum().P()
                  << " " << traj->GetInitialPosition().T()
                  << std::endl;
    }

    // Uncomment to save all of the events.
    // shouldSave = true;
    return shouldSave;
}

int main(int argc, char** argv) {
    int maxEntries = 1E+8; // Maximum to process.
    int firstEntry = 0;
    std::string outputName;

    while (true) {
        int c = getopt(argc,argv,"n:o:s:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
            break;
        }
        case 'o': {
            outputName = optarg;
            break;
        }
        case 's': {
            std::istringstream tmp(optarg);
            tmp >> firstEntry;
            break;
        }
        default: {
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
                      << "-o <number>  : Output file"
                      << std::endl
                      << "-s <number>  : Skip <number> entries"
                      << std::endl
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl;
            exit(1);
        }
        }
    }

    std::vector<std::string> inputNames;
    if (argc <= optind) throw std::runtime_error("Missing input file");
    while (optind < argc) {
        inputNames.push_back(argv[optind++]);
        std::cout << "Input Name " << inputNames.back() << std::endl;
    }

    if (outputName.empty()) {
        std::cout << "NO OUTPUT FILE!!!!" << std::endl;
    }

    // Attach to the input tree.
    std::unique_ptr<TFile> inputFile(
        new TFile(inputNames.back().c_str(),"old"));
    if (!inputFile->IsOpen()) throw std::runtime_error("Input file not open");

    TGeoManager* geom
        = dynamic_cast<TGeoManager*>(inputFile->Get("CubeReconGeometry"));
    if (!geom) {
        std::cout << "No geometry for file" << std::endl;
    }

    /// Attach to the input tree.
    std::unique_ptr<TChain> inputChain(new TChain("CubeEvents"));
    for (int i = 0; i<inputNames.size(); ++i) {
        inputChain->AddFile(inputNames[i].c_str());
    }
    Cube::Event *inputEvent = NULL;
    inputChain->SetBranchAddress("Event",&inputEvent);

    // Open the output file
    std::unique_ptr<TFile> outputFile;
    if (!outputName.empty()) {
        std::cout << "Open Output File: " << outputName << std::endl;
        outputFile.reset(new TFile(outputName.c_str(),"recreate"));
    }
    TTree *outputTree = new TTree("CubeEvents","Reconstructed Event");
    static Cube::Event *outputEvent = inputEvent;
    outputTree->Branch("Event",&outputEvent);

    // Loop through the events.
    int totalEntries = inputChain->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputChain->GetEntry(entry);
        outputEvent = inputEvent;
        std::cout << "Process event "
                  << entry
                  << "/" << inputEvent->GetRunId()
                  << "/" << inputEvent->GetEventId() << std::endl;
        bool save = AnalyzeEvent(*inputEvent);
        if (save) outputTree->Fill();
    }

    if (outputFile) {
        outputFile->Write();
        outputFile->Close();
    }

    return 0;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
