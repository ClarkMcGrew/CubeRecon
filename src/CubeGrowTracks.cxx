#include "CubeGrowTracks.hxx"
#include "CubeTrackFit.hxx"
#include "CubeCreateTrack.hxx"
#include "CubeCompareReconObjects.hxx"
#include "CubeMakeUsed.hxx"

#include <CubeHandle.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>

#include <CubeLog.hxx>
#include <CubeUnits.hxx>

#include <TMatrixD.h>

#include <memory>
#include <list>
#include <cmath>

Cube::GrowTracks::GrowTracks()
    : Cube::Algorithm("GrowTracks",
                 "Merge tracks that are end-to-end") {
    fMergeDistanceCut = 50.0*unit::mm;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.MergeDistance");

    fMinimumFollowDistanceCut = -10.0*unit::mm;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.MinimumDistance");

    fFollowCosineCut = 30.0*unit::deg;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.AllowedKink");
    fFollowCosineCut = std::cos(fFollowCosineCut);

    fGoodnessCut = 16.0;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.GoodnessCut");

    fMatchedCosineCut = 15.0*unit::deg;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.MatchedAngleCut");
    fMatchedCosineCut = std::cos(fMatchedCosineCut);

    fMatchedPositionCut = 15.0*unit::mm;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    //     "sfgRecon.GrowTracks.MatchedPositionCut");
}

Cube::GrowTracks::~GrowTracks() { }

double
Cube::GrowTracks::DistToApproach(const TVector3& dest,
                                   const TVector3& src,
                                   const TVector3& dir) {
    return (dest-src)*dir;
}

double
Cube::GrowTracks::ApproachDist(const TVector3& dest,
                                 const TVector3& src,
                                 const TVector3& dir) {
    double d = DistToApproach(dest,src,dir);
    return (dest - (src + d*dir)).Mag();
}

Cube::GrowTracks::Orientation
Cube::GrowTracks::TrackOrientation(Cube::Handle<Cube::ReconTrack> first,
                                     Cube::Handle<Cube::ReconTrack> second) {
    /// The positions of the front and back of the tracks.
    Cube::Handle<Cube::TrackState> front1
        = first->GetNodes().front()->GetState();
    Cube::Handle<Cube::TrackState> back1
        = first->GetNodes().back()->GetState();
    Cube::Handle<Cube::TrackState> front2
        = second->GetNodes().front()->GetState();
    Cube::Handle<Cube::TrackState> back2
        = second->GetNodes().back()->GetState();

    // Now figure out distances between the track ends.  The two closest
    // distance determine how the tracks will be arranged.
    double ffDist = (front1->GetPosition().Vect()
                     - front2->GetPosition().Vect()).Mag();
    double fbDist = (front1->GetPosition().Vect()
                     - back2->GetPosition().Vect()).Mag();
    double bfDist = (back1->GetPosition().Vect()
                     - front2->GetPosition().Vect()).Mag();
    double bbDist = (back1->GetPosition().Vect()
                     - back2->GetPosition().Vect()).Mag();

    // Check that the tracks are at least close!
    double minDist = std::min(ffDist,
                              std::min(fbDist,
                                       std::min(bfDist, bbDist)));
    if (minDist > fMergeDistanceCut) {
        return kNotClose;
    }

    // Figure out which ends are closest and use that to determine which order
    // the clusters get added.
    if (ffDist<fbDist && ffDist<bfDist && ffDist<bbDist) {
        // A front to front track.
        return kFrontFront;
    }
    else if (fbDist<ffDist && fbDist<bfDist && fbDist<bbDist) {
        // A front to back track.
        return kFrontBack;
    }
    else if (bfDist<fbDist && bfDist<ffDist && bfDist<bbDist) {
        // A back to front track.
        return kBackFront;
    }
    else if (bbDist<fbDist && bbDist<bfDist && bbDist<ffDist) {
        // A back to back track.
        return kBackBack;
    }

    // This can't happen!
    CUBE_ERROR << "This can't happen!" << std::endl;
    return kNotClose;
}

// t1, t2, pos, dir, cov, pos, dir, cov)

Cube::GrowTracks::Orientation
Cube::GrowTracks::MatchedStates(Cube::Handle<Cube::ReconTrack> t1,
                                  Cube::Handle<Cube::ReconTrack> t2,
                                  Cube::Handle<Cube::TrackState>& leading,
                                  double& reverseLeading,
                                  Cube::Handle<Cube::TrackState>& following,
                                  double& reverseFollowing) {
    leading = Cube::Handle<Cube::TrackState>();
    following = Cube::Handle<Cube::TrackState>();
    reverseLeading = 1;
    reverseFollowing = 1;

    Orientation orient = TrackOrientation(t1,t2);
    if (orient == kNotClose) return kNotClose;

    // For the purposes of the goodness calculation, the tracks are arranged
    // so that track1 is always "leading" and track2 is always following.  The
    // directions will be reversed so that they point from "leading" to
    // "following".
    switch (orient) {
    case kBackFront: // t1 leads
        leading = t1->GetNodes().back()->GetState();
        reverseLeading = 1;
        following = t2->GetNodes().front()->GetState();
        reverseFollowing = 1;
        break;
    case kBackBack: // t1 leads
        leading = t1->GetNodes().back()->GetState();
        reverseLeading = 1;
        following = t2->GetNodes().back()->GetState();
        reverseFollowing = -1;
        break;
    case kFrontFront: // t1 leads
        leading = t1->GetNodes().front()->GetState();
        reverseLeading = -1;
        following = t2->GetNodes().front()->GetState();
        reverseFollowing = 1;
        break;
    case kFrontBack: // t2 leads (don't reverse both)
        leading = t2->GetNodes().back()->GetState();
        reverseLeading = 1;
        following = t1->GetNodes().front()->GetState();
        reverseFollowing = 1;
        break;
    default:
        CUBE_ERROR << "This can't happen!" << std::endl;
        throw std::runtime_error("Invalid track orientation.");
    }

    return orient;
}

double
Cube::GrowTracks::MatchGoodness(Cube::Handle<Cube::ReconTrack> t1,
                                  Cube::Handle<Cube::ReconTrack> t2) {
    // The state for the track that will lead the pair of tracks.
    Cube::Handle<Cube::TrackState> leading;
    double reverseLeading = false;

    // The state for the track that will follow the heading track.
    Cube::Handle<Cube::TrackState> following;
    double reverseFollowing = false;

    Orientation orient = MatchedStates(t1, t2,
                                       leading, reverseLeading,
                                       following, reverseFollowing);

    TVector3 leadingPos = leading->GetPosition().Vect();
    TVector3 leadingDir = reverseLeading*leading->GetDirection();

    TVector3 followingPos = following->GetPosition().Vect();
    TVector3 followingDir = reverseFollowing*following->GetDirection();

    // Find the overall direction error matrix.
    TMatrixD dirErr(3,3);
    for (int i=0; i<3; ++i) {
        for (int j = 0; j<3; ++j) {
            dirErr(i,j) = leading->GetDirectionCovariance(i,j);
            dirErr(i,j) += following->GetDirectionCovariance(i,j);
        }
    }
    dirErr.InvertFast();

    // Find out how well the directions match.
    TVector3 dirDiff = followingDir - leadingDir;
    double dirGoodness = dirDiff*(dirErr*dirDiff);

    // Find the overall position error matrix.  This doesn't add in the extra
    // uncertainty from the direction since the two positions are suppose to
    // be close together so the additional variance is suppose to be small.
    // This means that the final chi2 will be microscopically larger than it
    // should be,
    TMatrixD posErr(3,3);
    for (int i=0; i<3; ++i) {
        for (int j = 0; j<3; ++j) {
            posErr(i,j) = leading->GetPositionCovariance(i,j);
            posErr(i,j) += following->GetPositionCovariance(i,j);
        }
    }
    posErr.InvertFast();

    double followingDist = DistToApproach(
        following->GetPosition().Vect(),
        leading->GetPosition().Vect(),
        reverseLeading*leading->GetDirection());

    TVector3 posDiff = followingPos - (leadingPos + followingDist*leadingDir);
    double posGoodness = posDiff*(posErr*posDiff);

    double result = dirGoodness + posGoodness;

    return result;
}


Cube::Handle<Cube::AlgorithmResult>
Cube::GrowTracks::Process(const Cube::AlgorithmResult& input,
                          const Cube::AlgorithmResult&,
                          const Cube::AlgorithmResult&) {
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input.GetObjectContainer();

    CUBE_LOG(0) << "Process Cube::GrowTracks" << std::endl;

    if (!inputObjects || !inputHits) {
        CUBE_ERROR << "No input objects or hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    // Create the output containers.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // Create a stack to keep tracks in.  When the stack is empty, all the
    // tracks that need to be merge have been merge.
    typedef std::list< Cube::Handle<Cube::ReconTrack> > TrackList;
    TrackList trackList;

    // Make a copy of all of the tracks in the input objects, and save the
    // non-tracks to the output.
    for (Cube::ReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        Cube::Handle<Cube::ReconTrack> track = *t;
        if (!track) {
            finalObjects->push_back(*t);
            continue;
        }
        trackList.push_back(track);
    }

    // Sort the tracks so that the longest ones are first.
    trackList.sort(Cube::CompareReconObjects());

    // Pop a track off the stack and see if it should be merged.  If the track
    // is merged, the result is pushed back on the stack.  If it doesn't get
    // merge, the track get's pushed into the final object container.
    while (!trackList.empty()) {
        Cube::Handle<Cube::ReconTrack> track1 = trackList.front();
        trackList.pop_front();  // This removes the track from the list!

        // Don't use a very short track as a base for merging.  The tracks
        // must be in order of number of nodes (decreasing), and a short track
        // will already be merged, or doesn't have a good partner.
        if (track1->GetNodes().size() < 3) {
            finalObjects->push_back(track1);
            continue;
        }

        int index = 0;
        for (TrackList::iterator t = trackList.begin();
             t!=trackList.end(); ++t) {
            ++index;
            Cube::Handle<Cube::ReconTrack> track2 = *t;

            bool goodMatch = false;
            do {
                // The state for the track that will lead the pair of tracks.
                Cube::Handle<Cube::TrackState> leading;
                double reverseLeading = false;

                // The state for the track that will follow the heading track.
                Cube::Handle<Cube::TrackState> following;
                double reverseFollowing = false;

                // Check that the orientation is valid!
                Orientation orient = MatchedStates(track1, track2,
                                                   leading, reverseLeading,
                                                   following, reverseFollowing);
                if (orient == kNotClose) {
                    break;
                }

                double followingDist = DistToApproach(
                    following->GetPosition().Vect(),
                    leading->GetPosition().Vect(),
                    reverseLeading*leading->GetDirection());

                // Make sure the tracks don't overlap (i.e. following actually
                // follows leading).
                if (followingDist < fMinimumFollowDistanceCut) {
                    goodMatch = false;
                    break;
                }

                // Make sure the directions are pointing in the same general
                // direction.  Opposing directions can happen if the tracks
                // connect at a "V" kink.
                double followingCos
                    = leading->GetDirection()*following->GetDirection();
                followingCos *= reverseLeading*reverseFollowing;
                if (followingCos < fFollowCosineCut) {
                    goodMatch = false;
                    break;
                }

                // If the two tracks have a good chi-squared, then always
                // merge.
                double match = MatchGoodness(track1,track2);
                if (match < fGoodnessCut) {
                    goodMatch = true;
                    break;
                }

                // If the directions don't match, then this is not a good
                // match.
                if (followingCos < fMatchedCosineCut) {
                    goodMatch = false;
                    break;
                }

                // If the projected end points are close (and the directions
                // match, already checked), then this is a good match
                TVector3 posDiff = following->GetPosition().Vect()
                    - (leading->GetPosition().Vect()
                       + reverseLeading*followingDist*leading->GetDirection());
                if (posDiff.Mag() < fMatchedPositionCut) {
                    goodMatch = true;
                    break;
                }

                goodMatch = false;
            } while (false);
            if (!goodMatch) continue;

            ///////////////////////////////////////////////////////
            // If we get here then the tracks should be merged.
            ///////////////////////////////////////////////////////

            // Remove the track iterator (the track is held in track2).
            trackList.erase(t);

            // Merge the tracks.  This preserves the direction of track1.  The
            // direction of track2 may be reversed.
            Orientation orient = TrackOrientation(track1,track2);
            Cube::Handle<Cube::ReconTrack> merged;
            if (orient == kBackFront || orient == kBackBack) {
                merged = Cube::CombineTracks(*track1,*track2);
            }
            else {
                merged = Cube::CombineTracks(*track2,*track1);
            }

            // Put the new track back into the front of the list (track1 was
            // the longest, so the merged track should also be the longest.
            trackList.push_front(merged);

            // Clear out the track variable.  This releases the handle.
            track1 = Cube::Handle<Cube::ReconTrack>();

            // Don't continue the inner loop.
            break;
        }

        // If we get to the bottom of the loop with a track, then push it on
        // to the final objects.
        if (track1) {
            finalObjects->push_back(track1);
        }

    }

    std::sort(finalObjects->begin(), finalObjects->end(),
              Cube::CompareReconObjects());

    Cube::TrackFit fitter;
    for (Cube::ReconObjectContainer::iterator o = finalObjects->begin();
         o != finalObjects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) continue;
        if (track->CheckStatus(Cube::ReconObject::kSuccess)) {
            continue;
        }
        *o = fitter(track);
    }

    result->AddObjectContainer(finalObjects);

    // Build the hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
