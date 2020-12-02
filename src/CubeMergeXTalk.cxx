#include "CubeMergeXTalk.hxx"
#include "CubeTrackFit.hxx"
#include "CubeCreateTrack.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeInfo.hxx"

#include <CubeHandle.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeLog.hxx>
#include <TTmplDensityCluster.hxx>

#include <TVector3.h>

#include <memory>
#include <vector>
#include <map>
#include <set>
#include <cmath>

Cube::MergeXTalk::MergeXTalk()
    : Cube::Algorithm("MergeXTalk",
                      "Merge candidate crosstalk with the tracks") {}

Cube::MergeXTalk::~MergeXTalk() {}

namespace {
    // Determine the proximity between cubes.  The proximity is defined as the
    // number of "cube to cube" jumps needed to get between the cubes.
    struct CubeProximity {
        double operator()(Cube::Handle<Cube::Hit> lhs,
                          Cube::Handle<Cube::Hit> rhs) {
            TVector3 delta = lhs->GetPosition() - rhs->GetPosition();
            double d = std::max(std::abs(delta.X()), std::abs(delta.Y()));
            return std::max(d,std::abs(delta.Z()));
        }
    };
}


Cube::Handle<Cube::AlgorithmResult>
Cube::MergeXTalk::Process(const Cube::AlgorithmResult& input,
                          const Cube::AlgorithmResult&,
                          const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Process Cube::MergeXTalk" << std::endl;
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input.GetObjectContainer();

    if (!inputObjects || !inputHits) {
        CUBE_ERROR << "No input objects or hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    // Create the output containers.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // A list of clusters that should be considered.
    typedef std::list<Cube::Handle<Cube::ReconCluster>> ClusterList;
    ClusterList clusterList;
    ClusterList leftOverClusters;

    // Make a copy of all of the tracks and clusters into the work lists, and
    // save the non-tracks to the output.
    AllTracks allTracks;
    AllHits allHits;
    for (Cube::ReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        Cube::Handle<Cube::ReconTrack> track = *t;
        if (track) {
            TrackNodeHits trackNodeHits;
            for (Cube::ReconNodeContainer::iterator n
                     = track->GetNodes().begin();
                 n != track->GetNodes().end(); ++n) {
                Cube::Handle<Cube::HitSelection> hits
                    = (*n)->GetObject()->GetHitSelection();
                NodeHits c(new HitSet);
                // Insert the hits into a set for each node.
                c->insert(hits->begin(), hits->end());
                // Add the vector of the set of hits in each node to the track
                trackNodeHits.push_back(c);
                // Make a map between the hits to the sets they are in.
                for (HitSet::iterator h = c->begin(); h != c->end(); ++h) {
                    allHits[(*h)].insert(c);
                }
            }
            allTracks.push_back(trackNodeHits);
            continue;
        }
        Cube::Handle<Cube::ReconCluster> cluster = *t;
        if (cluster) {
            clusterList.push_back(cluster);
            continue;
        }
        finalObjects->push_back(*t);
    }

    CUBE_LOG(0) << "xtalk :: Merge the clusters" << std::endl;
    // Check each cluster to see if it is consistent with a cluster made of
    // cross talk hits.
    for (ClusterList::iterator c = clusterList.begin();
         c != clusterList.end(); ++c) {
        int neighbors = CountClusterNeighbors(allHits, *c);
        int nHits = (*c)->GetHitSelection()->size();
        // If a cluster has hits that are not neighboring a track, then assume
        // it is not made up of cross talk.  Notice that a cluster with cross
        // talk will contain one hit in the track, and one or more hits
        // neighboring the track.
        if (nHits > neighbors) {
            leftOverClusters.push_back(*c);
            continue;
        }
        // All of the hits are neighbors to the track, so assume all the hits
        // are crosstalk, or members of a track.  Try to add them to all of
        // the tracks that contain the best neighbor cube.
        Cube::Handle<Cube::HitSelection> hits = (*c)->GetHitSelection();
        for (Cube::HitSelection::iterator h = hits->begin();
             h != hits->end(); ++h) {
            Cube::Handle<Cube::Hit> bestNeighbor = FindBestNeighbor(allHits,*h);
            if (!bestNeighbor) {
                CUBE_ERROR << "No neighbor, but there should be one."
                           << std::endl;
                continue;
            }
            // Add the current hit to all of the tracks that include
            // the bestNeighbor.
            for (std::set<NodeHits>::iterator s
                     = allHits[bestNeighbor].begin();
                 s != allHits[bestNeighbor].end(); ++s) {
                (*s)->insert(*h);
            }
        }
    }

    CUBE_LOG(0) << "xtalk: Rebuild the tracks " << allTracks.size()
                << std::endl;

    // Build the new tracks from the allTracks vector and add them to
    // finalObjects.  Keep track of all of the hits in the tracks.
    HitSet trackHits;
    Cube::TrackFit fitter;
    for (AllTracks::iterator t = allTracks.begin();
         t != allTracks.end(); ++t) {
        // Build the clusters for the nodes.
        std::vector<Cube::Handle<Cube::ReconCluster>> clusters;
        for (TrackNodeHits::iterator n = t->begin(); n != t->end(); ++n) {
            Cube::Handle<Cube::ReconCluster> cluster
                = Cube::CreateCluster("mergeXTalk",
                                     (*n)->begin(), (*n)->end());
            clusters.push_back(cluster);
            // Add to the set of all track hits.
            trackHits.insert((*n)->begin(),(*n)->end());
        }
        // Create the track
        Cube::Handle<Cube::ReconTrack> track
            = Cube::CreateTrackFromClusters(
                "xtalk",clusters.begin(), clusters.end());
        // Fit the track
        track = fitter(track);
        // Add to finalObjects
        finalObjects->push_back(track);
    }

#ifdef RECLUSTER_LEFTOVERS
    // build new clusters from the hits in leftOverClusters.
    HitSet clusterHits;
    for (ClusterList::iterator c = leftOverClusters.begin();
         c != leftOverClusters.end(); ++c) {
        for (Cube::HitSelection::iterator h = (*c)->GetHitSelection()->begin();
             h != (*c)->GetHitSelection()->end(); ++h) {
            if (trackHits.find(*h) != trackHits.end()) continue;
            clusterHits.insert(*h);
        }
    }

    // Apply the DBScan algorithm.
    // Work around template parsing bug in some GCC versions...
    typedef Cube::Handle<Cube::Hit> Arg;
    typedef TTmplDensityCluster<Arg, CubeProximity> ClusterCubes;
    ClusterCubes clusterCubes(1, 28.0); // Just less than 3 cubes.
    clusterCubes.Cluster(clusterHits.begin(), clusterHits.end());
    for (int i=0; i<clusterCubes.GetClusterCount(); ++i) {
        const ClusterCubes::Points& points = clusterCubes.GetCluster(i);
        Cube::Handle<Cube::ReconCluster> cluster
            = Cube::CreateCluster("mergeXTalk",
                                  points.begin(),points.end());
        finalObjects->push_back(cluster);
    }
#endif

    result->AddObjectContainer(finalObjects);

    // Build the hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

// Count the number of cubes in allHits that are neighboring to hit.
int Cube::MergeXTalk::CountHitNeighbors(AllHits& allHits,
                                        Cube::Handle<Cube::Hit>& hit) {
    CubeProximity proximity;
    int neighborHits = 0;
    for (AllHits::iterator h = allHits.begin(); h != allHits.end(); ++h) {
        Cube::Handle<Cube::Hit> trackHit = h->first;
        double diff = proximity(trackHit,hit);
        if (std::abs(diff) < 18.0) ++neighborHits;
    }
    return neighborHits;
}

// Count the number of cubes in a cluster that have a neighbor.
int Cube::MergeXTalk::CountClusterNeighbors(
    AllHits& allHits,
    Cube::Handle<Cube::ReconCluster>& cluster) {
    Cube::Handle<Cube::HitSelection> hits = cluster->GetHitSelection();
    int clusterNeighbors = 0;
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        int hitNeighbors = CountHitNeighbors(allHits,*h);
        if (hitNeighbors>0) ++clusterNeighbors;
    }
    return clusterNeighbors;
}

Cube::Handle<Cube::Hit> Cube::MergeXTalk::FindBestNeighbor(
    AllHits& allHits, Cube::Handle<Cube::Hit>& hit) {
    CubeProximity proximity;
    Cube::Handle<Cube::Hit> bestHit;
    for (AllHits::iterator h = allHits.begin(); h != allHits.end(); ++h) {
        Cube::Handle<Cube::Hit> trackHit = h->first;
        double dist = proximity(trackHit,hit);
        // If hit and trackHit are the same, it's the best hit.
        if (dist < 1.0) return trackHit;
        // The trackHit and hit are next to each other.  Take the biggest
        // charge.
        if (dist > 18.0) continue;
        if (!bestHit) bestHit = trackHit;
        else if (bestHit->GetCharge() < trackHit->GetCharge()) {
            bestHit = trackHit;
        }
    }
    return bestHit;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
