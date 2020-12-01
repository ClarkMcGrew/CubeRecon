#include "CubeSpanningTree.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeReconCluster.hxx>
#include <CubeAlgorithmResult.hxx>
#include <TTmplMinimalSpanningTree.hxx>
#include <set>

namespace {
    ///////////////////////////////////////////////////////////////////////
    // An evil way to control the way the edge weight is calculated.
    int EdgeWeight_DistanceType = 0;

    ///////////////////////////////////////////////////////////////////////
    // Base the weight between hits on the distance.  Apply a small correction
    // to break the degenerancy between hits that are the same distance apart.
    // The "distance" is raise by a small amount based on the charge of the
    // hits.  If the hits have a small charge it's raised by more than if they
    // have a large charge.  This favors connecting high charge hits first.
    template<typename HitHandle>
    class CubeEdgeWeight {
    public:
        double operator()(HitHandle lhs,HitHandle rhs){
            TVector3 diff = lhs->GetPosition() - rhs->GetPosition();
            double dist;
            switch (EdgeWeight_DistanceType)  {
            case 1:
                dist = diff.Mag();
                break;
            default:
                dist = std::max(std::abs(diff.X()),std::abs(diff.Y()));
                dist = std::max(dist, std::abs(diff.Z()));
                break;
            }
            double chargeSum = lhs->GetCharge() + rhs->GetCharge();
            const double epsilon = 0.01;
            double chargeCorr = epsilon/(1.0+chargeSum);
            return dist + chargeCorr;
        }
    };

    // A user data struct for keeping track of the number of children.
    struct CubeUserData {
        int Visited; // Not zero if visited.
        int Children; // The number of children below this vertex.
    };
}

Cube::SpanningTree::SpanningTree()
    : Cube::Algorithm("SpanningTree") {

    fDistanceType = 0;
    fOversizeCut = 3000;

}

Cube::Handle<Cube::AlgorithmResult>
Cube::SpanningTree::Process(const Cube::AlgorithmResult& input,
                            const Cube::AlgorithmResult&,
                            const Cube::AlgorithmResult&) {
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    CUBE_LOG(0) << "Process Cube::SpanningTree" << std::endl;

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        CUBE_ERROR << "No hits" << std::endl;
        return result;
    }
    if (inputHits->empty()) {
        CUBE_ERROR << "Hits empty" << std::endl;
        return result;
    }

    // Create the container for the unprocessed objects.
    Cube::Handle<Cube::ReconObjectContainer>
        unprocessedObjects(new Cube::ReconObjectContainer("unprocessed"));
    result->AddObjectContainer(unprocessedObjects);

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    // Work around template parsing bug in some GCC versions...
    typedef Cube::Handle<Cube::Hit> HitHandle;

    // Make a typedef for the Minimal Spanning Tree
    EdgeWeight_DistanceType = fDistanceType;
    typedef TTmplMinimalSpanningTree<
        HitHandle,
        CubeEdgeWeight<HitHandle>,
        CubeUserData> CubeTree;

    // Apply the minimal spanning tree to each input object.
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input.GetObjectContainer("final");
    for (Cube::ReconObjectContainer::iterator o = inputObjects->begin();
         o != inputObjects->end(); ++ o) {
        Cube::Handle<Cube::HitSelection> objectHits = (*o)->GetHitSelection();
        if (!objectHits) continue;

        if (objectHits->empty()) {
            continue;
        }

        // Make a set of hits to be completely sure there are no duplicates.
        std::set<Cube::Handle<Cube::Hit>> hitSet;
        for (Cube::HitSelection::iterator h = objectHits->begin();
             h != objectHits->end(); ++h) {
            hitSet.insert((*h));
        }

        if (hitSet.empty()) continue;

        if (hitSet.size() < 4) {
            Cube::Handle<Cube::ReconCluster> cluster
                = Cube::CreateCluster("spanningTreeSmall",
                                      hitSet.begin(), hitSet.end());
            unprocessedObjects->push_back(cluster);
            continue;
        }

        if (hitSet.size() > fOversizeCut) {
            Cube::Handle<Cube::ReconCluster> cluster
                = Cube::CreateCluster("spanningTreeLarge",
                                      hitSet.begin(), hitSet.end());
            CUBE_LOG(0) << "SpanningTree:: Reject large cluster w/ "
                        << hitSet.size() << " hits"
                        << std::endl;
            unprocessedObjects->push_back(cluster);
            continue;
        }

        CUBE_LOG(0) << "SpanningTree:: Build tree w/ "
                    << hitSet.size() << " hits"
                    << std::endl;

        // Build an MST and use it to find the deepest leaf.  The goal is to
        // find a hit that is at one end of a track.
        std::unique_ptr<CubeTree> makeTree1(new CubeTree);
        makeTree1->AddVertices(hitSet.begin(), hitSet.end());
        makeTree1->MakeTree(*(hitSet.begin()));

        const CubeTree::MinimalSpanningTree& tree1 = makeTree1->GetTree();

        int maxDepth = -1;
        int deepestVertex = -1;
        for (int i = 0; i < tree1.size(); ++i) {
            if (maxDepth < tree1[i].VertexDepth) {
                maxDepth = tree1[i].VertexDepth;
                deepestVertex = i;
            }
        }

        Cube::Handle<Cube::Hit> startingHit = tree1[deepestVertex].Object;

        // Build an MST from the old "deepestVertex" This makes sure that the
        // root is at one end of a track.
        std::unique_ptr<CubeTree> makeTree2(new CubeTree);
        makeTree2->AddVertices(hitSet.begin(), hitSet.end());
        makeTree2->MakeTree(startingHit);

        const CubeTree::MinimalSpanningTree& tree2 = makeTree2->GetTree();

        // Reset the user data.  It's going to be used to tell if a vertex has
        // been visited, and to count the number of children below
        for(int i = 0; i < tree2.size(); ++i) {
            tree2[i].UserData.Visited = 0;
            tree2[i].UserData.Children = 1;
        }

        // Count the number of children.
        while (true) {
            // Find the deepest vertex that has not been visited yet.
            maxDepth = -1;
            deepestVertex = -1;
            for(int i = 0; i < tree2.size(); ++i) {
                if (maxDepth < tree2[i].VertexDepth
                    && !tree2[i].UserData.Visited) {
                    maxDepth = tree2[i].VertexDepth;
                    deepestVertex = i;
                }
            }

            if (deepestVertex < 0) break;

            // Mark the deepestVertex as visited.
            tree2[deepestVertex].UserData.Visited = 1;

            // Add a count of the children to all of the parents.
            for (int vtx = tree2[deepestVertex].Parent;
                 vtx >= 0; vtx = tree2[vtx].Parent) {
                ++tree2[vtx].UserData.Children;
            }
        }

        // Reset the UserData.Visited flag.
        for(int i = 0; i < tree2.size(); ++i) {
            tree2[i].UserData.Visited = 0;
        }

        while (true) {
            // Find the deepest vertex that has not been visited yet.
            maxDepth = -1;
            deepestVertex = -1;
            for(int i = 0; i < tree2.size(); ++i) {
                if (maxDepth < tree2[i].VertexDepth
                    && !tree2[i].UserData.Visited) {
                    maxDepth = tree2[i].VertexDepth;
                    deepestVertex = i;
                }
            }

            // Stop if the next one has too many children, or if it is already
            // visited.


            // All of the vertices have been visited.
            if (deepestVertex < 0) break;

            Cube::HitSelection newHits;

            // Always add the deepest and mark it as visited.
            Cube::Handle<Cube::Hit> hit = tree2[deepestVertex].Object;
            newHits.push_back(hit);

            // Mark the deepest vertex as visited before looking for more hits.
            tree2[deepestVertex].UserData.Visited = 1;

            // Walk toward the root until there is a child.  Add hits for a
            // cluster as we go.  The first hit is the deepest one.  The last
            // hit is either the root of the tree, or at a vertex that has
            // children.
            int lastVtxSeen = deepestVertex;
            int vtx = tree2[deepestVertex].Parent;
            while (vtx >= 0) {

                // Add the hit to a hit selection to build a cluster...
                hit = tree2[vtx].Object;
                newHits.push_back(hit);

                // Check if the vertex has multiple children, and stop if it
                // is a clear branch.
                if (tree2[vtx].Children.size()>1) {
                    int maxChildren = 0;
                    for (int j = 0; j < tree2[vtx].Children.size(); ++j) {
                        int cvtx = tree2[vtx].Children[j];
                        if (cvtx == lastVtxSeen) continue;
                        maxChildren = std::max(maxChildren,
                                               tree2[cvtx].UserData.Children);
                    }
                    if (maxChildren > 2) break;
                }

                // If the vertex that was just added has already been visited,
                // then stop.
                if (tree2[vtx].UserData.Visited) break;

                // Mark the vertex as visited.
                tree2[vtx].UserData.Visited = 1;
                lastVtxSeen = vtx;

                // Move to the next vertex.
                vtx = tree2[vtx].Parent;
            }

            if (newHits.size() < 1) continue;

            // Build a cluster from the hits.
            Cube::Handle<Cube::ReconCluster> cluster(new Cube::ReconCluster);
            cluster->FillFromHits("cluster",newHits.begin(), newHits.end());
            finalObjects->push_back(cluster);
        }
    }

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
