#include <ToolMainTrajectory.hxx>
#include <ToolG4Hits.hxx>

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

#include <map>
#include <vector>

// Find the trajectory that contributed most to the track.  The longest
// trajectory wins.
int Cube::Tool::MainTrajectory(Cube::Event& event,
                               Cube::ReconObject& object) {
    std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
        = Cube::Tool::ObjectG4Hits(event,object);
    std::map<int,double> trajMap;
    for(std::vector<Cube::Handle<Cube::G4Hit>>::iterator
            g = g4Hits.begin(); g != g4Hits.end(); ++g) {
        trajMap[(*g)->GetPrimaryId()]
            += ((*g)->GetStart().Vect() - (*g)->GetStop().Vect()).Mag();
    }
    int maxTraj = -1;
    double maxLen = 0.0;
    for(std::map<int,double>::iterator t = trajMap.begin();
        t != trajMap.end(); ++t) {
        if (maxLen > t->second) continue;
        maxTraj = t->first;
        maxLen = t->second;
    }
    return maxTraj;
}

// Find the trajectory that contributed most to the track.  The longest
// trajectory wins.
double Cube::Tool::MainPurity(Cube::Event& event,
                              Cube::ReconObject& object) {
    std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
        = Cube::Tool::ObjectG4Hits(event,object);
    std::map<int,double> trajMap;
    for(std::vector<Cube::Handle<Cube::G4Hit>>::iterator
            g = g4Hits.begin(); g != g4Hits.end(); ++g) {
        trajMap[(*g)->GetPrimaryId()]
            += ((*g)->GetStart().Vect() - (*g)->GetStop().Vect()).Mag();
    }
    int maxTraj = -1;
    double maxLen = 0.0;
    double totalLen = 0.0;
    for(std::map<int,double>::iterator t = trajMap.begin();
        t != trajMap.end(); ++t) {
        totalLen += t->second;
        if (maxLen > t->second) continue;
        maxTraj = t->first;
        maxLen = t->second;
    }
    if (maxTraj < 0) return 0.0;
    return maxLen/totalLen;
}

// Find the trajectory that contributed most to the track.  The longest
// trajectory wins.
double Cube::Tool::MainCompleteness(Cube::Event& event,
                                    Cube::ReconObject& object) {
    int trackId = MainTrajectory(event,object);
    if (trackId < 0) return 0.0;
    std::vector<Cube::Handle<Cube::G4Hit>> trajSegs
        = Cube::Tool::TrajectoryG4Hits(event,trackId);
    double totalLen = 0.0;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator g = trajSegs.begin();
         g != trajSegs.end(); ++g) {
        totalLen += ((*g)->GetStart().Vect() - (*g)->GetStop().Vect()).Mag();
    }
    std::vector<Cube::Handle<Cube::G4Hit>> objSegs
        = Cube::Tool::ObjectG4Hits(event,object);
    double objectLen = 0.0;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator g = objSegs.begin();
         g != objSegs.end(); ++g) {
        if ((*g)->GetPrimaryId() != trackId) continue;
        objectLen += ((*g)->GetStart().Vect() - (*g)->GetStop().Vect()).Mag();
    }
    if (totalLen < 0.01) return 0.0;
    return objectLen / totalLen;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
