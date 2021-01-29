#ifndef ERepSimExport_hxx_seen
#define ERepSimExport_hxx_seen

#include <TFile.h>
#include <TTree.h>
#include <TRef.h>

#include <vector>
#include <map>

namespace ERepSim {
    class Input;
}

/// A singleton setup to read the ERepSim tree that is suitable to be copied
/// into other projects.  The format of the ERepSim tree is part of the
/// ERepSim API, so I'll try to keep it backward compatible.
class ERepSim::Input {
private:
    Input()
        : DataTree(NULL), RunId(0), EventId(0),
          HitSensorId(NULL),
          HitX(NULL), HitY(NULL), HitZ(NULL),
          HitTime(NULL), HitTimeWidth(NULL), HitCharge(NULL),
          HitSegmentBegin(NULL), HitSegmentEnd(NULL),
          SegmentIds(NULL), SegmentTrackId(NULL),
          SegmentPDG(NULL), SegmentEnergy(NULL),
          SegmentX1(NULL), SegmentY1(NULL), SegmentZ1(NULL),
          SegmentX2(NULL), SegmentY2(NULL), SegmentZ2(NULL),
          SegmentT(NULL),
          TrajectoryId(NULL), TrajectoryParent(NULL), TrajectoryPDG(NULL),
          TrajectoryX(NULL), TrajectoryY(NULL), TrajectoryZ(NULL),
          TrajectoryT(NULL),
          TrajectoryPx(NULL), TrajectoryPy(NULL), TrajectoryPz(NULL),
          TrajectoryPe(NULL),
          PropertyTree(NULL), Property(NULL), PropertyString(NULL) {}
    Input(const Input&) {}
    Input& operator=(const Input&) {return *ThisInput;}
    static Input *ThisInput;

public:
    static Input& Get(void)  {
        if (!ThisInput) {
            ThisInput = new ERepSim::Input;
        }
        return *ThisInput;
    }

    TTree* DataTree;
    int RunId;
    int EventId;
    std::vector<int> *HitSensorId;
    std::vector<float> *HitX;
    std::vector<float> *HitY;
    std::vector<float> *HitZ;
    std::vector<float> *HitTime;
    std::vector<float> *HitTimeWidth;
    std::vector<float> *HitCharge;
    std::vector<int> *HitSegmentBegin;
    std::vector<int> *HitSegmentEnd;

    // The identifier of the hit segment.  This is also the index of the hit
    // segment in the edep-sim TG4HitSegment in the TG4HitSegmentContainer.
    std::vector<int> *SegmentIds;

    // The trajectory id of the main particle for this segment (i.e. the first
    // particle that contributed to the segment).
    std::vector<int> *SegmentTrackId;

    // The parameters of the segment.
    std::vector<int>   *SegmentPDG;
    std::vector<float> *SegmentEnergy;
    std::vector<float> *SegmentX1;
    std::vector<float> *SegmentY1;
    std::vector<float> *SegmentZ1;
    std::vector<float> *SegmentX2;
    std::vector<float> *SegmentY2;
    std::vector<float> *SegmentZ2;
    std::vector<float> *SegmentT;

    // The identifier for the trajectory from geant.  The value should be
    // equal to the index.
    std::vector<int> *TrajectoryId;

    // The identifier of the parent particle.  The parent will be negative
    // when this is a primary particle (check, it might be zero).
    std::vector<int> *TrajectoryParent;

    // The particle information.
    std::vector<int> *TrajectoryPDG;
    std::vector<float> *TrajectoryX;
    std::vector<float> *TrajectoryY;
    std::vector<float> *TrajectoryZ;
    std::vector<float> *TrajectoryT;
    std::vector<float> *TrajectoryPx;
    std::vector<float> *TrajectoryPy;
    std::vector<float> *TrajectoryPz;
    std::vector<float> *TrajectoryPe;

    TTree* PropertyTree;
    std::map<std::string,double> *Property;
    std::map<std::string,std::string> *PropertyString;

    void Attach(TFile* inFile) {
        PropertyTree = dynamic_cast<TTree*>(inFile->Get("ERepSimProperties"));
        if (PropertyTree) {
            PropertyTree->SetBranchAddress("Property",&Property);
            PropertyTree->SetBranchAddress("PropertyString",&PropertyString);
            PropertyTree->GetEntry(0);
        }

        DataTree = dynamic_cast<TTree*>(inFile->Get("ERepSimEvents"));
        if (DataTree) {
            DataTree->SetBranchAddress("RunId",&RunId);
            DataTree->SetBranchAddress("EventId",&EventId);
            DataTree->SetBranchAddress("HitSensorId",&HitSensorId);
            DataTree->SetBranchAddress("HitX",&HitX);
            DataTree->SetBranchAddress("HitY",&HitY);
            DataTree->SetBranchAddress("HitZ",&HitZ);
            DataTree->SetBranchAddress("HitTime",&HitTime);
            DataTree->SetBranchAddress("HitTimeWidth",&HitTimeWidth);
            DataTree->SetBranchAddress("HitCharge",&HitCharge);
            DataTree->SetBranchAddress("HitSegmentBegin",&HitSegmentBegin);
            DataTree->SetBranchAddress("HitSegmentEnd",&HitSegmentEnd);

            DataTree->SetBranchAddress("SegmentIds",&SegmentIds);
            DataTree->SetBranchAddress("SegmentTrackId",&SegmentTrackId);
            DataTree->SetBranchAddress("SegmentEnergy",&SegmentEnergy);
            DataTree->SetBranchAddress("SegmentPDG",&SegmentPDG);
            DataTree->SetBranchAddress("SegmentX1",&SegmentX1);
            DataTree->SetBranchAddress("SegmentY1",&SegmentY1);
            DataTree->SetBranchAddress("SegmentZ1",&SegmentZ1);
            DataTree->SetBranchAddress("SegmentX2",&SegmentX2);
            DataTree->SetBranchAddress("SegmentY2",&SegmentY2);
            DataTree->SetBranchAddress("SegmentZ2",&SegmentZ2);
            DataTree->SetBranchAddress("SegmentT",&SegmentT);

            DataTree->SetBranchAddress("TrajectoryId",&TrajectoryId);
            DataTree->SetBranchAddress("TrajectoryParent",&TrajectoryParent);
            DataTree->SetBranchAddress("TrajectoryPDG",&TrajectoryPDG);
            DataTree->SetBranchAddress("TrajectoryX",&TrajectoryX);
            DataTree->SetBranchAddress("TrajectoryY",&TrajectoryY);
            DataTree->SetBranchAddress("TrajectoryZ",&TrajectoryZ);
            DataTree->SetBranchAddress("TrajectoryT",&TrajectoryT);
            DataTree->SetBranchAddress("TrajectoryPx",&TrajectoryPx);
            DataTree->SetBranchAddress("TrajectoryPy",&TrajectoryPy);
            DataTree->SetBranchAddress("TrajectoryPz",&TrajectoryPz);
            DataTree->SetBranchAddress("TrajectoryPe",&TrajectoryPe);
        }
    }

    void GetEntry(int i) {
        DataTree->GetEntry(i);
        PropertyTree->GetEntry(0);
    }
};

ERepSim::Input* ERepSim::Input::ThisInput = NULL;
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/erep-build.sh force"
// End:
