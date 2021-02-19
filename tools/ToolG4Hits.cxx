#include "ToolG4Hits.hxx"

#include <CubeInfo.hxx>

#include <set>

namespace {
    bool G4HitTimeCompare(Cube::Handle<Cube::G4Hit>& lhs,
                          Cube::Handle<Cube::G4Hit>& rhs) {
        return (lhs->GetStart().T() < rhs->GetStart().T());
    }
}

// General comment: This code is pretty inefficient, because I've put
// "clarity" over efficiency.  For instance, a lot of vectors are passed (or
// returned) by value, and are sorted multiple times.

std::vector<Cube::Handle<Cube::G4Hit>>
Cube::Tool::HitG4Hits(Cube::Event& event, Cube::Handle<Cube::Hit> hit) {
    // Collect which trajectories added something to this hit
    std::vector<Cube::Handle<Cube::G4Hit>> result;
    if (Cube::Info::Is3DST(hit->GetIdentifier())
        && hit->GetConstituentCount() > 0) {
        std::set<Cube::Handle<Cube::G4Hit>> segXY;
        std::set<Cube::Handle<Cube::G4Hit>> segXZ;
        std::set<Cube::Handle<Cube::G4Hit>> segYZ;
        for (int c = 0; c<hit->GetConstituentCount(); ++c) {
            Cube::Handle<Cube::Hit> hh = hit->GetConstituent(c);
            if (!hh) continue;
            int projection
                = Cube::Info::IdentifierProjection(hh->GetIdentifier());
            for (int g = 0; g < hh->GetContributorCount(); ++g) {
                int seg = hh->GetContributor(g);
                Cube::Handle<Cube::G4Hit> g4Hit = event.G4Hits[seg];
                if (!g4Hit) continue;
                switch (projection) {
                case Cube::Info::kXYProj: {
                    segXY.insert(g4Hit);
                    break;
                }
                case Cube::Info::kXZProj: {
                    segXZ.insert(g4Hit);
                    break;
                }
                case Cube::Info::kYZProj: {
                    segYZ.insert(g4Hit);
                    break;
                }
                default:
                    throw std::runtime_error("Inconceivable");
                }
            }
        }
        for (std::set<Cube::Handle<Cube::G4Hit>>::iterator xy = segXY.begin();
             xy != segXY.end(); ++xy) {
            if (segXZ.find(*xy) == segXZ.end()) continue;
            if (segYZ.find(*xy) == segYZ.end()) continue;
            result.push_back(*xy);
        }
    }
    else if (Cube::Info::Is3DST(hit->GetIdentifier())
        && hit->GetConstituentCount() < 1) {
        for (int g = 0; g < hit->GetContributorCount(); ++g) {
            int seg = hit->GetContributor(g);
            Cube::Handle<Cube::G4Hit> g4Hit = event.G4Hits[seg];
            if (!g4Hit) continue;
            result.push_back(g4Hit);
        }
    }
    else if (Cube::Info::IsTPC(hit->GetIdentifier())) {
        if (hit->GetConstituentCount() > 0) {
            for (int c = 0; c<hit->GetConstituentCount(); ++c) {
                Cube::Handle<Cube::Hit> hh = hit->GetConstituent(c);
                std::vector<Cube::Handle<Cube::G4Hit>> tmp(HitG4Hits(event,hh));
                std::copy(tmp.begin(),tmp.end(),std::back_inserter(result));
            }
        }
        else if (hit->GetContributorCount() > 0) {
            for (int g = 0; g < hit->GetContributorCount(); ++g) {
                int seg = hit->GetContributor(g);
                Cube::Handle<Cube::G4Hit> g4Hit = event.G4Hits[seg];
                if (!g4Hit) continue;
                result.push_back(g4Hit);
            }
        }
        else {
            std::cout << "Inconcievable."
                      << " Can't get hit segments for this TPC hit"
                      << std::endl;
            throw std::runtime_error("TPC Hit is invalid");
        }
    }
    else if (Cube::Info::IsECal(hit->GetIdentifier())) {
        std::cout << "Can't get hit segments for ECal yet" << std::endl;
    }
    return result;
}

std::vector<Cube::Handle<Cube::G4Hit>>
Cube::Tool::SelectionG4Hits(Cube::Event& event, Cube::HitSelection& hits) {
    // Collect which trajectories added something to this hit selection
    std::set<Cube::Handle<Cube::G4Hit>> segSet;
    for (Cube::HitSelection::iterator h = hits.begin();
         h != hits.end(); ++h) {
        std::vector<Cube::Handle<Cube::G4Hit>> hitResult
            = HitG4Hits(event,*h);
        segSet.insert(hitResult.begin(), hitResult.end());
    }
    std::vector<Cube::Handle<Cube::G4Hit>> result;
    std::copy(segSet.begin(), segSet.end(), std::back_inserter(result));
    std::sort(result.begin(), result.end());
    std::vector<Cube::Handle<Cube::G4Hit>>::iterator end
        = std::unique(result.begin(),result.end());
    result.erase(end,result.end());
    std::sort(result.begin(), result.end(), G4HitTimeCompare);
    return result;
}

std::vector<Cube::Handle<Cube::G4Hit>>
Cube::Tool::ObjectG4Hits(Cube::Event& event, Cube::ReconObject& object) {
    std::vector<Cube::Handle<Cube::G4Hit>> result;
    Cube::Handle<Cube::HitSelection> hits = object.GetHitSelection();
    if (!hits) return result;
    result = SelectionG4Hits(event,*hits);
    return result;
}

std::vector<Cube::Handle<Cube::G4Hit>>
Cube::Tool::TrajectoryG4Hits(Cube::Event& event, int trackId) {
    std::vector<Cube::Handle<Cube::G4Hit>> result;
    for (Cube::Event::G4HitContainer::iterator g = event.G4Hits.begin();
         g != event.G4Hits.end(); ++g) {
        if (g->second->GetPrimaryId() != trackId) continue;
        result.push_back(g->second);
    }
    std::sort(result.begin(), result.end());
    std::vector<Cube::Handle<Cube::G4Hit>>::iterator end
        = std::unique(result.begin(),result.end());
    result.erase(end,result.end());
    std::sort(result.begin(), result.end(), G4HitTimeCompare);
    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
