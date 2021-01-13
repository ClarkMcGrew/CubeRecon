#include "CubeInfo.hxx"

#include <CubeLog.hxx>

#include <iostream>

int Cube::Info::SubDetector(int id) {
    int id0 = (id>>27)&0x1F;
    return id0;
}

bool Cube::Info::Is3DST(int id) {
    return (Cube::Info::SubDetector(id) == 13);
}

int Cube::Info::CubeNumber(int id) {
    if (!Cube::Info::Is3DST(id)) {
        throw std::runtime_error("Not the 3DST");
    }
    int idc = (id>>18) & 0x000001FF;
    if (idc > 500) return -1;
    return idc;
}

int Cube::Info::CubeBar(int id) {
    if (!Cube::Info::Is3DST(id)) {
        throw std::runtime_error("Not the 3DST");
    }
    int idb = (id>>9) & 0x000001FF;
    if (idb > 500) return -1;
    return idb;
}

int Cube::Info::CubePlane(int id) {
    if (!Cube::Info::Is3DST(id)) {
        throw std::runtime_error("Not the 3DST");
    }
    int idp = id & 0x000001FF;
    if (idp > 500) return -1;
    return idp;
}

bool Cube::Info::IsTPC(int id) {
    if (Cube::Info::SubDetector(id) == 25) return true;
    if (Cube::Info::SubDetector(id) == 26) return true;
    if (Cube::Info::SubDetector(id) == 27) return true;
    return false;
}

int Cube::Info::TPCNumber(int id) {
    if (!Cube::Info::IsTPC(id)) {
        throw std::runtime_error("Not a TPC");
    }
    return Cube::Info::SubDetector(id) - 25;
}

int Cube::Info::TPCAnode(int id) {
    if (!Cube::Info::IsTPC(id)) {
        throw std::runtime_error("Not a TPC");
    }
    int ida = (id >> 24) & 0x00000007;
    return ida;
}

int Cube::Info::TPCPadY(int id) {
    if (!Cube::Info::IsTPC(id)) {
        throw std::runtime_error("Not a TPC");
    }
    int idy = (id >> 12) & 0x00000FFF;
    return idy;
}

int Cube::Info::TPCPadZ(int id) {
    if (!Cube::Info::IsTPC(id)) {
        throw std::runtime_error("Not a TPC");
    }
    int idz = 0x00000FFF;
    return idz;
}

bool Cube::Info::IsECal(int id) {
    if (Cube::Info::SubDetector(id) == 12) return true;
    return false;
}

int Cube::Info::ECalModule(int id) {
    if (!Cube::Info::IsECal(id)) {
        throw std::runtime_error("Not the ECal");
    }
    id = (id >> 21) & 0x0000003F;
    return id;
}

int Cube::Info::ECalLayer(int id) {
    if (!Cube::Info::IsECal(id)) {
        throw std::runtime_error("Not the ECal");
    }
    id = (id >> 17) & 0x00000F;
    return id;
}

int Cube::Info::ECalCell(int id) {
    if (!Cube::Info::IsECal(id)) {
        throw std::runtime_error("Not the ECal");
    }
    id = (id >> 9) & 0x00000FF;
    return id;
}

int Cube::Info::ECalEnd(int id) {
    if (!Cube::Info::IsECal(id)) {
        throw std::runtime_error("Not the ECal");
    }
    id = id & 0x00000001;
    return id;
}

int Cube::Info::IdentifierProjection(int id) {
    int proj = 0;
    if (!Cube::Info::Is3DST(id)) return proj;
    if (Cube::Info::CubeNumber(id) > -1) proj += Cube::Info::kXAxis;
    if (Cube::Info::CubeBar(id) > -1) proj += Cube::Info::kYAxis;
    if (Cube::Info::CubePlane(id) > -1) proj += Cube::Info::kZAxis;
    return proj;
}

int Cube::Info::Identifier3DST(int num, int bar, int pln) {
    if (num<0) num = 511+num;
    if (bar<0) bar = 511+bar;
    if (pln<0) pln = 511+pln;
    int id = 13;
    id = num + (id << 9);
    id = bar + (id << 9);
    id = pln + (id << 9);
    return id;
}

int Cube::Info::Combine3DST(int id1, int id2, int id3) {
    int ss = 13;
    int cc = -1;
    int bb = -1;
    int pp = -1;

    int s1 = Cube::Info::SubDetector(id1);
    int c1 = Cube::Info::CubeNumber(id1);
    int b1 = Cube::Info::CubeBar(id1);
    int p1 = Cube::Info::CubePlane(id1);
    int s2 = Cube::Info::SubDetector(id2);
    int c2 = Cube::Info::CubeNumber(id2);
    int b2 = Cube::Info::CubeBar(id2);
    int p2 = Cube::Info::CubePlane(id2);

    // Get the number, bar, and plane for the cube.
    if (c1 < 0) cc = c2;
    if (c2 < 0) cc = c1;
    if (c1 == c2) cc = c1;
    if (b1 < 0) bb = b2;
    if (b2 < 0) bb = b1;
    if (b1 == b2) bb = b1;
    if (p1 < 0) pp = p2;
    if (p2 < 0) pp = p1;
    if (p1 == p2) pp = p1;

    // Sanity checking!
    if (ss != 13) return 0;
    if (s1 != ss) return 0;
    if (s2 != ss) return 0;
    if (cc < 0) return 0;
    if (bb < 0) return 0;
    if (pp < 0) return 0;

    // build the id.
    int sid = Cube::Info::Identifier3DST(cc,bb,pp);

    // The third id wasn't provided, so return the current id.
    if (id3 < 0) return sid;

    // The third id was provided, so make sure it's consistent.
    int s3 = Cube::Info::SubDetector(id3);
    int c3 = Cube::Info::CubeNumber(id3);
    int b3 = Cube::Info::CubeBar(id3);
    int p3 = Cube::Info::CubePlane(id3);

    if (s3 != ss) return 0;

    if (c3 < 0) {
        if (bb != b3) return 0;
        if (pp != p3) return 0;
    }

    if (b3 < 0) {
        if (cc != c3) return 0;
        if (pp != p3) return 0;
    }

    if (p3 < 0) {
        if (cc != c3) return 0;
        if (bb != b3) return 0;
    }

    return sid;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
