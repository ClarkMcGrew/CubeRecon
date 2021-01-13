#ifndef CubeInfo_hxx_seen
#define CubeInfo_hxx_seen

namespace Cube {
    class Info;
}

/// A collection of static functions and singleton classes that collect
/// detector related information.
class Cube::Info {
public:
    /// Return the sub-detector for the sensor id.  The sensor id has the
    /// following bit definitions.
    ///
    ///  DDDDD XXXXXXXXXXXXXXXXXXXXXXXXXXX
    ///
    /// DDDDD (31-27) -- The subdetector id.
    /// XX..X (26-0)  -- Subdetector specific.

    static int SubDetector(int i);

    /// @{ Return the number of the 3DST plane, bar or cube based on the
    /// sensor id.  If the id is not for the cube detector, this will throw an
    /// exception.  If the id is for a sensor, then the two relevant axis will
    /// be positive, and the sensor axis is -1.  If this is expanded for
    /// double ended readout, the "far" sensor will be "-2".  The "number"
    /// increments along the "X" axis.  The "bar" increments along the "Y"
    /// axis.  The "plane" increments along the "Z" axis.  The 3DST bit
    /// definitions are
    ///
    /// The 3DST sensor id is DDDDDCCCCCCCCCBBBBBBBBBPPPPPPPPP
    ///
    /// DDDDD     (31-27) -- The subdetector (3DST == 13).
    /// CCCCCCCCC (26-18) -- The 3DST cube
    /// BBBBBBBBB (17-9)  -- The 3DST bar
    /// PPPPPPPPP (8-0)   -- The 3DST plane
    ///
    static bool Is3DST(int id);
    static int CubeNumber(int id);
    static int CubeBar(int id);
    static int CubePlane(int id);
    /// @}


    /// @{ Return the number of the TPC anode (left or right) and pad (Y and
    /// Z) based on the sensor id.  If the id is not for the cube detector,
    /// this will throw an exception.  If the id is for a sensor, then the two
    /// relevant axis will be positive, and the sensor axis is -1.  If this is
    /// expanded for double ended readout, the "far" sensor will be "-2".  The
    /// "number" increments along the "X" axis.  The "bar" increments along
    /// the "Y" axis.  The "plane" increments along the "Z" axis.  The 3DST
    /// bit definitions are
    ///
    /// The TPC sensor id is DDDDDPPPYYYYYYYYYYYYZZZZZZZZZZZZ
    ///
    /// DDDDD        (31-27) -- The subdetector (TPC == 25, 26, 27).
    /// PPP          (26-24) -- The anode
    /// YYYYYYYYYYYY (23-12)  -- The TPC pad (Y)
    /// ZZZZZZZZZZZZ (11-0)   -- The TPC pad (Z)
    ///
    static bool IsTPC(int id);
    static int TPCNumber(int id);
    static int TPCAnode(int id);
    static int TPCPadY(int id);
    static int TPCPadZ(int id);

    /// @{ Return the information out of the ECal sensor id.  If the id is not
    /// for the ecal detector, this will throw an exception.
    ///
    /// The ECal sensor id is DDDDDMMMMMMLLLLCCCCCCCCXXXXXXXXE
    ///
    /// DDDDD        (31-27) -- The subdetector (TPC == 25, 26, 27).
    /// MMMMMM       (26-21) -- The module
    /// LLLL         (20-17) -- The layer
    /// CCCCCCCC     (16-9)  -- The cell
    /// XXXXXXXX     (8-1)   -- Reserved
    /// E            (0)     -- The end.
    static bool IsECal(int id);
    static int ECalModule(int id);
    static int ECalLayer(int id);
    static int ECalCell(int id);
    static int ECalEnd(int id);
    /// @}

    /// Return the projection for a hit based on the identifier.  The
    /// projections are X == 1 (0b001), Y == 2 (0b010), Z = 4 (0b100),
    /// YZ == 6 (0b110), XZ == 5 (0b101), XY == 3 (0b011),
    ///  or XYZ == 7 (0b111).
    enum { kXAxis = 1, kYAxis = 2, kZAxis = 4,
           kYZProj = 6, kXZProj = 5, kXYProj = 3,
           kXYZProj = 7};

    static int IdentifierProjection(int id1);

    /// Build a cube identifer.
    static int Identifier3DST(int num, int bar, int pln);

    /// Take two or three hit identifiers and return the identifier of where
    /// the intersect.  If the ids do not correspond to a valid intersection,
    /// then return 0.
    static int Combine3DST(int id1, int id2, int id3 = -1);


};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
