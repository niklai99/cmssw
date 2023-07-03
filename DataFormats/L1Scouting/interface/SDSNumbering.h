#ifndef L1Scouting_SDSNumbering_h
#define L1Scouting_SDSNumbering_h

/** 
  *
  * This class holds the Scouting Data Source (SDS)
  * numbering scheme for the Level 1 scouting system
  *
  */

class SDSNumbering {
  public:
    static constexpr int lastSDSId() { return MAXSDSID; }

    enum {
      NOT_A_SDSID = -1,
      MAXSDSID = 128,
      GmtSDSID = 0x1,
      CaloSDSID = 0x2,
      GtSDSID = 0x4
    };
    };
};

#endif // L1Scouting_SDSNumbering_h
