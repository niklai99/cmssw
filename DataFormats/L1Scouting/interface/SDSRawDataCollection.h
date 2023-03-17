#ifndef L1Scouting_SRawDataCollection_h
#define L1Scouting_SRawDataCollection_h

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/Common/interface/traits.h"
#include "FWCore/Utilities/interface/GCCPrerequisite.h"


/** 
  *
  * This collection holds the raw data for all the 
  * scouting sources. It is a collection of FEDRawData
  *
  */

class SDSRawDataCollection: public edm::DoNotRecordParents {
  public:
    SDSRawDataCollection();

    virtual ~SDSRawDataCollection();

    // retrive data for the scouting source @params srcId
    const FEDRawData& FEDData(int sourceId) const;
    
    // retrive data for the scouting source @params srcId
    FEDRawData& FEDData(int sourceId);

    SDSRawDataCollection(const SDSRawDataCollection&);

    void swap(SDSRawDataCollection& other) { data_.swap(other.data_); }

  private:
    std::vector<FEDRawData> data_;  //< the raw data
};

inline void swap(SDSRawDataCollection& a, SDSRawDataCollection& b) { a.swap(b); }

#endif // L1Scouting_SRawDataCollection_h