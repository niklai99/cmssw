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

class SRDCollection: public edm::DoNotRecordParents {
  public:
    SRDCollection();

    virtual ~SRDCollection();

    // retrive data for the scouting source @params srcId
    const FEDRawData& FEDData(int sourceId) const;
    
    // retrive data for the scouting source @params srcId
    FEDRawData& FEDData(int sourceId);

    SRDCollection(const SRDCollection&);

    void swap(SRDCollection& other) { data_.swap(other.data_); }

    void setSourceSize(int sourceId, size_t size) {len_[sourceId] = size;}
    size_t getSourceSize(int sourceId) const {return len_[sourceId];}

  private:
    std::vector<FEDRawData> data_;  // vector of raw data
    std::vector<size_t> len_;       // (Real) length of each buffer
};

inline void swap(SRDCollection& a, SRDCollection& b) { a.swap(b); }

#endif // L1Scouting_SRawDataCollection_h