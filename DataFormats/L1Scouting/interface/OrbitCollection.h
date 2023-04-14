#ifndef DataFormats_L1Scouting_OrbitCollection_h
#define DataFormats_L1Scouting_OrbitCollection_h


#include "DataFormats/Common/interface/traits.h"
#include "FWCore/Utilities/interface/GCCPrerequisite.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include <vector>

namespace scoutingRun3 {

  template <class T>
  class OrbitCollection {
  
    public: 
      OrbitCollection(): bxData_(3565), nObjects_(0) {}

      void push_back(int bx, T &object) {
          bxData_[bx].push_back(object);
          nObjects_ ++;
      }

      void flatten() {
        index_.reserve(3565);
        flatData_.reserve(nObjects_);
        index_[0] = 0;
        int idx = 1;
        for (auto &bxVec: bxData_) {
            flatData_.insert(flatData_.end(), bxVec.begin(), bxVec.end());
            // increase index position
            index_[idx] = index_[idx-1] + bxVec.size();
            idx++;
        }
        //bxData_.clear();
      }

    private:
      std::vector<int> index_;
      std::vector<T> flatData_;
      mutable std::vector<std::vector<T>> bxData_;
      int nObjects_;
  };

  typedef OrbitCollection<l1t::Muon>   MuonOrbitCollection;
  typedef OrbitCollection<l1t::Jet>    JetOrbitCollection;
  typedef OrbitCollection<l1t::EGamma> EGammaOrbitCollection;
  typedef OrbitCollection<l1t::Tau>    TauOrbitCollection;
  typedef OrbitCollection<l1t::EtSum>  EtSumOrbitCollection;

}
#endif // DataFormats_L1Scouting_OrbitCollection_h
