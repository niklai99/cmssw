#include "FWCore/Framework/interface/MakerMacros.h"

// system include files
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

// user include files
//   base class
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"


// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
// #include "TDirectory.h"



// ----------------------------- CLASS DECLARATION  ----------------------------
class kBmtfReadFromRaw : public edm::stream::EDAnalyzer<> {


    public:

        // constructor and destructor
        explicit kBmtfReadFromRaw(const edm::ParameterSet&); 
        ~kBmtfReadFromRaw() override{};

        // method for analyzing the events
        void analyze(const edm::Event&, const edm::EventSetup&) override;


    private:

        // structure to store the stubs variables
        struct BmtfStub {
            int bx, sc, st, wh, phi, phiB, qual, eta1, eta2;
            bool tag;
        };
        
        // structure to store the BMTF muon variables
        struct BmtfMuon {
            int bx, idx, dxy, qual, sign, signValid, processor;
            double pt, ptUnconstrained, phi, eta;
        };

        // structure to store the GMT muon variables
        struct GmtMuon {
            int bx, idx, tfidx, dxy, charge, chargeValid, quality;
            double pt, ptUnconstrained, phi, eta, phiAtVtx, etaAtVtx;
        };


        // methods to initialize the analysis type
        void initCollectData(bool, bool, bool);     // collectStubs, collectBmtfMuons, collectGmtMuons


        // conversion methods from hardware to physical coordinates
        int calcGlobalPhi(const l1t::RegionalMuonCand&) const;
        double calcPhysPhi(int) const;
        double calcPhysEta(int) const;
        double calcPtGeV(unsigned int) const;

        // methods to process the stubs and muons
        void processBmtfStubs(const edm::Handle<scoutingRun3::BmtfStubOrbitCollection>&, int, bool);
        void processBmtfMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection>&, int, bool);
        bool processGmtMuons(const edm::Handle<scoutingRun3::MuonOrbitCollection>&, int); 

        // select the GMT muons in the barrel
        bool isBarrel(const l1t::Muon*) const;



        // the tokens to access the data
        edm::EDGetTokenT<scoutingRun3::BmtfStubOrbitCollection> bmtfStubToken_;
        edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection>     bmtfMuonToken_;
        edm::EDGetTokenT<scoutingRun3::MuonOrbitCollection>     gmtMuonToken_;

        // the min and max BX to be analyzed
        int minBx_;
        int maxBx_;

        // phi and eta conversion factors
        double phiMult_;
        double etaMult_;

        // debug and verbose flags
        bool debug_;
        bool vverbose_;

        // the root file service to handle the output file
        edm::Service<TFileService> fs;

        // the histograms for the stubs 
        TH1D* h_stubOccupancyBX;
        TH1D* h_stubMultiplicityBX;

        // the histograms for the BMTF muons
        TH1D* h_bmtfMuonOccupancyBX;
        TH1D* h_bmtfMuonMultiplicityBX;

        // the histograms for the GMT muons
        TH1D* h_gmtMuonOccupancyBX;
        TH1D* h_gmtMuonMultiplicityBX;

        // the tree for the stubs
        TTree* stubsTree_;
        BmtfStub stubData_;

        // the tree for the BMTF muons
        TTree* bmtfMuonsTree_;
        BmtfMuon bmtfMuonData_;

        // the tree for the GMT muons
        TTree* gmtMuonsTree_;
        GmtMuon gmtMuonData_;

        // the orbit number
        unsigned int orbitNum_;

};
// -----------------------------------------------------------------------------



// -------------------------------- constructor  -------------------------------

kBmtfReadFromRaw::kBmtfReadFromRaw(const edm::ParameterSet& iConfig):
    bmtfStubToken_(consumes<scoutingRun3::BmtfStubOrbitCollection>(iConfig.getParameter<edm::InputTag>("bmtfStubTag"))),            // the token to access the stubs
    bmtfMuonToken_(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("bmtfMuonTag"))),                // the token to access the muons from the BMTF
    gmtMuonToken_(consumes<scoutingRun3::MuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtMuonTag"))),                  // the token to access the muons from the GMT
    minBx_(iConfig.getParameter<int>("minBx")),                                                                                     // the min BX to be analyzed [0]
    maxBx_(iConfig.getParameter<int>("maxBx")),                                                                                     // the max BX to be analyzed [3564]
    phiMult_(iConfig.getParameter<double>("phiMult")),                                                                              // the phi conversion factor from hardware to physical units 
    etaMult_(iConfig.getParameter<double>("etaMult")),                                                                              // the eta conversion factor from hardware to physical units
    debug_(iConfig.getParameter<bool>("debug")),                                                                                    // the debug flag (currently not used)
    vverbose_(iConfig.getParameter<bool>("vverbose"))                                                                               // the very verbose flag (currently not used)
{   

    bool collectStubs     = true;
    bool collectBmtfMuons = true;
    bool collectGmtMuons  = true;

    // initialize the analysis type
    initCollectData(collectStubs, collectBmtfMuons, collectGmtMuons);

    

}

// -----------------------------------------------------------------------------





// ----------------------- method called for each orbit  -----------------------
void kBmtfReadFromRaw::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

    // get the stubs
    edm::Handle<scoutingRun3::BmtfStubOrbitCollection> bmtfStubHandle; 
    iEvent.getByToken(bmtfStubToken_, bmtfStubHandle);

    // get the muons
    edm::Handle<l1t::RegionalMuonCandBxCollection> bmtfMuonHandle; // a BXVector of RegionalMuonCand
    iEvent.getByToken(bmtfMuonToken_, bmtfMuonHandle);

    // get the GMT muons
    edm::Handle<scoutingRun3::MuonOrbitCollection> gmtMuonHandle; 
    iEvent.getByToken(gmtMuonToken_, gmtMuonHandle);

    // update the orbit number
    orbitNum_++;

    // loop over the BXs
    for (int bx = minBx_; bx < maxBx_; ++bx) {

        // process the GMT muons
        bool thereIsMuon = processGmtMuons(gmtMuonHandle, bx);

        // process the kBMTF muons
        processBmtfMuons(bmtfMuonHandle, bx, thereIsMuon);

        // process the stubs 
        processBmtfStubs(bmtfStubHandle, bx, thereIsMuon);

        if (debug_ && thereIsMuon) {
            std::cout << "------------------------------------" << std::endl;
        }

    }

}
// -----------------------------------------------------------------------------



// ------------------------ methods called for each BX  ------------------------

// method to process stubs
void kBmtfReadFromRaw::processBmtfStubs(const edm::Handle<scoutingRun3::BmtfStubOrbitCollection>& bmtfStubHandle, int bx, bool thereIsMuon) {
    
    // get the index vector
    const std::vector<int>* bmtfStubIndex = bmtfStubHandle->getIndex();

    // get the number of stubs in this BX
    int current_bx_index = bmtfStubIndex->at(bx);
    int next_bx_index    = bmtfStubIndex->at(bx+1);
    int nStubs           = next_bx_index - current_bx_index;

    if (debug_ && thereIsMuon) {
        std::cout << std::endl << "BX = " << bx << std::endl;
        std::cout << "STUBS" << std::endl;
        std::cout << "nStubs = " << nStubs << std::endl;
    }

    // fill the histograms
    h_stubOccupancyBX->Fill(bx, nStubs);
    h_stubMultiplicityBX->Fill(nStubs);

    // if there are no stubs in this BX, return
    if (nStubs == 0) return;

    // loop over the stubs in this BX
    for (int iStub = current_bx_index; iStub < next_bx_index; ++iStub) {
        
        // get the stub
        const L1MuKBMTCombinedStub* stub = bmtfStubHandle->getFlatData(iStub);

        // fill the tree branches
        stubData_.bx   = bx;
        stubData_.sc   = stub->scNum();
        stubData_.st   = stub->stNum();
        stubData_.wh   = stub->whNum();
        stubData_.phi  = stub->phi();
        stubData_.phiB = stub->phiB();
        stubData_.qual = stub->quality();
        stubData_.eta1 = stub->eta1();
        stubData_.eta2 = stub->eta2();
        stubData_.tag  = stub->tag();

        if (debug_ && thereIsMuon) {
            std::cout << "ID = " << iStub-current_bx_index << std::endl;
            std::cout << "SC = " << stub->scNum() << " ST = " << stub->stNum() << " WH = " << stub->whNum() << std::endl;
            std::cout << "PHI = " << stub->phi() << " PHI_B = " << stub->phiB() << std::endl;
            std::cout << "QUAL = " << stub->quality() << std::endl;
            std::cout << "ETA1 = " << stub->eta1() << " ETA2 = " << stub->eta2() << std::endl;
        }

        // fill the tree
        stubsTree_->Fill();  
    }
}



// method to process muons
void kBmtfReadFromRaw::processBmtfMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection>& bmtfMuonHandle, int bx, bool thereIsMuon) {

    // get the number of muons in this BX
    int nMuons = bmtfMuonHandle->size(bx);

    if (debug_ && thereIsMuon) {
        std::cout << std::endl << "BX = " << bx << std::endl;
        std::cout << "BMTF MUONS" << std::endl;
        std::cout << "nMuons = " << nMuons << std::endl;
    }

    // fill the histograms
    h_bmtfMuonOccupancyBX->Fill(bx, nMuons);
    h_bmtfMuonMultiplicityBX->Fill(nMuons);

    // if there are no muons in this BX, return
    if (nMuons == 0) return;

    // loop over the muons in this BX
    for (int iMuon = 0; iMuon < nMuons; ++iMuon) {
    
        // get the muon
        const l1t::RegionalMuonCand& muon = bmtfMuonHandle->at(bx, iMuon);  

        // fill the tree branches
        bmtfMuonData_.bx              = bx;
        bmtfMuonData_.pt              = calcPtGeV(muon.hwPt());
        bmtfMuonData_.ptUnconstrained = calcPtGeV(muon.hwPtUnconstrained());
        bmtfMuonData_.dxy             = muon.hwDXY();
        bmtfMuonData_.phi             = calcPhysPhi(calcGlobalPhi(muon));
        bmtfMuonData_.eta             = calcPhysEta(muon.hwEta());
        bmtfMuonData_.qual            = muon.hwQual();
        bmtfMuonData_.sign            = muon.hwSign();
        bmtfMuonData_.signValid       = muon.hwSignValid();
        bmtfMuonData_.idx             = iMuon;
        bmtfMuonData_.processor       = muon.processor();

        if (debug_ && thereIsMuon) {
            std::cout << "ID = " << iMuon << std::endl;
            std::cout << "PT = " << calcPtGeV(muon.hwPt()) << " PT_UNCONSTR = " << calcPtGeV(muon.hwPtUnconstrained()) << std::endl;
            std::cout << "DXY = " << muon.hwDXY() << std::endl;
            std::cout << "PHI = " << calcPhysPhi(muon.hwPhi()) << " ETA = " << calcPhysEta(muon.hwEta()) << std::endl;
            std::cout << "QUAL = " << muon.hwQual() << std::endl;
            std::cout << "PROCESSOR = " << muon.processor() << std::endl;
        }

        // fill the tree
        bmtfMuonsTree_->Fill();  
    }
}



// method to process GMT muons (similar to the stubs one with the Index and FlatData)
bool kBmtfReadFromRaw::processGmtMuons(const edm::Handle<scoutingRun3::MuonOrbitCollection>& gmtMuonHandle, int bx) {

    // get the index vector 
    const std::vector<int>* gmtMuonIndex = gmtMuonHandle->getIndex();

    // get the number of muons in this BX
    int current_bx_index = gmtMuonIndex->at(bx);
    int next_bx_index    = gmtMuonIndex->at(bx+1);
    int nMuons           = next_bx_index - current_bx_index;

    // fill the histograms
    h_gmtMuonOccupancyBX->Fill(bx, nMuons);
    h_gmtMuonMultiplicityBX->Fill(nMuons);

    bool thereIsMuon = false;

    // if there are no muons in this BX, return
    if (nMuons == 0) return thereIsMuon;

    // loop over the muons in this BX
    for (int iMuon = current_bx_index; iMuon < next_bx_index; ++iMuon) {
    
        // get the muon
        const l1t::Muon* muon = gmtMuonHandle->getFlatData(iMuon);

        if (!isBarrel(muon)) continue;

        thereIsMuon = true;

        // fill the tree branches
        gmtMuonData_.bx              = bx;
        gmtMuonData_.pt              = calcPtGeV(muon->hwPt());
        gmtMuonData_.ptUnconstrained = calcPtGeV(muon->hwPtUnconstrained());
        gmtMuonData_.dxy             = muon->hwDXY();
        gmtMuonData_.phi             = calcPhysPhi(muon->hwPhi());
        gmtMuonData_.eta             = calcPhysEta(muon->hwEta());
        gmtMuonData_.phiAtVtx        = muon->phiAtVtx();
        gmtMuonData_.etaAtVtx        = muon->etaAtVtx();
        gmtMuonData_.charge          = muon->hwCharge();
        gmtMuonData_.chargeValid     = muon->hwChargeValid();
        gmtMuonData_.tfidx           = muon->tfMuonIndex();
        gmtMuonData_.idx             = iMuon - current_bx_index;
        gmtMuonData_.quality         = muon->hwQual();

        if (debug_ && thereIsMuon) {
            std::cout << std::endl << "BX = " << bx << std::endl;
            std::cout << "GMT MUONS" << std::endl;
            std::cout << "nMuons = " << nMuons << std::endl;
            std::cout << "ID = " << iMuon - current_bx_index << std::endl;
            std::cout << "PT = " << calcPtGeV(muon->hwPt()) << " PT_UNCONSTR = " << calcPtGeV(muon->hwPtUnconstrained()) << std::endl;
            std::cout << "DXY = " << muon->hwDXY() << std::endl;
            std::cout << "PHI = " << calcPhysPhi(muon->hwPhi()) << " ETA = " << calcPhysEta(muon->hwEta()) << std::endl;
            std::cout << "PHI_AT_VTX = " << muon->phiAtVtx() << " ETA_AT_VTX = " << muon->etaAtVtx() << std::endl;
            std::cout << "CHARGE = " << muon->hwCharge() << std::endl;
        }

        // fill the tree
        gmtMuonsTree_->Fill();  
    }

    return thereIsMuon;
}

// -----------------------------------------------------------------------------


// --------------------------- analysis initializers  --------------------------

// method to initialize the data collection
void kBmtfReadFromRaw::initCollectData(bool collectStubs, bool collectBmtfMuons, bool collectGmtMuons) {

    // initialize the orbit number
    orbitNum_ = 0;

    if (collectStubs) {
        // initialize the histograms for the stubs
        h_stubOccupancyBX    = fs->make<TH1D>("h_stubOccupancyBX",    "h_stubOccupancyBX",    maxBx_ - minBx_, minBx_, maxBx_);
        h_stubMultiplicityBX = fs->make<TH1D>("h_stubMultiplicityBX", "h_stubMultiplicityBX", 50, 0, 50);

         // initialize the tree and call it "stubs"
        stubsTree_ = fs->make<TTree>("stubs", "stubs");

        // initialize the tree branches
        stubsTree_->Branch("stubBx",   &stubData_.bx,   "stubBx/I");
        stubsTree_->Branch("stubSc",   &stubData_.sc,   "stubSc/I");
        stubsTree_->Branch("stubSt",   &stubData_.st,   "stubSt/I");
        stubsTree_->Branch("stubWh",   &stubData_.wh,   "stubWh/I");
        stubsTree_->Branch("stubPhi",  &stubData_.phi,  "stubPhi/I");
        stubsTree_->Branch("stubPhiB", &stubData_.phiB, "stubPhiB/I");
        stubsTree_->Branch("stubQual", &stubData_.qual, "stubQual/I");
        stubsTree_->Branch("stubEta1", &stubData_.eta1, "stubEta1/I");
        stubsTree_->Branch("stubEta2", &stubData_.eta2, "stubEta2/I");
        stubsTree_->Branch("stubTag",  &stubData_.tag,  "stubTag/O");

        stubsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");

    }

    if (collectBmtfMuons) {
        // initialize the histograms for the muons
        h_bmtfMuonOccupancyBX    = fs->make<TH1D>("h_bmtfMuonOccupancyBX",    "h_bmtfMuonOccupancyBX",    maxBx_ - minBx_, minBx_, maxBx_);
        h_bmtfMuonMultiplicityBX = fs->make<TH1D>("h_bmtfMuonMultiplicityBX", "h_bmtfMuonMultiplicityBX", 50, 0, 50);

        // initialize the tree and call it "bmtfMuons"
        bmtfMuonsTree_ = fs->make<TTree>("bmtfMuons", "bmtfMuons");

        // initialize the tree branches
        bmtfMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
        bmtfMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
        bmtfMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
        bmtfMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
        bmtfMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
        bmtfMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
        bmtfMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
        bmtfMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
        bmtfMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
        bmtfMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
        bmtfMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");

        bmtfMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");
    }

    
    if (collectGmtMuons) {
        // initialize the histograms for the GMT muons
        h_gmtMuonOccupancyBX    = fs->make<TH1D>("h_gmtMuonOccupancyBX",    "h_gmtMuonOccupancyBX",    maxBx_ - minBx_, minBx_, maxBx_);
        h_gmtMuonMultiplicityBX = fs->make<TH1D>("h_gmtMuonMultiplicityBX", "h_gmtMuonMultiplicityBX", 50, 0, 50);

        // initialize the tree and call it "gmtMuons"
        gmtMuonsTree_ = fs->make<TTree>("gmtMuons", "gmtMuons");

        // initialize the tree branches
        gmtMuonsTree_->Branch("gmtMuonBx",              &gmtMuonData_.bx,              "gmtMuonBx/I");
        gmtMuonsTree_->Branch("gmtMuonIdx",             &gmtMuonData_.idx,             "gmtMuonIdx/I");
        gmtMuonsTree_->Branch("gmtMuonPt",              &gmtMuonData_.pt,              "gmtMuonPt/D");
        gmtMuonsTree_->Branch("gmtMuonPtUnconstrained", &gmtMuonData_.ptUnconstrained, "gmtMuonPtUnconstrained/D");
        gmtMuonsTree_->Branch("gmtMuonDxy",             &gmtMuonData_.dxy,             "gmtMuonDxy/I");
        gmtMuonsTree_->Branch("gmtMuonPhi",             &gmtMuonData_.phi,             "gmtMuonPhi/D");
        gmtMuonsTree_->Branch("gmtMuonEta",             &gmtMuonData_.eta,             "gmtMuonEta/D");
        gmtMuonsTree_->Branch("gmtMuonPhiAtVtx",        &gmtMuonData_.phiAtVtx,        "gmtMuonPhiAtVtx/D");
        gmtMuonsTree_->Branch("gmtMuonEtaAtVtx",        &gmtMuonData_.etaAtVtx,        "gmtMuonEtaAtVtx/D");
        gmtMuonsTree_->Branch("gmtMuonTfIdx",           &gmtMuonData_.tfidx,           "gmtMuonTfIdx/I");
        gmtMuonsTree_->Branch("gmtMuonCharge",          &gmtMuonData_.charge,          "gmtMuonCharge/I");
        gmtMuonsTree_->Branch("gmtMuonChargeValid",     &gmtMuonData_.chargeValid,     "gmtMuonChargeValid/I");
        gmtMuonsTree_->Branch("gmtMuonQuality",         &gmtMuonData_.quality,         "gmtMuonQuality/I");
        
        gmtMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");

    }

}

// -----------------------------------------------------------------------------




// ------------------------------ utility methods  -----------------------------

// method that computes the globalPhi coordinate of the muon in the hardware units
int kBmtfReadFromRaw::calcGlobalPhi(const l1t::RegionalMuonCand& muon) const {
    int globalPhi = muon.processor()*48 + muon.hwPhi();
    globalPhi += 576 - 24;
    globalPhi = globalPhi % 576;
    return globalPhi;
}

// method that computes the physical phi coordinate of the muon
double kBmtfReadFromRaw::calcPhysPhi(int globalPhi) const {
    double physPhi = 1.0 * globalPhi / phiMult_;
    if (physPhi > M_PI) {
        physPhi = physPhi - 2*M_PI;
    }
    return physPhi;
}

// method that computes the physical eta coordinate of the muon
double kBmtfReadFromRaw::calcPhysEta(int hwEta) const {
    double physEta = 1.0 * hwEta / etaMult_;
    return physEta;
}

// method that computes the physical pt of the muon
double kBmtfReadFromRaw::calcPtGeV(unsigned int pt) const {
    return pt * 0.5;
}


// method to check if the GMT muon is in the barrel to match it to a BMTF muon
bool kBmtfReadFromRaw::isBarrel(const l1t::Muon* gmtMuon) const {
    // the GMT muon is in the barrel it has tfIndex between 36 and 70 (included)
    return (gmtMuon->tfMuonIndex() >= 36 && gmtMuon->tfMuonIndex() <= 70);
}

// -----------------------------------------------------------------------------


DEFINE_FWK_MODULE(kBmtfReadFromRaw);


