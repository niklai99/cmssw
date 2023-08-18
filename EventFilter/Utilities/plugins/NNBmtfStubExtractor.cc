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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

// track extrapolation
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

//
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFParameters.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFMasks.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFMasksRcd.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "CondFormats/L1TObjects/interface/L1TMuonBarrelParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonBarrelParamsRcd.h"






// using statements
using reco::MuonCollection;
using l1t::MuonBxCollection;
using l1t::RegionalMuonCandBxCollection;






// ----------------------------- CLASS DECLARATION  ----------------------------
class NNBmtfStubExtractor : public edm::stream::EDAnalyzer<> {

    public:
        // constructor and destructor
        explicit NNBmtfStubExtractor(const edm::ParameterSet&);
        ~NNBmtfStubExtractor() override {};

        // method for analyzing the events
        void analyze(const edm::Event&, const edm::EventSetup&) override;


    private:

        // conversion methods from hardware to physical coordinates
        int calcGlobalPhi(const l1t::RegionalMuonCand&) const;
        double calcPhysPhi(int) const;
        double calcPhysEta(int) const;
        double calcPtGeV(unsigned int) const;

        double calcDr(const l1t::RegionalMuonCand&, const l1t::Muon*) const;
        double calcDr(const l1t::Muon*, const TrajectoryStateOnSurface*) const;
        bool isBarrel(const l1t::Muon*) const;

        // muon EDM Tokens
        edm::EDGetTokenT<MuonCollection>                recoMuonsToken_;
        edm::EDGetTokenT<MuonBxCollection>              simGmtStage2DigisToken_;
        edm::EDGetTokenT<RegionalMuonCandBxCollection>  simKBmtfDigisToken_;
        edm::EDGetTokenT<L1MuKBMTrackBxCollection>      simKBmtfTracksToken_;
        // edm::EDGetTokenT<L1MuDTChambPhContainer>        L1MuDTChambPhsToken_;
        // edm::EDGetTokenT<L1MuDTChambThContainer>        L1MuDTChambThsToken_;

        // magnetic field
        //edm::ESHandle<MagneticField> theMagField;
        const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
        

        // const edm::ESGetToken<L1TMuonBarrelParams, L1TMuonBarrelParamsRcd> bmtfParamsToken_;

        // stub processor
        L1TMuonBarrelKalmanStubProcessor stubProcessor_;

        // propagators
        PropagateToMuonSetup const muPropagator1stSetup_;
        PropagateToMuonSetup const muPropagator2ndSetup_;

        // cuts and constants
        double drCut_;      // For muon matching
        double phiMult_;
        double etaMult_;
        std::vector<int> cotTheta_1_;
        std::vector<int> cotTheta_2_;
        std::vector<int> cotTheta_3_;


        // the min and max BX to be analyzed
        int minBx_;
        int maxBx_;

        // output files
        string filenameNNBmtfStubs_;
        std::ofstream fileNNBmtfStubs_;

        // debug and verbose flags
        bool debug_;
        bool verbose_;

};
// -----------------------------------------------------------------------------





// -------------------------------- constructor  -------------------------------
NNBmtfStubExtractor::NNBmtfStubExtractor(const edm::ParameterSet& iConfig):
    recoMuonsToken_(consumes<MuonCollection>(iConfig.getParameter<edm::InputTag>("recoMuons"))),                                    // reco muons token
    simGmtStage2DigisToken_(consumes<MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1tMuonBXVector"))),                    // emulated l1t GMT BMTF muon token
    simKBmtfDigisToken_(consumes<RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("l1tRegionalMuonCandBXVector"))),// emulated l1t BMTF regional candidate token
    simKBmtfTracksToken_(consumes<L1MuKBMTrackBxCollection>(iConfig.getParameter<edm::InputTag>("L1MuKBMTrackBXVector"))),          // emulated l1t KBMTF track token
    magfieldToken_(consumesCollector().esConsumes()),
    muPropagator1stSetup_(iConfig.getParameter<edm::ParameterSet>("muProp1st"),consumesCollector()),                                // reco muon propagator to St1
    muPropagator2ndSetup_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"),consumesCollector()),                                // reco muon propagator to St2
    drCut_(iConfig.getParameter<double>("drCut")),                                                                                  // DR cut
    phiMult_(iConfig.getParameter<double>("phiMult")),                                                                              // GMT muon phi conversion
    etaMult_(iConfig.getParameter<double>("etaMult")),                                                                              // GMT muon eta conversion
    minBx_(iConfig.getParameter<int>("minBx")),                                                                                     // the min BX to be analyzed [0]
    maxBx_(iConfig.getParameter<int>("maxBx")),                                                                                     // the max BX to be analyzed [3564]
    filenameNNBmtfStubs_(iConfig.getParameter<string>("filenameNNBmtfStubs")),                                                      // csv filename
    fileNNBmtfStubs_(filenameNNBmtfStubs_),
    debug_(iConfig.getParameter<bool>("debug")),                                                                                    // debug flag (currently not used)
    verbose_(iConfig.getParameter<bool>("verbose"))                                                                                 // verbose flag (currently not used)
{

    fileNNBmtfStubs_ << "run,"
                     << "ls,"
                     << "event,"
                     << "ptReco,"
                     << "etaExtRecoSt2,"
                     << "phiExtRecoSt2,"
                     << "chargeReco,"
                     << "etaVtxReco,"
                     << "phiVtxReco,"
                     << "dXYReco,"
                     << "isTight,"
                     // L1 quantities
                     << "ptL1,"
                     << "pt2L1,"
                     << "hwDXYL1,"
                     << "phiL1,"
                     << "etaL1,"
                     << "hwSignL1,"
                     << "hwSignValidL1,"
                     << "hwQualityL1,"
                     //    << "hwPhiAtVtx,"
                     //    << "hwEtaAtVtx,"
                     //    << "phiAtVtx,"
                     //    << "etaAtVtx,"
                     //    << "tfMuonIndex,"
                     // KBMTF quantities
                     << "kbmtfRVtx,kbmtfPhiVtx,kbmtfRMu,kbmtfPhiMu,kbmtfPhiBendMu,kbmtfFloatPTUnconstrained,kbmtfDXY,kbmtfApproxChi2,kbmtfCoarseEta,kbmtfRank,kbmtfFineEta,kbmtfHasFineEta,"
                     << "kbmtfSector,kbmtfHitPattern,kbmtfQuality,kbmtfWheel,kbmtfStep,"
                     // GMT quantities
                     << "GmtHwEtaAtVtx,GmtHwPhiAtVtx,GmtEtaAtVtx,GmtPhiAtVtx,GmtPt,GmtPtU,GmtEta,GmtPhi,"
                     << "n_stubs,"
                     << "s1_stNum,s1_scNum,s1_whNum,s1_eta,s1_phi,s1_phiB,s1_quality,"
                     << "s2_stNum,s2_scNum,s2_whNum,s2_eta,s2_phi,s2_phiB,s2_quality,"
                     << "s3_stNum,s3_scNum,s3_whNum,s3_eta,s3_phi,s3_phiB,s3_quality,"
                     << "s4_stNum,s4_scNum,s4_whNum,s4_eta,s4_phi,s4_phiB,s4_quality,"
                     << std::endl;

}
// -----------------------------------------------------------------------------





// ----------------------- method called for each orbit  -----------------------
void NNBmtfStubExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;

    //--------------------------------------------------------------------------
    //
    // Initializations
    //
    //--------------------------------------------------------------------------
    // muon propagators
    PropagateToMuon muPropagator1st_ = muPropagator1stSetup_.init(iSetup);  // propagator to 1st muon station
    PropagateToMuon muPropagator2nd_ = muPropagator2ndSetup_.init(iSetup);  // propagator to 2nd muon station

    // matched and unmatched muons vectors
    // std::vector<const ScouterMuon*>   vMatchedMuons;            // vector of matched muons in event
    // std::vector<const ScouterMuon*>   vMatchedMuonsWithKBMTF;   // vector of matched KBMTF muons in event
    // std::vector<const L1ScouterMuon*> vUnMatchedMuons;          // vector of UNmatched muons in event
    L1MuKBMTCombinedStubCollection    vCombinedStubs;           // vector of combined stubs in event

    // BMTF parameters
    // const L1TMuonBarrelParams& bmtfParams = iSetup.getData(bmtfParamsToken_);
    // L1MuDTTFMasks msks = bmtfParams.l1mudttfmasks;

    //
    BXVector<l1t::Muon> vGmtMuons = iEvent.get(simGmtStage2DigisToken_);
    vector<reco::Muon> vRecoMuons = iEvent.get(recoMuonsToken_);
    // vUGMTmuons.setBXRange(0,0);

    auto const theMagField = iSetup.getHandle(magfieldToken_);
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    //
    // Match KBMTF
    //
    //--------------------------------------------------------------------------
    // loop over reco muons
    int ii = -1;
    for (const auto& recoMuon : vRecoMuons) {

        ++ii;

        double l1dr_min = 100.0;
        int l1match_jj = -1;
        int l1match_kk = -1;
        int l1match_ll = -1;

        // create FreeTrajectoryState for propagation
        FreeTrajectoryState* ftsRecoMuon = new FreeTrajectoryState(
            GlobalPoint(
                recoMuon.vx(),
                recoMuon.vy(),
                recoMuon.vz()
            ),
            GlobalVector(
                recoMuon.px(),
                recoMuon.py(),
                recoMuon.pz()
            ),
            recoMuon.charge(),
            theMagField.product()
        );

        // std::cout << recoMuon.vx() << " " << recoMuon.vy() << " " << recoMuon.vz() << " " << recoMuon.px() << " " << recoMuon.py() << " " << recoMuon.pz() << std::endl;
        

        // propagate to St1 and St2
        TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*ftsRecoMuon);
        TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*ftsRecoMuon);

        // continue only if reco track is global muon
        bool recoMuonCheck = true; // recoMuon.isGlobalMuon();
        if (recoMuonCheck) {
            bool tight_ = 1; // track.passed(3); // track.muonID("GlobalMuonPromptTight");

            // std::cout << "0" << std::endl;

            // if propagation is valid...
            if (stateAtMuSt1.isValid() && stateAtMuSt2.isValid()) {

            	// std::cout << "-1" << std::endl;

                // loop over KBMTF l1tRegionalMuonCands
                int jj = -1;
                for (const auto& l1tRegMuon : iEvent.get(simKBmtfDigisToken_)) {
                    ++jj;

                    // std::cout << "--2" << std::endl;

                    // loop over KBMTF L1MuKBMTrack
                    int kk = -1;
                    for (const auto& l1tKBmtfTrack : iEvent.get(simKBmtfTracksToken_)) {
                        ++kk;

                        // std::cout << "---3" << std::endl;

                        // loop over uGMT muons
                        int ll = -1;
                        for(const auto& l1tGmtMuon : vGmtMuons) {
                            ++ll;

                            // std::cout << "----4" << std::endl;

                            // matching BMTF KF track with internal KF track + matching GMT KF track with BMTF kf track
                            bool fineEtaCheck   = (l1tKBmtfTrack.hasFineEta()==false) && (l1tRegMuon.hwEta()==l1tKBmtfTrack.coarseEta());
                            bool coarseEtaCheck = (l1tKBmtfTrack.hasFineEta()==true)  && (l1tRegMuon.hwEta()==l1tKBmtfTrack.fineEta());
                            bool etaCheck       = fineEtaCheck || coarseEtaCheck;
                            bool uncPtCheck     = l1tRegMuon.hwPtUnconstrained()==l1tGmtMuon.hwPtUnconstrained();
                            bool dxyCheck       = l1tRegMuon.hwDXY()==l1tGmtMuon.hwDXY();
                            if (etaCheck && uncPtCheck && dxyCheck) {

                                // std::cout << "----4 if" << std::endl;

                                // compute DR between l1t and reco muon candidates
                                double dr = calcDr(&l1tGmtMuon, &stateAtMuSt2);
                                if (dr < l1dr_min) {
                                    l1dr_min = dr;
                                    l1match_jj = jj;
                                    l1match_kk = kk;
                                    l1match_ll = ll;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (l1dr_min < drCut_) {
            const auto& l1tRegMuon    = iEvent.get(simKBmtfDigisToken_)[l1match_jj];
            const auto& l1tKBmtfTrack = iEvent.get(simKBmtfTracksToken_)[l1match_kk];
            const auto& l1tGmtMuon    = vGmtMuons[l1match_ll];


            fileNNBmtfStubs_    << iEvent.id().run() << "," << iEvent.luminosityBlock() << "," << iEvent.id().event() << ","
                                // RECO quantities
                                << recoMuon.pt() << ","                                            // "ptReco,"
                                << stateAtMuSt2.globalPosition().eta() << ","                      // "etaExtRecoSt2,"
                                << stateAtMuSt2.globalPosition().phi() << ","                      // "phiExtRecoSt2,"
                                << recoMuon.charge() << ","                                        // "chargeReco,"
                                << recoMuon.eta() << ","                                           // "etaVtxReco,"
                                << recoMuon.phi() << ","                                           // "phiVtxReco,"
                                << "-1" << ","                                                     // "dXYReco,"
                                << "1" << ","                                                      // "isTight,"
                                // L1 quantities
                                << calcPtGeV(l1tRegMuon.hwPt()) << ","                             // "ptL1,"
                                << calcPtGeV(l1tRegMuon.hwPtUnconstrained()) << ","                // "pt2L1,"
                                << l1tRegMuon.hwDXY() << ","                                       // "hwDXYL1,"
                                << calcPhysPhi(calcGlobalPhi(l1tRegMuon)) << ","                   // "phiL1,"
                                << calcPhysEta(l1tRegMuon.hwEta()) << ","                          // "etaL1,"
                                << l1tRegMuon.hwSign() << ","                                      // "hwSignL1,"
                                << l1tRegMuon.hwSignValid() << ","                                 // "hwSignValidL1,"
                                << l1tRegMuon.hwQual() << ","                                      // "hwQualityL1,"
                                // KBMTF quantities
                                << l1tKBmtfTrack.curvatureAtVertex() << ","                        // "kbmtfRVtx,"
                                << l1tKBmtfTrack.phiAtVertex() << ","                              // "kbmtfPhiVtx,"
                                << l1tKBmtfTrack.curvature() << ","                                // "kbmtfRMu,"
                                << l1tKBmtfTrack.positionAngle() << ","                            // "kbmtfPhiMu,"
                                << l1tKBmtfTrack.bendingAngle() << ","                             // "kbmtfPhiBendMu,"
                                << l1tKBmtfTrack.ptUnconstrained() << ","                          // "kbmtfFloatPTUnconstrained,"
                                << l1tKBmtfTrack.dxy() << ","                                      // "kbmtfDXY,"
                                << l1tKBmtfTrack.approxChi2() << ","                               // "kbmtfApproxChi2,"
                                << l1tKBmtfTrack.coarseEta() << ","                                // "kbmtfCoarseEta,"
                                << l1tKBmtfTrack.rank() << ","                                     // "kbmtfRank,"
                                << l1tKBmtfTrack.fineEta() << ","                                  // "kbmtfFineEta,"
                                << l1tKBmtfTrack.hasFineEta() << ","                               // "kbmtfHasFineEta,"
                                << l1tKBmtfTrack.sector() << ","                                   // "kbmtfSector,"
                                << l1tKBmtfTrack.hitPattern() << ","                               // "kbmtfHitPattern,"
                                << l1tKBmtfTrack.quality() << ","                                  // "kbmtfQuality,"
                                << l1tKBmtfTrack.wheel() << ","                                    // "kbmtfWheel,"
                                << l1tKBmtfTrack.step() << ","                                     // "kbmtfStep,"
                                //  << l1tKBmtfTrack.curvatureAtMuon() << ","
                                //  << l1tKBmtfTrack.phiAtMuon() << ","
                                //  << l1tKBmtfTrack.phiBAtMuon() << ","
                                // GMT quantities
                                << l1tGmtMuon.hwEtaAtVtx() << ","                                  // "GmtHwEtaAtVtx,"
                                << l1tGmtMuon.hwPhiAtVtx() << ","                                  // "GmtHwPhiAtVtx,"
                                << l1tGmtMuon.etaAtVtx() << ","                                    // "GmtEtaAtVtx,"
                                << l1tGmtMuon.phiAtVtx() << ","                                    // "GmtPhiAtVtx,"
                                << calcPtGeV(l1tGmtMuon.hwPt()) << ","                             // "GmtPt,"
                                << calcPtGeV(l1tGmtMuon.hwPtUnconstrained()) << ","                // "GmtPtU,"
                                << calcPhysEta(l1tGmtMuon.hwEta()) << ","                          // "GmtEta,"
                                << calcPhysPhi(l1tGmtMuon.hwPhi()) << ",";                         // "GmtPhi,"


            // write the number of stubs
            fileNNBmtfStubs_ << l1tKBmtfTrack.stubs().size() << ",";


            // STUBS INFO
            for (const auto& stub : l1tKBmtfTrack.stubs()) {
                fileNNBmtfStubs_ << stub->stNum() << ","
                                 << stub->scNum() << ","
                                 << stub->whNum() << ",";

                if (stub->tag()==0) {
                    fileNNBmtfStubs_ << stub->eta1() << ",";
                } else {
                    fileNNBmtfStubs_ << stub->eta2() << ",";
                }

                fileNNBmtfStubs_ << stub->phi() << ","
                                 << stub->phiB() << ","
                                 << stub->quality() << ",";
            }

            // fill dummy stubs if there are less than 4 stubs
            for (size_t i = 0; i < (4-l1tKBmtfTrack.stubs().size()); ++i) {
                fileNNBmtfStubs_ << "-1,-1,-3,-1,0,0,0,";
            }

            fileNNBmtfStubs_ << std::endl;
        }
    }
    //--------------------------------------------------------------------------

}
// -----------------------------------------------------------------------------





// ------------------------------ utility methods  -----------------------------

// method that computes the globalPhi coordinate of the muon in the hardware units
int NNBmtfStubExtractor::calcGlobalPhi(const l1t::RegionalMuonCand& muon) const {
    int globalPhi = muon.processor()*48 + muon.hwPhi();
    globalPhi += 576 - 24;
    globalPhi = globalPhi % 576;
    return globalPhi;
}

// method that computes the physical phi coordinate of the muon
double NNBmtfStubExtractor::calcPhysPhi(int globalPhi) const {
    double physPhi = 1.0 * globalPhi / phiMult_;
    if (physPhi > M_PI) {
        physPhi = physPhi - 2*M_PI;
    }
    return physPhi;
}

// method that computes the physical eta coordinate of the muon
double NNBmtfStubExtractor::calcPhysEta(int hwEta) const {
    double physEta = 1.0 * hwEta / etaMult_;
    return physEta;
}

// method that computes the physical pt of the muon
double NNBmtfStubExtractor::calcPtGeV(unsigned int pt) const {
    return pt * 0.5;
}

// method that computes the deltaR between a BMTF muon and a GMT muon
double NNBmtfStubExtractor::calcDr(const l1t::RegionalMuonCand& bmtfMuon, const l1t::Muon* gmtMuon) const {
    double dPhi = fabs(calcPhysPhi(calcGlobalPhi(bmtfMuon)) - calcPhysPhi(gmtMuon->hwPhi()));
    double dEta = fabs(calcPhysEta(bmtfMuon.hwEta()) - calcPhysEta(gmtMuon->hwEta()));
    return sqrt(dPhi*dPhi + dEta*dEta);
}

// method that computes the deltaR between a GMT muon and a reco track propagated to St2
double NNBmtfStubExtractor::calcDr(const l1t::Muon* l1tGmtMuon, const TrajectoryStateOnSurface* track) const {
    double l1tPhysPhi = l1tGmtMuon->hwPhi() / phiMult_;
    if (l1tPhysPhi > M_PI) {
        l1tPhysPhi = l1tPhysPhi - 2*M_PI;
    }

    double dphi = abs(l1tPhysPhi - track->globalPosition().phi());
    if (dphi > M_PI) {
        dphi = std::abs(dphi - 2*M_PI);
    }

    double deta = l1tGmtMuon->hwEta() / etaMult_ - track->globalPosition().eta();
    double dr = std::sqrt(dphi*dphi + deta*deta);

    return dr;
}

// method to check if the GMT muon is in the barrel to match it to a BMTF muon
bool NNBmtfStubExtractor::isBarrel(const l1t::Muon* gmtMuon) const {
    // the GMT muon is in the barrel it has tfIndex between 36 and 70 (included)
    return (gmtMuon->tfMuonIndex() >= 36) && (gmtMuon->tfMuonIndex() <= 70);
}

// // method to save reco muon data
// void NNBmtfStubExtractor::saveL1TRegMuonData(std::ofstream& ofile, const reco::Muon& recoMuon, ) { , const L1MuKBMTrack& l1tKBmtfTrack, const l1t::Muon& l1tGmtMuon) {
//     ofile <<
// }

// // method to save l1t regional muon candidate data
// void NNBmtfStubExtractor::saveL1TRegMuonData(std::ofstream& ofile, const l1t::RegionalMuonCand& l1tRegMuon) {
//     ofile <<
// }

// // method to save L1MuKBMTrack data
// void NNBmtfStubExtractor::saveL1MuKBMTrackData(std::ofstream& ofile, const L1MuKBMTrack& l1tKBmtfTrack) {
//     ofile <<
// }

// // method to save l1t gmt muon data
// void NNBmtfStubExtractor::saveL1TGmtMuonData(std::ofstream& ofile, const l1t::Muon& l1tGmtMuon) {
//     ofile <<
// }





// define this as a plug-in
DEFINE_FWK_MODULE(NNBmtfStubExtractor);
