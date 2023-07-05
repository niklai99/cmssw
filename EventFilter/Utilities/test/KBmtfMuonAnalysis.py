import FWCore.ParameterSet.Config as cms
import math

process = cms.Process( "DUMP" )

process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        enable = cms.untracked.bool(True)
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    )
)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:PoolOutputTest.root")
)

process.dump = cms.EDAnalyzer("KBmtfMuonAnalysis",
    gmtMuonInputTag = cms.InputTag("GmtUnpacker", "", "SCPU"),
    bmtfMuonInputTag = cms.InputTag("simKBmtfDigis", "BMTF", "SCPU"),
    drCut = cms.double(1.0),
    phiMult = cms.double(576./(2*math.pi)),
    etaMult = cms.double(1./0.010875),
    minBx = cms.int32(0),
    maxBx = cms.int32(3564),
    debug = cms.bool(False),
    filenameResults = cms.string("results.csv")
)


process.p = cms.Path(
    process.dump
)