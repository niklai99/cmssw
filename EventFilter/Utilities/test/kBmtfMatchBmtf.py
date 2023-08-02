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
    fileNames = cms.untracked.vstring("file:kBMTFGMTdata.root")
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('BmtfDuplicates.root')
)

process.dump = cms.EDAnalyzer("kBmtfMatchBmtf",
    bmtfMuonTag     = cms.InputTag("simKBmtfDigis", "BMTF", "SCPU"),
    minBx           = cms.int32(0),
    maxBx           = cms.int32(3564),
    phiMult         = cms.double(576./(2*math.pi)),
    etaMult         = cms.double(1./0.010875),
    drThreshold     = cms.double(0.4),
    debug           = cms.bool(False),
    vverbose        = cms.bool(True)
)

process.p = cms.Path(
    process.dump
)