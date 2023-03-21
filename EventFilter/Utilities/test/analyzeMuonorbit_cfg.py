import FWCore.ParameterSet.Config as cms

process = cms.Process( "DEMO" )

process.MessageLogger = cms.Service("MessageLogger",
  cerr = cms.untracked.PSet(
    enable = cms.untracked.bool(False)
  ),
  cout = cms.untracked.PSet(
    enable = cms.untracked.bool(True),
    threshold = cms.untracked.string('INFO')
  )
)


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:PoolOutputTest.root")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histodemo.root')
)

process.histo = cms.EDAnalyzer("DemoSCAnalyzer",
  muInputTag = cms.InputTag("GmtUnpacker"),
  minBx = cms.int32(0),
  maxBx = cms.int32(3565)
)

process.options.numberOfThreads = 4
process.options.numberOfStreams = 1

process.p = cms.Path(
  process.histo
)
