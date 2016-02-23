import FWCore.ParameterSet.Config as cms

#part_id = cms.untracked.vint32()
#
#for i in range(200):
#    part_id.append(cms.vint32(211));
#
#

generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(12.30),
        MinPt = cms.double(0.30),
        PartID = cms.vint32(211),
        MaxEta = cms.double(2.4),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.4),
        MinPhi = cms.double(-3.14159265359) ## in radians
    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('pi pt 0.3 12.3'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)
