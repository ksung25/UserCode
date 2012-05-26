#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

RUNLIST=Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-02
ELE2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/merged_${RUNLIST}.root
MU2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012APromptV1/merged_${RUNLIST}.root

#
# options
#

MuonFakeRate=2
ElectronFakeRate=4

#
# electrons
#

ElectronFO=32
ElectronFO15=320
ElectronFO50=321
ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronIsoNumerator2011=400
ElectronIDNumerator2011=401
ElectronIDIsoNumerator2011=402

# default jet threshold 35
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FR_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_FR_IDISO/extra\",\"null\",\"HWW_Electron2012A_FR_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron2012A_FR_IDISO/extra\",\"HWW_Electron2012A_FR_IDISO/eff.root\",\"HWW_Electron2012A_FR_IDISO/eff.root\"\)

# syst jet 15
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FRJet15_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO15,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_FRJet15_IDISO/extra\",\"null\",\"HWW_Electron2012A_FRJet15_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron2012A_FRJet15_IDISO/extra\",\"HWW_Electron2012A_FRJet15_IDISO/eff.root\",\"HWW_Electron2012A_FRJet15_IDISO/eff.root\"\)

# syst jet 50
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FRJet50_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO50,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_FRJet50_IDISO/extra\",\"null\",\"HWW_Electron2012A_FRJet50_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron2012A_FRJet50_IDISO/extra\",\"HWW_Electron2012A_FRJet50_IDISO/eff.root\",\"HWW_Electron2012A_FRJet50_IDISO/eff.root\"\)

#
# muons
#

MuonFO=12
MuonFO5=120
MuonFO30=121
MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonIsoNumerator2011=200
MuonIDNumerator2011=201
MuonIDIsoNumerator2011=202

# default jet threshold 15
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FR_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_FR_IDISO/extra\",\"null\",\"HWW_Muon2012A_FR_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon2012A_FR_IDISO/extra\",\"HWW_Muon2012A_FR_IDISO/eff.root\",\"HWW_Muon2012A_FR_IDISO/eff.root\"\)

# syst jet 5
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FRJet5_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO5,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_FRJet5_IDISO/extra\",\"null\",\"HWW_Muon2012A_FRJet5_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon2012A_FRJet5_IDISO/extra\",\"HWW_Muon2012A_FRJet5_IDISO/eff.root\",\"HWW_Muon2012A_FRJet5_IDISO/eff.root\"\)

# syst jet 30
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FRJet30_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO30,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_FRJet30_IDISO/extra\",\"null\",\"HWW_Muon2012A_FRJet30_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon2012A_FRJet30_IDISO/extra\",\"HWW_Muon2012A_FRJet30_IDISO/eff.root\",\"HWW_Muon2012A_FRJet30_IDISO/eff.root\"\)

#rm *.so *.d

