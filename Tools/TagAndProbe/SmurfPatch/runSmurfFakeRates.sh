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
ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronIsoNumerator2011=400
ElectronIDNumerator2011=401
ElectronIDIsoNumerator2011=402

#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FR_ISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIsoNumerator2012\)
#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FR_ID\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIDNumerator2012\)
#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_FR_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_FR_IDISO/extra\",\"null\",\"HWW_Electron2012A_FR_IDISO/eff.root\"\)
root -b -q plotFakeRate.C+\(\"HWW_Electron2012A_FR_IDISO/extra\",\"HWW_Electron2012A_FR_IDISO/eff.root\"\)

#2011
#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_2011FR_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIDIsoNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Electron2012A_2011FR_IDISO/extra\",\"HWW_Electron2012A_2011FR_IDISO/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron2012A_FR_IDISO/extra\",\"HWW_Electron2012A_2011FR_IDISO/eff.root\",\"HWW_Electron2012A_FR_IDISO/eff.root\"\)

#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_2011FR_ISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIsoNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Electron2012A_2011FR_ISO/extra\",\"HWW_Electron2012A_2011FR_ISO/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron2012A_FR_ISO/extra\",\"HWW_Electron2012A_2011FR_ISO/eff.root\",\"HWW_Electron2012A_FR_ISO/eff.root\"\)

#root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron2012A_2011FR_ID\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO,$ElectronIDNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Electron2012A_2011FR_ID/extra\",\"HWW_Electron2012A_2011FR_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron2012A_FR_ID/extra\",\"HWW_Electron2012A_2011FR_ID/eff.root\",\"HWW_Electron2012A_FR_ID/eff.root\"\)

#
# muons
#

MuonFO=12
MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonIsoNumerator2011=200
MuonIDNumerator2011=201
MuonIDIsoNumerator2011=202

#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FR_ISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIsoNumerator2012\)
#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FR_ID\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIDNumerator2012\)
#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_FR_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIDIsoNumerator2012\)
#root -b -q plotFakeRate.C+\(\"HWW_Muon2012A_FR_IDISO/extra\",\"HWW_Muon2012A_FR_IDISO/eff.root\"\)

# 2011
#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_2011FR_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIDIsoNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Muon2012A_2011FR_IDISO/extra\",\"HWW_Muon2012A_2011FR_IDISO/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon2012A_FR_IDISO/extra\",\"HWW_Muon2012A_2011FR_IDISO/eff.root\",\"HWW_Muon2012A_FR_IDISO/eff.root\"\)

#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_2011FR_ISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIsoNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Muon2012A_2011FR_ISO/extra\",\"HWW_Muon2012A_2011FR_ISO/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon2012A_FR_ISO/extra\",\"HWW_Muon2012A_2011FR_ISO/eff.root\",\"HWW_Muon2012A_FR_ISO/eff.root\"\)

#root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon2012A_2011FR_ID\",\"png\",1,0,0,$MuonFakeRate,$MuonFO,$MuonIDNumerator2011\)
#root -b -q plotFakeRate.C+\(\"HWW_Muon2012A_2011FR_ID/extra\",\"HWW_Muon2012A_2011FR_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon2012A_FR_ID/extra\",\"HWW_Muon2012A_2011FR_ID/eff.root\",\"HWW_Muon2012A_FR_ID/eff.root\"\)

#rm *.so *.d

