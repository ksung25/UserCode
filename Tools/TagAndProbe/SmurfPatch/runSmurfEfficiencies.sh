#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

RUNLIST=Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-02
ELE2012A=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/merged_${RUNLIST}.root
MU2012A=/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/merged_${RUNLIST}.root
MC=/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged.root

#
# options
#

MuonTagAndProbe=1
MuonFakeRate=2
ElectronTagAndProbe=3
ElectronFakeRate=4
MuonTagAndProbeMC=5
ElectronTagAndProbeMC=6

#
# electrons
#

ElectronIsoDenominator2012=30
ElectronIDDenominator2012=31
ElectronIDIsoDenominator=33
ElectronIDIsoRndTagDenominator=34
ElectronFO=32

ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronTrigSglNumerator2012=43
ElectronTrigDblLeadNumerator2012=44
ElectronTrigDblTrailNumerator2012=45


# ID wrt Iso
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,2,2,2,\"${ELE2012A}\",\"HWW_Electron2012A_NM1Eff_ID\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator2012,$ElectronIDNumerator2012\,60,120,\"$MC\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MC}\",\"HWW_ElectronMC_NM1Eff_ID\",\"png\",1,4,0,$ElectronTagAndProbeMC,$ElectronIsoDenominator2012,$ElectronIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron2012A_NM1Eff_ID/extra\",\"HWW_ElectronMC_NM1Eff_ID/eff.root\",\"HWW_Electron2012A_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron2012A_NM1Eff_ID/extra\",\"HWW_ElectronMC_NM1Eff_ID/eff.root\",\"HWW_Electron2012A_NM1Eff_ID/eff.root\"\)

# Iso wrt ID
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,2,2,2,\"${ELE2012A}\",\"HWW_Electron2012A_NM1Eff_Iso\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDDenominator2012,$ElectronIsoNumerator2012\,60,120,\"$MC\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MC}\",\"HWW_ElectronMC_NM1Eff_Iso\",\"png\",1,4,0,$ElectronTagAndProbeMC,$ElectronIDDenominator2012,$ElectronIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron2012A_NM1Eff_Iso/extra\",\"HWW_ElectronMC_NM1Eff_Iso/eff.root\",\"HWW_Electron2012A_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron2012A_NM1Eff_Iso/extra\",\"HWW_ElectronMC_NM1Eff_Iso/eff.root\",\"HWW_Electron2012A_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"HWW_Electron2012A_NM1Eff_TrigSgl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoDenominator,$ElectronTrigSglNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Electron2012A_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"HWW_Electron2012A_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblLeadNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Electron2012A_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"HWW_Electron2012A_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblTrailNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Electron2012A_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Electron2012A_NM1Eff_TrigTrailDbl/eff.root\"\)


#
# muons
#

MuonIsoDenominator2012=10
MuonIDDenominator2012=11
MuonFO=12
MuonIDIsoDenominator=13
MuonIDIsoRndTagDenominator=14

MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonTrigSglNumerator2012=23
MuonTrigDblLeadNumerator2012=24
MuonTrigDblTrailNumerator2012=25

MC=/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/merged.root

#root -b -q plotEff.C+\(\"MuonOffline.bins\",1,2,1,2,\"${MU2012A}\",\"HWW_Muon2012A_NM1Eff_ID\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIsoDenominator2012,$MuonIDNumerator2012\,60,120,\"$MC\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MC}\",\"HWW_MuonMC_NM1Eff_ID\",\"png\",1,4,0,$MuonTagAndProbeMC,$MuonIsoDenominator2012,$MuonIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon2012A_NM1Eff_ID/extra\",\"HWW_MuonMC_NM1Eff_ID/eff.root\",\"HWW_Muon2012A_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon2012A_NM1Eff_ID/extra\",\"HWW_MuonMC_NM1Eff_ID/eff.root\",\"HWW_Muon2012A_NM1Eff_ID/eff.root\"\)

#root -b -q plotEff.C+\(\"MuonOffline.bins\",1,2,1,2,\"${MU2012A}\",\"HWW_Muon2012A_NM1Eff_Iso\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDDenominator2012,$MuonIsoNumerator2012\,60,120,\"$MC\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MC}\",\"HWW_MuonMC_NM1Eff_Iso\",\"png\",1,4,0,$MuonTagAndProbeMC,$MuonIDDenominator2012,$MuonIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon2012A_NM1Eff_Iso/extra\",\"HWW_MuonMC_NM1Eff_Iso/eff.root\",\"HWW_Muon2012A_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon2012A_NM1Eff_Iso/extra\",\"HWW_MuonMC_NM1Eff_Iso/eff.root\",\"HWW_Muon2012A_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"HWW_Muon2012A_NM1Eff_TrigSgl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoDenominator,$MuonTrigSglNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Muon2012A_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"HWW_Muon2012A_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblLeadNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Muon2012A_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"HWW_Muon2012A_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblTrailNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon2012A_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Muon2012A_NM1Eff_TrigTrailDbl/eff.root\"\)

#rm *.so *.d

