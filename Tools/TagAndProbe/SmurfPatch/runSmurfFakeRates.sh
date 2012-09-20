#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

#RUNLIST=Cert_190456-195658_8TeV_PromptReco_Collisions12_JSON
#TAG=V00-02-04
#ELE2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012/merged_${RUNLIST}.root
#MU2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012/merged_${RUNLIST}.root
#OUTNAME="2012June12V0"

RUNLIST=Cert_190456-196509_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-06
ELE2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012/merged_${RUNLIST}.root
MU2012=/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012/merged_${RUNLIST}.root
OUTNAME="2012June22V0"


#
# options
#

MuonFakeRate=2
ElectronFakeRate=4

#
# electrons
#

ElectronFO15=320
ElectronFO20=321
ElectronFO25=322
ElectronFO30=323
ElectronFO35=324
ElectronFO40=325
ElectronFO45=326
ElectronFO50=327
ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronIsoNumerator2011=400
ElectronIDNumerator2011=401
ElectronIDIsoNumerator2011=402

# syst jet 15
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet15_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO15,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet15_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet15_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet15_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet15_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet15_IDISO/eff.root\"\)

# syst jet 20
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet20_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO20,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet20_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet20_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet20_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet20_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet20_IDISO/eff.root\"\)

# syst jet 25
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet25_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO25,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet25_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet25_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet25_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet25_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet25_IDISO/eff.root\"\)

# syst jet 30
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet30_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO30,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet30_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet30_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet30_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet30_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet30_IDISO/eff.root\"\)

# default syst jet 35
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet35_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO35,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet35_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet35_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet35_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet35_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet35_IDISO/eff.root\"\)

# syst jet 40
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet40_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO40,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet40_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet40_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet40_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet40_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet40_IDISO/eff.root\"\)

# syst jet 45
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet45_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO45,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet45_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet45_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet45_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet45_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet45_IDISO/eff.root\"\)

# syst jet 50
root -b -q plotEff.C+\(\"ElectronFR.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_FRJet50_IDISO\",\"png\",1,0,0,$ElectronFakeRate,$ElectronFO50,$ElectronIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_FRJet50_IDISO/extra\",\"null\",\"HWW_Electron${OUTNAME}_FRJet50_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_FRJet50_IDISO/extra\",\"HWW_Electron${OUTNAME}_FRJet50_IDISO/eff.root\",\"HWW_Electron${OUTNAME}_FRJet50_IDISO/eff.root\"\)

#
# muons
#

MuonFO0=120
MuonFO5=121
MuonFO10=122
MuonFO15=123
MuonFO20=124
MuonFO25=125
MuonFO30=126
MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonIsoNumerator2011=200
MuonIDNumerator2011=201
MuonIDIsoNumerator2011=202

# syst jet 5
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet5_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO5,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet5_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet5_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet5_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet5_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet5_IDISO/eff.root\"\)

# syst jet 10
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet10_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO10,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet10_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet10_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet10_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet10_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet10_IDISO/eff.root\"\)

# default syst jet 15
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet15_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO15,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet15_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet15_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet15_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet15_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet15_IDISO/eff.root\"\)

# syst jet 20
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet20_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO20,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet20_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet20_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet20_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet20_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet20_IDISO/eff.root\"\)

# syst jet 25
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet25_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO25,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet25_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet25_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet25_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet25_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet25_IDISO/eff.root\"\)

# syst jet 30
root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_FRJet30_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO30,$MuonIDIsoNumerator2012\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_FRJet30_IDISO/extra\",\"null\",\"HWW_Muon${OUTNAME}_FRJet30_IDISO/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_FRJet30_IDISO/extra\",\"HWW_Muon${OUTNAME}_FRJet30_IDISO/eff.root\",\"HWW_Muon${OUTNAME}_FRJet30_IDISO/eff.root\"\)

#rm *.so *.d

