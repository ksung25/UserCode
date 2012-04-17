#! /bin/bash

INPUTDIR=/data/blue/ksung/EWKAna/8TeV/Selection
OUTPUTDIR=/data/blue/ksung/EWKAna/8TeV/Efficiency

#
# Select probes for muon efficiencies
#
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",1,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",2,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff\",3,1\)

root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/data_select.root\",\"${OUTPUTDIR}/Data_MuHLTEff\",0,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/data_select.root\",\"${OUTPUTDIR}/Data_MuSelEff\",1,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/data_select.root\",\"${OUTPUTDIR}/Data_MuTrkEff\",2,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/data_select.root\",\"${OUTPUTDIR}/Data_MuStaEff\",3,0\)


#
# Select probes for electron efficiencies
#
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/zee_select.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",0,1\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/zee_select.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",1,1\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/zee_select.root\",\"${OUTPUTDIR}/Zee_EleRecoEff\",2,1\)

root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/data_select.root\",\"${OUTPUTDIR}/Data_EleHLTEff\",0,0\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/data_select.root\",\"${OUTPUTDIR}/Data_EleSelEff\",1,0\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/data_select.root\",\"${OUTPUTDIR}/Data_EleRecoEff\",2,0\)

rm *.so *.d
