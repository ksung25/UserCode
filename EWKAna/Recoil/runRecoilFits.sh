#! /bin/bash

INPUTDIR=/data/blue/ksung/EWKAna/8TeV/Selection

root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"ZmmData\"\)
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",2,2,1,\"ZmmMC\"\)
#root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/wm_select.root\",2,2,1,\"WmMC\"\)

root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",2,2,1,\"ZeeData\"\)
#root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",2,2,1,\"ZeeMC\"\)
#root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Wenu/ntuples/we_select.root\",2,2,1,\"WeMC\"\)

rm *.so *.d
