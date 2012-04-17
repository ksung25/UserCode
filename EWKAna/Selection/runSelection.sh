#! /bin/bash

NTUPDIR=/data/blue/ksung/EWKAna/8TeV/Selection
LUMI=0.013

root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\"\)
root -l -q plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)

root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\"\)
root -l -q plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)

root -l -q selectWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu\"\)
root -l -q plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

root -l -q selectWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu\"\)
root -l -q plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

rm *.so *.d
