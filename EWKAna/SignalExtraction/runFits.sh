#! /bin/bash

LUMI=0.018618

ECM=8

root -l -q fitWm.C+\(\"Wmunu\",${LUMI},${ECM},0\)
root -l -q fitZmm.C+\(\"Zmumu\",${LUMI},${ECM},0\)

root -l -q fitWe.C+\(\"Wenu\",${LUMI},${ECM},0\)
root -l -q fitZee.C+\(\"Zee\",${LUMI},${ECM},0\)

rm *.so *.d
