#! /bin/bash

HFLUMI=0.018729
PIXLUMI=0.018479
LUMI=${HFLUMI}

root -l -q fitWm.C+\(\"May23/Wmunu\",${LUMI},0\)
#root -l -q fitZmm.C+\(\"May23/Zmumu\",${LUMI},0\)

root -l -q fitWe.C+\(\"May23/Wenu\",${LUMI},0\)
#root -l -q fitZee.C+\(\"May23/Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

root -l -q plotZmm.C+\(\"Zplots/Zmumu\",${LUMI}\)
root -l -q plotZee.C+\(\"Zplots/Zee\",${LUMI}\)

rm *.so *.d
