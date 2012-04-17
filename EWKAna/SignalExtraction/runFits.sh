#! /bin/bash

LUMI=0.049
ECM=7

#root -l -q fitWe.C+\(\"test\",${LUMI},${ECM},1\)
#root -l -q fitWm.C+\(\"test\",${LUMI},${ECM},1\)

rm *.so *.d
