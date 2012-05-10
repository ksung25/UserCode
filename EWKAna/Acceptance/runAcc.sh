#! /bin/bash

#
# W->munu
#
root -l -q computeAccGenWm.C+\(\"wm.conf\",\"Wmunu/plus\",1\)
root -l -q computeAccGenWm.C+\(\"wm.conf\",\"Wmunu/minus\",-1\)
root -l -q computeAccGenWm.C+\(\"wm.conf\",\"Wmunu/incl\",0\)
root -l -q computeAccSelWm.C+\(\"wm.conf\",\"Wmunu/plus\",1,1\)
root -l -q computeAccSelWm.C+\(\"wm.conf\",\"Wmunu/minus\",-1,1\)
root -l -q computeAccSelWm.C+\(\"wm.conf\",\"Wmunu/incl\",0,1\)

#
# W->enu
#
root -l -q computeAccGenWe.C+\(\"we.conf\",\"Wenu/plus\",1\)
root -l -q computeAccGenWe.C+\(\"we.conf\",\"Wenu/minus\",-1\)
root -l -q computeAccGenWe.C+\(\"we.conf\",\"Wenu/incl\",0\)
root -l -q computeAccSCWe.C+\(\"we.conf\",\"Wenu/plus\",1,1\)
root -l -q computeAccSCWe.C+\(\"we.conf\",\"Wenu/minus\",-1,1\)
root -l -q computeAccSCWe.C+\(\"we.conf\",\"Wenu/incl\",0,1\)
root -l -q computeAccSelWe.C+\(\"we.conf\",\"Wenu/plus\",1,1\)
root -l -q computeAccSelWe.C+\(\"we.conf\",\"Wenu/minus\",-1,1\)
root -l -q computeAccSelWe.C+\(\"we.conf\",\"Wenu/incl\",0,1\)

#
# Z->mumu
#
root -l -q computeAccGenZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
root -l -q computeAccSelZmm.C+\(\"zmm.conf\",\"Zmumu\",1\)

#
# Z->ee
#
root -l -q computeAccGenZee.C+\(\"zee.conf\",\"Zee\"\)
root -l -q computeAccSCZee.C+\(\"zee.conf\",\"Zee\",1\)
root -l -q computeAccSelZee.C+\(\"zee.conf\",\"Zee\",1\)

rm *.so *.d
