#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/data/smurf/dlevans/LeptonTree/forKevin

#
# Cut-and-count
#
root -l -q plotEff.C+\(\"example.bins\",0,0,0,0,\"${NTUPLEDIR}/leptonTree.root\",\"count\",\"png\",1,0,0\)

#
# Fit
#
root -l -q plotEff.C+\(\"example.bins\",1,0,1,1,\"${NTUPLEDIR}/leptonTree.root\",\"fit\",\"png\",1,0,0\)

rm *.so *.d
