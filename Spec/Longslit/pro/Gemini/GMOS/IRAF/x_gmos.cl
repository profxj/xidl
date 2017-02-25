# Package script for MY_GMOS package

# load necessary packages

#set   x_gmos        = "./x_gmos/"
package x_gmos

# Put task definitions here


task    xgmosBias       = "x_gmos$xgmosBias.cl"
task    xgmosFlat       = "x_gmos$xgmosFlat.cl"
task    xgmosProc       = "x_gmos$xgmosProc.cl"
task    xgmosArcProc    = "x_gmos$xgmosArcProc.cl"
task	xgsflat	= "x_gmos$xgsflat.cl"
#task	xgbias	= "x_gmos$xgbias.cl"
task	xgsreduce	= "x_gmos$xgsreduce.cl"
task    xgemcombine    = "x_gmos$xgemcombine.cl"


clbye()
