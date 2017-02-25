procedure xgmosArcProc(list, rawpath,superbias, superflat)

string  list          {prompt="e.g. @allproc_R150.lst"}
string  rawpath       {prompt="e.g. /Users/jhennawi/extern/Gemini/"}
string  superbias     {prompt="e.g. N20040423_bias"}
string  superflat     {prompt="e.g. N20040423_flat"}

begin

gemini
gmos

# set up the logfile for this reduction
#gmos.logfile="gmosBias.log"

# Make the bias image with gbias
gsreduce(list, fl_over=yes, fl_trim=yes, fl_bias=no, fl_gscr=no, fl_flat=no, fl_gmos=yes, fl_fixp=no, fl_gsap=no, bias=superbias, flat=superflat,geointer="nearest", rawpath=rawpath, fl_cut=no)

end  
