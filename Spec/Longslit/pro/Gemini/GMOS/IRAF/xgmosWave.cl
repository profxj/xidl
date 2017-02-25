procedure gmosWave(rawpath,arc,superbias)

string  rawpath       {prompt="e.g., /Users/jhennawi/RAW_DATA/GEMINI_apr23_2003/"}
string  arc	      {prompt="e.g., N20040423_arc400"}
string  superbias     {prompt="e.g., N20040423_bias"}
begin

string	gsarc
gsarc = "gs" // arc


gemini
gmos

# set up the logfile for this reduction
#gmos.logfile="gmosWave.log"

# Reduce the CuAr arc.  The CuAr is not flat fielded, 
# gaps are not fixpixed.
gsreduce(arc,rawpath=rawpath,fl_flat=no,fl_fixpix=no,bias=superbias)

# For CuAr taken through the 2.0 arcsec slit must smooth the data so that
# accurate line centers can be found by GSWAVELENGTH.  Also need to adjust
# the line width, centering radius and minimum separation parameters.

#changed interactive feature to no
# Establish the wavelength calibration
gswavelength(gsarc,fl_inter=no,verbose=yes)

# Transform the CuAr spectrum, for checking that the transformation is OK.
# Output image name is defined by the default outpref="t"
gstransform(gsarc,wavtran=gsarc) 

end  
