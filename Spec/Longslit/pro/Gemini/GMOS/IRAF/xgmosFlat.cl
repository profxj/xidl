procedure xgmosFlat(list,rawpath,superflat,superbias,order)

string  list          {prompt="e.g. @gmos_R150_flats.lst"}
string  rawpath       {prompt="e.g. /Users/jhennawi/RAW_DATA/GEMINI_apr23_2003/"}
string  superflat     {prompt="e.g. N20040423_flat400"}
string  superbias     {prompt="e.g. N20040423_bias"}
# 13-Apr-2011 GW: change order type from int to string to allow for CCD-specific entries, e.g. "15,15,25" for CCDs 1-3 in that order
# 13-Apr-2011 GW: increase default from 7 to 15 to reduce artifacts
string	order = 15    {prompt="order of polynomial for flat field"}

begin

string	combname, tmpfil

#tmpfil = mktemp("tmpoutfil")
#tmpfil = "test1"

gemini
gmos



# set up the logfile for this reduction
#gmos.logfile="gmosFlat.log"

# Make the flat with gsflat
gsflat(list,superflat,order=order,rawpath=rawpath,bias=superbias, fl_fixpix=no, fl_dete=yes)

###################################
#  These lines were a mistake...  Saving to remind me of this!
#Mosaic
#gmosaic(tmpfil, outimag=superflat, fl_past=no, fl_vard=no, fl_clea=yes, fl_fixp=no, geointe="nearest", gap=37)
#imdelete(tmpfil, verif-)

end  
