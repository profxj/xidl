procedure xgmosBias(list,rawpath,superbias)

string  list          {prompt="e.g. @gmosBias.lst"}
string  rawpath       {prompt="e.g. /Users/jhennawi/extern/Gemini/"}
string  superbias     {prompt="e.g. N20040423_bias"}

begin

gemini
gmos

# set up the logfile for this reduction
#gmos.logfile="gmosBias.log"

# Make the bias image with gbias
gbias(list,superbias,rawpath=rawpath)

end  
