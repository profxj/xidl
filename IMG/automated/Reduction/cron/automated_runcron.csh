source ~/.cshrc
source ~/.local_cshrc
date
echo Beginning data reduction via cron job
#set AUTOM_DIR="/b/Blazars/superlotis/"
#set AUTOM_DIR="/c/Blazars/data/testbed/"
cd $AUTOM_DIR
setenv IDL_STARTUP /u/xavier/idl/xidl/IMG/automated/Reduction/cron/idlReductionStartup

#Set x-windows server to run in background
#/usr/X11R6/bin/Xvfb :1 -screen 0 1280x1024x24 -ac -terminate &
#setenv DISPLAY :1.0                          #Set X-windows display

idl
#idl $XIDL_DIR/IMG/automated/Reduction/pro/automated_cron_wrapper

echo Done with data reduction via cron job
date
setenv IDL_STARTUP /u/xavier/idl/.idlstartup

#killall -9 Xvfb                              #Kill X-windows server
