source ~/.cshrc
source ~/.local_cshrc
date
echo Update the webpage via cron job
set AUTOM_DIR="/b/Blazars/superlotis/"
cd $AUTOM_DIR

#Set x-windows server to run in background
#/usr/X11R6/bin/Xvfb :1 -screen 0 1280x1024x24 -ac -terminate &
#setenv DISPLAY :1.0                          #Set X-windows display

idl $XIDL_DIR/IMG/automated/Reduction/pro/automated_updateweb

echo Done with webpage update
date

#killall -9 Xvfb                              #Kill X-windows server
