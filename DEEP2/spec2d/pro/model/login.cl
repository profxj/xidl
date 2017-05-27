# LOGIN.CL -- User login file for the IRAF command language.

# Identify login.cl version (checked in images.cl).
if (defpar ("logver"))
    logver = "IRAF V2.11 May 1997"

set	home		= "/coma0/jnewman/deimos/model/"
set	imdir		= "HDR$pix/"
set	uparm		= "home$uparm/"
set	userid		= "jnewman"
set     raw             = "/cdrom/keckk34_990309_raw/"
set     idata           = "/dazed/leonard/Keck/Practice/Morepractice/"
set 	model		= "/coma0/jnewman/deimos/model/"


# Set the terminal type.
#if (envget("TERM") == "sun") {
#    if (!access (".hushiraf"))
#	print "setting terminal type to gterm..."
#    stty gterm
#} else {
#    if (!access (".hushiraf"))
#	print "setting terminal type to xterm..."
#    stty xterm
#}

stty xgterm

# Uncomment and edit to change the defaults.
set	editor		= emacs
set	printer		= lw
set	stdimage	= imt2048
set	stdimcur	= stdimage
set	stdplot		= lw
#set	clobber		= no
#set	filewait	= yes
#set	cmbuflen	= 512000
#set	min_lenuserarea	= 64000
#set	imtype		= "imh"

# IMTOOL/XIMAGE stuff.  Set node to the name of your workstation to
# enable remote image display.  The trailing "!" is required.
#set	node		= "my_workstation!"

# CL parameters you mighth want to change.
#ehinit   = "nostandout eol noverify"
#epinit   = "standout showall"
showtype = yes

# Default USER package; extend or modify as you wish.  Note that this can
# be used to call FORTRAN programs from IRAF.

package user

task	$adb $bc $cal $cat $comm $cp $csh $date $dbx $df $diff	= "$foreign"
task	$du $find $finger $ftp $grep $lpq $lprm $ls $mail $make	= "$foreign"
task	$man $mon $mv $nm $od $ps $rcp $rlogin $rsh $ruptime	= "$foreign"
task	$rwho $sh $spell $sps $strings $su $telnet $tip $top	= "$foreign"
task	$touch $vi $emacs $w $wc $less $rusers $sync $pwd $gdb	= "$foreign"

task	$xc $mkpkg $generic $rtar $wtar $buglog			= "$foreign"
#task	$fc = "$xc -h $* -limfort -lsys -lvops -los"
task	$fc = ("$" // envget("iraf") // "unix/hlib/fc.csh" //
	    " -h $* -limfort -lsys -lvops -los")
task	$nbugs = ("$(setenv EDITOR 'buglog -e';" //
	    "less -Cqm +G " // envget ("iraf") // "local/bugs.*)")
task	$cls = "$clear;ls"

# Local scripts that should be loaded at startup

task    iterstat = home$scripts/iterstat.cl
task    minv = home$scripts/minv.cl
task    szap = home$scripts/szap.cl
task    fileroot = home$scripts/fileroot.cl
task    makemask = home$scripts/makemask.cl
task    qzap = home$scripts/qzap.cl
task    xzap = home$scripts/xzap.cl
task    fgain = home$scripts/fgain.cl
task    nickelred = home$scripts/nickelred.cl
task    airmass = home$scripts/airmass.cl
task    tempsub = home$scripts/tempsub.cl
task    fitgauss = home$scripts/fitgauss.cl
task    nickelbase = home$scripts/nickelbase.cl
task    nickelscan = home$scripts/nickelscan.cl
task    backsub = home$scripts/backsub.cl
task    backsub1 = home$scripts/backsub1.cl
task    setup = home$scripts/setup.cl
task    crremove = home$scripts/crremove.cl
task    flatcor = home$scripts/flatcor.cl
task    meccdproc = home$scripts/meccdproc.cl
task    mehedit   = home$scripts/mehedit.cl
task    meapall = home$scripts/meapall.cl
task    crremove_old = home$scripts/crremove_old.cl
task    meapall2 = home$scripts/meapall2.cl
task    arcfits = home$scripts/arcfits.cl
task    objectfits = home$scripts/objectfits.cl
task    standardfits = home$scripts/standardfits.cl
task    medoitall = home$scripts/medoitall.cl
task    meflux  = home$scripts/meflux.cl
task    meadd   = home$scripts/meadd.cl
task    medispcor = home$scripts/medispcor.cl
task    meadd3  = home$scripts/meadd3.cl
task    merfits = home$scripts/merfits.cl
task    medoitallagain = home$scripts/medoitallagain.cl
task    meadd4 = home$scripts/meadd4.cl
task    mechangenames = home$scripts/mechangenames.cl
task    twoampbias = home$scripts/twoampbias.cl
task    mebackground = home$scripts/mebackground.cl
task    meflat2 = home$scripts/meflat2.cl
task    meflat3 = home$scripts/meflat3.cl
task    mewspectext = home$scripts/mewspectext.cl
task    meimcopy = home$scripts/meimcopy.cl
task    meapallarc = home$scripts/meapallarc.cl
task    medoitallarc = home$scripts/medoitallarc.cl
task    medoitallstandard = home$scripts/medoitallstandard.cl
task    lrisfixhead = home$scripts/lrisfixhead.cl
task 	deimos = model$deimos.cl

if (access ("home$loginuser.cl"))
    cl < "home$loginuser.cl"
;

keep;   clpackage

prcache directory
cache   directory page type help

# Print the message of the day.
if (access (".hushiraf"))
    menus = no
else {
    clear; type hlib$motd
}

# Delete any old MTIO lock (magtape position) files.
if (deftask ("mtclean"))
    mtclean
else
    delete uparm$mt?.lok,uparm$*.wcs verify-

# List any packages you want loaded at login time, ONE PER LINE.
images          # general image operators
plot            # graphics tasks
dataio          # data conversions, import export
lists           # list processing
#deimos

# The if(deftask...) is needed for V2.9 compatibility.
if (deftask ("proto"))
    proto       # prototype or ad hoc tasks

tv              # image display
utilities       # miscellaneous utilities
noao            # optical astronomy packages
imred
ccdred
twodspec
longslit
kpnoslit
onedspec


keep
