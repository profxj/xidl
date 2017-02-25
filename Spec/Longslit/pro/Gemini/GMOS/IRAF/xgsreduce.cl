# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure gsreduce(inimages) 

# Reduction of spectroscopic data - a wrapper for gireduce
# Difference in spectroscopy mode is that the flatfield generated 
# by "gsflat" will have been processed through gmosaic/gpaste, and thus
# it may be either a single gmosaiced science extension, or three non-mosaiced
# extensions.  So, in gsreduce we need to make the basic 
# calibration steps and then potentially gmosaic before applying the flatfield. 
# When the flatfield is already mosaiced, the flatfielding and error 
# propogation is done within this task, not by gireduce. 
#
# Science frame for full reduction is:
#       sci = (a-b-c-d)/e
# where a = original image, b = overscan fit (using colbias); 
#       c = bias image; d = dark image; e = flatfield image
#
# Variance frame for full reduction is:
#       var = (a+b+c+d+e*e*f)/(g*g)
# where a = var(original image); b = (rms of overscan fit)**2; 
#       c = var(bias); d = var(dark); e = final science frame
#       f = var(flatfield); g = flatfield image
#
# DQ frame for full reduction is the combination of the individual
# data quality frames using addmasks
#
# Version  Feb 28, 2002  ML,BM  v1.3 release
#          Jul 07, 2002  IJ fixed filter error msg bug
#          Aug 15, 2002  BM update for new gscut (x/yoffset parameters)
#          Aug 18, 2002  BM fix writing object titles, use MDFROW to match
#                            extension with MDF row
#          Aug 19, 2002  IJ exit gracefully when bias, dark or flat field 
#                           missing, handle the flags consistently.
#          Aug 27, 2002  IJ new parameter geointer to use with gmosaic
#          Sept 19, 2002 IJ flag for gscrrej
#          Sept 20, 2002    v1.4 release
#          Oct 15, 2002  IJ allow open filter position to match any open position
#          Sep 03, 2003  MB allow flatfields with 3 science extensions.
#                        MB added yadd parameter for gscut calls
#	   Nov 25, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#          Feb 06, 2004  KL Fix gain division bug
#          Feb 09, 2004  BM Fix imgets bug when determining nsci on GPREPARED spectra
#	   Feb 09, 2004  KL Moved mktemp call into loops (FITS kernel cache problem)
#          Sep 16, 2004  IJ Support for GSCUT with edge finding
 
char    inimages        {"",prompt="Input GMOS images or list"}
char    outimages       {"",prompt="Output images or list"}
char    outpref         {"gs",prompt="Prefix for output images"}
bool    fl_over         {no,prompt="Subtract overscan level"}
bool    fl_trim         {yes,prompt="Trim off the overscan section"}
bool    fl_bias         {yes,prompt="Subtract bias image"}
bool    fl_gscrrej      {no,prompt="Clean images for cosmic rays"}
bool    fl_dark         {no,prompt="Subtract (scaled) dark image"}
bool    fl_flat         {yes,prompt="Apply flat field correction"}
bool    fl_gmosaic      {yes,prompt="Mosaic science extensions"}
bool    fl_fixpix       {yes,prompt="Interpolate across chip gaps if mosaicing"}
bool    fl_gsappwave    {yes,prompt="Run gsappwave on reduced image"}
bool    fl_cut          {yes,prompt="Cut slits into separate spectra if mosaicing"}
bool    fl_title        {yes,prompt="Put object id in title of cut spectra (MOS only)"}
bool    fl_vardq        {no,prompt="Create variance and data quality frames"}
char    bias            {"",prompt="Bias image name"}
char    dark            {"",prompt="Dark image name"}
char    flatim          {"",prompt="Flatfield (output of GSFLAT) image"}
char    geointer        {"linear",min="linear|nearest|poly3|poly5|spline3|sinc",prompt="Interpolation to use if mosaicing"}
char    gradimage       {"",prompt="Image to use for finding slit edges using the gradient"}
char    refimage        {"",prompt="Reference image for slit positions"}
char    key_exptime     {"EXPTIME",prompt="Header keyword for exposure time"}
char    key_biassec     {"BIASSEC",prompt="Header keyword for overscan strip image section."} 
char    key_datasec     {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}
# gscrrej
bool    fl_inter        {no,prompt="Interactive cosmic ray cleaning fitting"}

#gprepare
char    rawpath         {"",prompt="GPREPARE: Path for input raw images"}
char    sci_ext         {"SCI",prompt="Name of science extension"}
char    var_ext         {"VAR",prompt="Name of variance extension"}
char    dq_ext          {"DQ",prompt="Name of data quality extension"}
char    key_mdf         {"MASKNAME",prompt="Header keyword for the Mask Definition File filename"}
char    mdffile         {"",prompt="MDF file to use if keyword not found"} 
char    mdfdir          {"gmos$data/", prompt="MDF database directory"}
char    bpm             {"",prompt="Bad Pixel Mask filename"}
char    gaindb          {"default",prompt="Database with gain data"}
char    gratingdb       {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
char    filterdb        {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
real    xoffset         {INDEF,prompt="X offset in wavelength [nm]"}
real    yoffset         {INDEF,prompt="Y offset in unbinned pixels"}
real    yadd            {0,prompt="Additional pixels to add to each end of MOS slitlet lengths"}
char    bpmfile         {"gmos$data/chipgaps.dat",prompt="Info on location of chip gaps"}
char    key_ron         {"RDNOISE",prompt="Header keyword for readout noise"}
char    key_gain        {"GAIN",prompt="Header keyword for gain (e-/ADU"}
real    ron             {3.5,prompt="Readout noise value to use if keyword not found"} 
real    gain            {2.2,prompt="Gain value to use if keyword not found"}
int     sat             {65000,prompt="Saturation level in raw images"}

#colbias
bool    ovs_flinter     {no,prompt="Interactive overscan fitting"}
bool    ovs_med         {no,prompt="Use median instead of average in column bias"} 
char    ovs_func        {"chebyshev",min="spline3|legendre|chebyshev|spline1",  prompt="Overscan fitting function"}
int     ovs_order       {1,prompt="Order of overscan fitting function"}
real    ovs_lowr        {3.,prompt="Low sigma rejection factor"}
real    ovs_highr       {3.,prompt="High sigma rejection factor"}
int     ovs_niter       {2,prompt="Number of rejection iterations"}


#general 
char    logfile         {"",prompt="Logfile"}
bool    verbose         {yes,prompt="Verbose output?"}
int     status          {0,prompt="Exit status (0=good)"}
struct  *scanfile1      {"",prompt="Internal use only"}
struct  *scanfile2      {"",prompt="Internal use only"}


begin

    file    temp1
    char    tempim1, tempim2, tmpmdf, tmp4
    char    img, outimg, gratid, filtid1, filtid2, igratid, ifiltid1, ifiltid2
    char    flatorig 
    char    l_geointer
    char    l_inimages, l_outimages, l_outpref, l_logfile, l_key_biassec
    char    l_gradimage, l_refimage
    char    l_key_datasec
    char    l_bias, l_dark, l_flatim, junk, tmpout, tmpvar
    char    l_sci_ext, l_var_ext, l_dq_ext, l_filtname, l_key_exptime
    char    l_bpm, l_key_ron, l_key_gain, l_ovs_func, l_kw_filter, l_filter1
    char    l_key_mdf, l_mdffile, l_mdfdir, inlist, outlist, inlistf, outlistf 
    char    dqexp, varexp, suf, mask1, mask2, mask3, tmpim
    char    l_gratingdb, l_filterdb, l_gaindb, l_rawpath, l_bpmfile

    bool    l_fl_over, l_fl_bias, l_fl_dark, l_fl_flat, l_fl_vardq
    bool    l_ovs_flinter
    bool    l_fl_gscrrej, l_ovs_med
    bool    l_verbose, useprefix, l_fl_gmosaic, mosmode, l_fl_trim, l_fl_inter
    bool    l_fl_fixpix, l_fl_title, l_fl_gsappwave
    bool    fl_rawpath, l_fl_cut
    bool    l_giflat

    int     l_sat, l_ovs_order, l_ovs_niter, nbad, grtilt, igrtilt
    int     nim, nout, obsmode, nsci, xbin, ybin #,l_xbinout,l_ybinout
    int     nx, ny, n_i, nobjects, nsciflat
    char    objectid, gaindbname

    real    l_ron, l_gain, l_ovs_lowr, l_ovs_highr, l_xoffset, l_yoffset, l_yadd
    struct  l_struct

    #Keep some task parameters from changing from the outside
    cache ("imgets", "keypar", "gimverify", "fparse", "tinfo", "gemdate")

    #Will use gireduce for trimming, overscan, and bias subtraction.  In the
    #case of a 3 science extension (i.e. non-gmosaiced) flatfield, we will
    #also use gireduce to do the flatfielding.
    #Before calling gireduce need to check that the flatfield given was really 
    #produced by gsflat, determine how many science extensions it has, and that 
    #the grating ID and filter ID's match the input images.
    #
    # Set values of local variables 

    l_inimages=inimages; l_outimages=outimages; l_outpref=outpref
    l_logfile=logfile
    l_fl_over=fl_over; l_fl_bias=fl_bias; l_fl_dark=fl_dark; l_fl_flat=fl_flat
    l_fl_gscrrej=fl_gscrrej ; l_fl_inter=fl_inter
    l_geointer=geointer
    l_bias=bias; l_dark=dark ; l_flatim=flatim; l_key_exptime=key_exptime  
    l_key_biassec=key_biassec; l_key_datasec=key_datasec; l_fl_vardq=fl_vardq 
    l_sci_ext=sci_ext; l_var_ext=var_ext; l_dq_ext=dq_ext
    l_key_mdf=key_mdf; l_mdffile=mdffile; l_mdfdir=mdfdir ; l_bpm=bpm
    l_key_ron=key_ron; l_key_gain=key_gain ; l_ron=ron ; l_gain=gain  
    l_sat=sat ; l_ovs_flinter=ovs_flinter; l_ovs_med=ovs_med
    l_refimage=refimage ; l_gradimage=gradimage
    l_ovs_func=ovs_func 
    l_ovs_order=ovs_order; l_ovs_lowr=ovs_lowr; l_ovs_highr=ovs_highr; 
    l_ovs_niter=ovs_niter; ; l_verbose=verbose ; l_fl_gmosaic=fl_gmosaic 
    l_gratingdb=gratingdb ; l_filterdb=filterdb ; l_gaindb=gaindb
    l_fl_trim=fl_trim ; l_rawpath=rawpath ; l_fl_cut=fl_cut
    l_fl_fixpix=fl_fixpix ; l_fl_gsappwave=fl_gsappwave ; l_fl_title=fl_title
    l_bpmfile=bpmfile
    l_xoffset=xoffset ; l_yoffset=yoffset ; l_yadd=yadd
    #l_xbinout=int(xbinout) ; l_ybinout=int(ybinout)
    l_giflat=no

    status = 0
    inlist = ""
    outlist = ""
    temp1 = ""

    # Tests the logfile:
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    } else if (l_logfile=="" || l_logfile==" ") {
        l_logfile = gmos.logfile
        if (l_logfile=="" || l_logfile==" ") {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSREDUCE: Both gsreduce.logfile and \
                gmos.logfile are empty", l_logfile, yes)
            printlog ("                    Using default; gmos.log",
                l_logfile, yes)
        }
    }
#    l_logfile=""

    date | scan(l_struct) 
    printlog ("", logfile=l_logfile, verbose=l_verbose)
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)
    printlog ("GSREDUCE -- "//l_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Logs the relevant parameters:
    printlog ("Input image or list  = "//l_inimages, l_logfile, l_verbose)
    printlog ("Output image or list = "//l_outimages, l_logfile, l_verbose)
    printlog ("Output prefix        = "//l_outpref, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("Overscan subtraction  fl_over      = "//l_fl_over,
        l_logfile, l_verbose)
    printlog ("Trim image            fl_trim      = "//l_fl_trim,
        l_logfile, l_verbose)
    printlog ("Bias subtraction      fl_bias      = "//l_fl_bias,
        l_logfile, l_verbose)
    printlog ("Cosmic ray cleaning   fl_gscrrej   = "//l_fl_gscrrej,
        l_logfile, l_verbose)
    printlog ("Dark subtraction      fl_dark      = "//l_fl_dark,
        l_logfile, l_verbose)
    printlog ("Flat field correction fl_flat      = "//l_fl_flat,
        l_logfile, l_verbose)
    printlog ("Mosaic the CCDs       fl_gmosaic   = "//l_fl_gmosaic,
        l_logfile, l_verbose)
    printlog ("Fixpix chip gaps      fl_fixpix    = "//l_fl_fixpix,
        l_logfile, l_verbose)
    printlog ("Cut MOS spectra       fl_cut       = "//l_fl_cut,
        l_logfile, l_verbose)
    printlog ("Add titles to MOS     fl_title     = "//l_fl_title,
        l_logfile, l_verbose)
    printlog ("Run gsappwave         fl_gsappwave = "//l_fl_gsappwave,
        l_logfile, l_verbose)
    printlog ("VAR & DQ planes       fl_vardq     = "//l_fl_vardq,
        l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("geointer = "//l_geointer, l_logfile, l_verbose)
    printlog ("gradimage= "//l_gradimage, l_logfile, l_verbose)
    printlog ("refimage = "//l_refimage, l_logfile, l_verbose)
    printlog ("bias     = "//l_bias, l_logfile, l_verbose)
    printlog ("dark     = "//l_dark, l_logfile, l_verbose)
    printlog ("flatim   = "//l_flatim, l_logfile, l_verbose)
    printlog ("sci_ext  = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("var_ext  = "//l_var_ext, l_logfile, l_verbose)
    printlog ("dq_ext   = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("key_mdf  = "//l_key_mdf, l_logfile, l_verbose)
    printlog ("mdffile  = "//l_mdffile, l_logfile, l_verbose)
    printlog ("mdfdir   = "//l_mdfdir, l_logfile, l_verbose)
    printlog ("bpm      = "//l_bpm, l_logfile, l_verbose)
    printlog ("Chip gaps= "//l_bpmfile, l_logfile, l_verbose)
    printlog ("", logfile=l_logfile, verbose=l_verbose)


    # Check that bias,dark and flat images are defined 
    nbad = 0
    if (l_fl_bias && l_bias=="") {
        printlog ("ERROR - GSREDUCE: No bias image supplied",
            l_logfile, l_verbose)
        nbad+=1
    }
    if (l_fl_dark && l_dark=="") {
        printlog ("ERROR - GSREDUCE: No dark image supplied",
            l_logfile, l_verbose)
        nbad+=1
    }
    if (l_fl_flat && l_flatim=="") {
        printlog ("ERROR - GSREDUCE: No flat field image supplied",
            l_logfile, l_verbose)
        nbad+=1
    }
    if (nbad > 0) {
        printlog ("ERROR - GSREDUCE: Exiting with "//nbad//" errors",
            l_logfile, yes)
        goto error
    }

    nbad = 0

    # check existence of list files
    if (substr(l_inimages,1,1) == "@") {
        inlistf = substr (l_inimages, 2, strlen(l_inimages))
        if (!access(inlistf)) {
            printlog ("ERROR - GSREDUCE: Input list "//inlistf//" not found",
                l_logfile, yes)
            nbad = nbad+1
        }
    }

    if (substr(l_outimages,1,1)=="@") {
        useprefix = no 
        outlistf = substr(l_outimages, 2, strlen(l_outimages))
        if (!access(outlistf)) {
            printlog ("ERROR - GSREDUCE: Output list "//outlistf//" not found",
                l_logfile, yes)
            nbad = nbad+1
        }
    }

    if ((l_outimages =="" || l_outimages ==" ") && (l_outpref =="" || l_outpref ==" ")) {
        printlog ("ERROR - GSREDUCE: Neither outimages nor outpref is \
            specified.", l_logfile, yes)
        nbad = nbad+1
    } else { 
        if (l_outimages =="" || l_outimages ==" ")
            useprefix = yes
        else
            useprefix = no
    }

    #If VAR propogation is requested, make sure the extension names are given 
    if (l_fl_vardq) {
        if (l_dq_ext=="" || l_dq_ext ==" ") {
            printlog ("ERROR - GSREDUCE: extension name dq_ext is missing",
                l_logfile, yes)
            nbad = nbad+1
        } else if (l_var_ext=="" || l_var_ext ==" ") {
            printlog ("ERROR - GSREDUCE: extension name var_ext is missing",
                l_logfile, yes)
            nbad = nbad+1
        }
    }

    #check no commas in sci_ext, var_ext and dq_ext
    if (stridx(",",l_sci_ext)>0 || stridx(",",l_var_ext)>0 || stridx(",",l_dq_ext)>0 ) {
        printlog ("ERROR - GSREDUCE: sci_ext, var_ext or dq_ext contains \
            commas, give root name only", l_logfile, yes)
        nbad = nbad+1
    }

    #check on existence of grating and filter databases
    if (!access(l_gratingdb)) { 
        printlog ("ERROR - GSREDUCE: Gratings database file : "//l_gratingdb//\
            " does not exist", l_logfile, yes)
        nbad = nbad+1
    }
    if (!access(l_filterdb)) { 
        printlog ("ERROR - GSREDUCE: Filters database file : "//l_filterdb//\
            " does not exist", l_logfile, yes)
        nbad = nbad+1
    }
    gaindbname = ""
    if ((l_gaindb != "") && (l_gaindb != "default")) {
        gaindbname = osfn(l_gaindb)
        if (access(gaindbname) == no) {
            printlog ("ERROR - GSREDUCE: Gain database file : "//l_gaindb//\
                " does not exist", l_logfile, yes)
            nbad = nbad+1
        }
    }    

    inlist = mktemp("tmpinlist")
    files (l_inimages, > inlist)
    scanfile1 = inlist
    count (inlist) | scan (nim)
    outlist = mktemp("tmpoutlist")

    #Loop over input images and make sure they are MEFs and also that they are 
    #spectroscopy mode frames

    while (fscan(scanfile1,img) != EOF) {
        fl_rawpath = no
        gimverify (img)
        img = gimverify.outname//".fits"
        if (gimverify.status>=1 && l_rawpath!="") {
            gimverify (l_rawpath//img)
            img = gimverify.outname//".fits"
            fl_rawpath = yes
        }
        if (gimverify.status != 0) { 
            if (gimverify.status == 1) { 
                printlog ("ERROR - GSREDUCE: Input image : "//img//" does \
                    not exist", l_logfile, yes)
                nbad = nbad+1 
            } else if (gimverify.status == 2) { 
                printlog ("ERROR - GSREDUCE: Input image : "//img//" is in \
                    OIF format, not MEF", l_logfile, yes)
                nbad = nbad+1
            } else if (gimverify.status == 3) { 
                printlog ("ERROR - GSREDUCE: Input image : "//img//" is a \
                    GEIS file, not MEF", l_logfile, yes)
                nbad = nbad+1
            } else if (gimverify.status == 4) { 
                printlog ("ERROR - GSREDUCE: Input image : "//img//" is a \
                    simple FITS file, not MEF", l_logfile, yes)
                nbad = nbad+1
            }
        }

        if (nbad == 0) { 
            imgets (img//"[0]", "MASKTYP", >& "dev$null")
            obsmode = int(imgets.value)
            if (obsmode == 0 || obsmode == -1) { 
                printlog ("ERROR - GSREDUCE: "//img//" has MASKTYP other than \
                    spectroscopy mode", l_logfile, yes)
                printlog ("                  This task is not meant for image \
                    or IFU mode images", l_logfile, yes)
                nbad = nbad+1
            } 

            if (useprefix) {
                fparse (img, verbose-)
                files (l_outpref//fparse.root//".fits", >> outlist)  
                scanfile2 = outlist 
            }
            imgets (img//"[0]", "GSREDUCE", >& "dev$null")
            if (imgets.value != " " && imgets.value != "" && imgets.value != "0") {
                printlog ("WARNING - GSREDUCE: "//img//" has already been \
                    processed with GSREDUCE", l_logfile, yes)
                printlog ("                    Make sure this is what you want \
                    to do", l_logfile, yes)
            } 
            imgets (img//"[0]", "GPREPARE", >& "dev$null")
            if ((imgets.value == " " || imgets.value == "0") && (!l_fl_over && !l_fl_trim && !l_fl_bias && !l_fl_dark)){
                printlog ("ERROR - GSREDUCE: Image "//img//" has not been \
                    processed with GPREPARE and no", l_logfile, yes)
                printlog ("                  basic reductions are requested. \
                    Run GPREPARE first", l_logfile, yes)
                nbad = nbad+1
            }
        } #end of if(nbad==0) loop
    } #end of while loop over input images 
    scanfile1 = ""

    if (nbad > 0) {
        printlog ("ERROR - GSREDUCE: Exiting with "//nbad//" errors",
            l_logfile, yes)
        goto error
    }

    if (!useprefix) { 
        files (l_outimages, >> outlist)
        scanfile2 = outlist
        count (outlist) | scan (nout)
        if (nim != nout) {
            printlog ("ERROR - GSREDUCE: Number of input and out output \
                files do not match", l_logfile, yes)
            nbad = nbad+1
        }
    }

    while (fscan(scanfile2,img) != EOF) {
        gimverify (img)
        img = gimverify.outname
        if (gimverify.status != 1) { 
            printlog ("ERROR - GSREDUCE: Output image "//img//" exists",
                l_logfile, yes)
            nbad = nbad+1
        }
    }
    scanfile2 = ""


    # If flatfielding is requested, check that flatfield exists and that it
    # is of the right type and consistent with the grating and filter settings
    # of the input images 
    if (l_fl_flat) {
        imgets (l_flatim//"[0]", "GSFLAT", >& "dev$null")
        flatorig = (imgets.value)
        if (flatorig == "" || flatorig == "0") {
            printlog ("ERROR - GSREDUCE: Flatfield was not produced from \
                GSFLAT", l_logfile, yes)
            nbad = nbad+1
        }
        imgets (l_flatim//"[0]", "NSCIEXT", >& "dev$null")
        nsciflat = int(imgets.value)
        imgets (l_flatim//"[0]", "FILTER1", >& "dev$null")
        filtid1 = (imgets.value)
        imgets (l_flatim//"[0]", "FILTER2", >& "dev$null")
        filtid2 = (imgets.value) 
        imgets (l_flatim//"[0]", "GRATING", >& "dev$null")
        gratid = (imgets.value)
        imgets (l_flatim//"[0]", "GRTILT", >& "dev$null")
        grtilt = real(imgets.value)
        printlog (" ", l_logfile, l_verbose) 
        printlog ("Filter(s)/Grating ID for flatfield: ", l_logfile, l_verbose)
        printlog ("     "//l_flatim//" = "//filtid1//" + "//\
            filtid2//" / "//gratid, l_logfile, l_verbose)
        if (!l_fl_gmosaic && nsciflat==1) {
            printlog ("WARNING - GSREDUCE: fl_flat=yes with a gmosaiced \
                flatfield requires fl_gmosaic=yes", l_logfile, l_verbose)
            printlog ("Setting fl_gmosaic=yes", l_logfile, l_verbose)
            l_fl_gmosaic = yes
        }
        # Changed this so open is open no matter the open position, IJ
        if (substr(filtid1,1,4)=="open") filtid1 = "open"
        if (substr(filtid2,1,4)=="open") filtid2 = "open"

        printlog (" ", l_logfile, l_verbose)

        scanfile1 = inlist
        while (fscan(scanfile1, img) != EOF) {
            if (fl_rawpath)
                img = l_rawpath//img

            imgets (img//"[0]", "GPREPARE", >& "dev$null")
            if (imgets.value=="0") {
                fxhead (img, format_file="", long_head-, count_lines-,
                    >& "dev$null")
                fxhead (img, format_file="", long_head-, count_lines-) | \
                    count ("STDIN") | scan (nsci)
                fxhead (img, format_file="", long_head-, count_lines-) | \
                    fields ("STDIN", "1", lines=str(nsci), quit_if_miss-,
                    print_file-) | scan (nsci)
            } else {
                imgets (img//"[0]", "NSCIEXT", >& "dev$null")
                nsci = int(imgets.value)
            }
            imgets (img//"[0]", "FILTER1", >& "dev$null")
            ifiltid1 = (imgets.value)
            imgets (img//"[0]", "FILTER2", >& "dev$null")
            ifiltid2 = (imgets.value) 
            imgets (img//"[0]", "GRATING", >& "dev$null")
            igratid = (imgets.value)
            imgets (img//"[0]", "GRTILT", >& "dev$null")
            igrtilt = real(imgets.value)
            # Changed this so open is open no matter the open position, IJ
            if (substr(ifiltid1,1,4)=="open") ifiltid1 = "open"
            if (substr(ifiltid2,1,4)=="open") ifiltid2 = "open"
            if (nsci!=3 && nsci!=1 && nsci!=6) {
                printlog ("ERROR - GSREDUCE: Input image "//img//" does not \
                    have one, three, or six science extensions",
                    l_logfile,l_verbose)
                nbad = nbad+1
            } 
            if (nsciflat>1  && (nsciflat!=nsci)) {
                printlog ("ERROR - GSREDUCE: Input flatfield has not been \
                    gmosaiced, and image "//img//" has a different number of \
                    science extensions.", l_logfile, yes)
                nbad = nbad+1
            } 
            if ((ifiltid1 != filtid1) || (ifiltid2 != filtid2)) {
                printlog ("ERROR - GSREDUCE: Filter IDs for image: "//img//\
                    " do not match flatfield", l_logfile, yes)
                nbad = nbad+1
            }
            if ((igratid != gratid) || (igrtilt != grtilt)) {
                printlog ("ERROR - GSREDUCE: Grating ID or tilt for image: "//\
                    img//" does not match flatfield", l_logfile, yes)
                nbad = nbad+1
            }
        }
        scanfile1 = ""
    }
    # Check the flags are set consistently
    if (l_fl_cut && !l_fl_gmosaic) {
        printlog ("WARNING - GSREDUCE: fl_cut=yes requires fl_gmosaic=yes",
            l_logfile, l_verbose)
        printlog ("Setting fl_gmosaic=yes", l_logfile, l_verbose)
        l_fl_gmosaic = yes
    }
    if (l_fl_fixpix && !l_fl_gmosaic) {
        printlog ("WARNING - GSREDUCE: fl_fixpix=yes requires fl_gmosaic=yes",
            l_logfile, l_verbose)
        printlog ("Setting fl_fixpix=no", l_logfile, l_verbose)
        l_fl_fixpix = no
    }

    #If something is wrong then exit with an error status 
    if (nbad > 0)
        goto error

    #Call gireduce. If flatfielding is requested and input flatfield has 
    #already been gmosaiced, then call gireduce for std processing, 
    #and then gmosaic the output and flatten afterwards. 
    #Otherwise, we will use gireduce to apply the flatfield

    if (l_fl_flat) {
        if (nsciflat==3 || nsciflat==6) {
            ##need to trick gireduce into using a non-mosaiced spectral flat
            hedit (l_flatim//"[0]", "GIFLAT", "produced by gsflat", add+,
                addonly-, delete-, show-, verify-, update+, >& "dev$null")
            l_giflat = yes
        }
    }
    temp1 = mktemp("tmp")
    join (inlist//","//outlist, output=temp1, delim=" ", missing="Missing",
        maxchars=161, shortest=no, verbose=no)
    scanfile1 = temp1
    while (fscan(scanfile1,img,outimg) != EOF) {
        # Create tmp FITS file names used within this loop
        tempim1 = mktemp("tmpimage1")//".fits"
        tempim2 = mktemp("tmpimage2")//".fits"

        suf = substr (outimg, strlen(outimg)-3, strlen(outimg))
        if (suf !="fits")
            outimg = outimg//".fits"
        suf = substr (img, strlen(img)-3, strlen(img))
        if (suf !="fits")
            img = img//".fits"
        print ("Input image = "//img//"; Output image = "//outimg)

        nsci = 0
        if (!fl_rawpath) {
            imgets (img//"[0]", "NSCIEXT", >& "dev$null")
            nsci = int(imgets.value)
        }
        if (nsci != 1) {
        
            if (l_gaindb == "default") {
                if (fl_rawpath) ggdbhelper (l_rawpath//img//"[0]")
                else            ggdbhelper (img//"[0]")
                if (ggdbhelper.status != 0) goto error
                gaindbname = osfn(ggdbhelper.gaindb)
                if (access(gaindbname) == no) {
                    printlog ("ERROR - GSREDUCE: Gain database not found",
                        l_logfile, yes)
                    goto error
                }
            }
        
            if (l_fl_over || l_fl_trim || l_fl_bias || l_fl_dark || l_giflat) { 
                printlog ("GSREDUCE: Calling GIREDUCE to process image",
                    l_logfile, l_verbose)            
                gireduce (img, outpref="", outimages=tempim1, sci_ext=l_sci_ext,
                    var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=l_fl_vardq,
                    fl_over=l_fl_over, key_biassec=l_key_biassec,
                    median=l_ovs_med, fl_inter=l_ovs_flinter,
                    function=l_ovs_func, order=l_ovs_order,
                    low_reject=l_ovs_lowr, high_reject=l_ovs_highr,
                    niterate=l_ovs_niter, fl_trim=l_fl_trim,
                    key_datasec=l_key_datasec, fl_bias=l_fl_bias, bias=l_bias,
                    fl_dark=l_fl_dark, dark=l_dark, fl_flat=l_giflat,
                    key_filter="GRATING", flat1=l_flatim, flat2="", flat3="",
                    flat4="", gp_outpref="g", bpm=l_bpm, key_ron=l_key_ron,
                    key_gain=l_key_gain, fl_addmdf=yes, key_mdf=l_key_mdf,
                    mdffile=l_mdffile, mdfdir=l_mdfdir, gaindb=gaindbname,
                    rawpath=l_rawpath, logfile=l_logfile, verbose=l_verbose)

                #Check that gireduce finished ok
                if (gireduce.status != 0) { 
                    printlog ("ERROR - GSREDUCE: There was an apparent fatal \
                        error with GIREDUCE", l_logfile, yes)
                    nbad = nbad+1 
                    goto error
                } 
            } else
                copy (img, tempim1)

            # Clean for cosmic rays
            if (l_fl_gscrrej) {
                gscrrej (tempim1, "x"//tempim1, verbose=l_verbose,
                    fl_inter=l_fl_inter, logfile=l_logfile)
                if (gscrrej.status != 0)
                    goto error
                delete (tempim1, verify-)

                # KL - Commented out to avoid re-using a FITS file name.
                #      Beside, there is no need to rename the image, we just
                #      need to re-assign the variable tempim1 to the new
                #      'x' file name.
                #
                #rename("x"//tempim1,tempim1,field="all")

                tempim1 = "x"//tempim1
            }

            # Multiply by the gains, but check first, since we could get here 
            # with a previously gsreduced image
            # ****
            #  This needs to be done only if the flat has already been 
            #  gmosaiced.  Otherwise, the image has been flatfielded in 
            #  gireduce and therefore the gain has already been applied. 
            #    (KL - 5 feb 2004)
            if (!l_giflat) {
                imgets (tempim1//"[0]", "GGAIN", >& "dev$null")
                if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
                    printlog ("GSREDUCE: Multiplying by the gains, using GGAIN",
                        l_logfile, l_verbose)
                    ggain (tempim1, gaindb=gaindbname, logfile=l_logfile,
                        key_gain=l_key_gain, key_ron=l_key_ron, gain=l_gain,
                        ron=l_ron, fl_update+, fl_mult+, verbose=no,
                        sci_ext=l_sci_ext, var_ext=l_var_ext)
                }
            }

            if (l_fl_gmosaic) {         
                printlog ("GSREDUCE: Calling GMOSAIC", l_logfile, l_verbose)
#                    fl_paste=no, geointer=l_geointer, gap=37, logfile=l_logfile,
                gmosaic (inimages=tempim1, outimages=tempim2, outpre="",
                    fl_paste=no, geointer=l_geointer, gap=37, 
                    sci_ext=l_sci_ext, fl_fixpix=l_fl_fixpix, bpmfile=l_bpmfile,
                    key_datsec=l_key_datasec, fl_real-, verbose=l_verbose) 

                printlog ("GSREDUCE: returned from GMOSAIC.",
                    l_logfile, l_verbose)
                delete (tempim1, verify-, >& "dev$null")
            } else
                tempim2 = tempim1
        } else {  #nsci == 1
            tempim2 = img
            printlog ("", l_logfile, yes)
            if (l_fl_flat) {
                printlog ("WARNING - GSREDUCE: "//img//" has only one \
                    science extension. Assume it has already",
                    l_logfile, yes)
                printlog ("                    been processed as desired - \
                    only flatfielding will be done.", l_logfile, yes)
            } else {
                printlog ("WARNING - GSREDUCE: "//img//" has only one science \
                    extension. Assume it has already", l_logfile, yes)
                printlog ("                    been processed as desired.",
                    l_logfile,yes)
            }
            printlog ("", l_logfile, yes)
        }

        #Now we actually flatten the image.  Since there is only one SCI 
        #extension we can just use imarith. First, tcopy MDF file to a 
        #temporary file, then create the output image with the MDF file.
        #In the temporary image tempim2, the MDF file is in extension 2

        if (l_fl_flat && !l_giflat) {
            tmpmdf = mktemp("tmpmdf")
            tcopy (tempim2//"[MDF]", tmpmdf//".fits", verbose-)
            wmef (input=tmpmdf//".fits", output=outimg, extnames="MDF",
                phu=tempim2, verbose-, >& "dev$null")
            if (wmef.status != 0) { 
                printlog ("ERROR - GSREDUCE: Problem writing output MEF",
                    l_logfile, yes)
                nbad = nbad+1
                goto error 
            } 
            imdelete (tmpmdf//".fits", verif-, > "dev$null")
            junk = mktemp("tmpjunk")
            junk = junk//".fits"

            printlog ("GSREDUCE: Flatfield image : "//img, l_logfile, l_verbose)

            imarith (tempim2//"["//l_sci_ext//"]","/",
                l_flatim//"["//l_sci_ext//"]", junk, verb-)
            fxinsert (junk, outimg//"[1]", groups="", verbose-, >& "dev$null")
            gemhedit (outimg//"[2]", "EXTNAME", l_sci_ext, "Extension name")
            gemhedit (outimg//"[2]", "EXTVER", 1, "Extension version")
            imdelete (junk, verif-, >& "dev$null")

            if (l_fl_vardq) { 
                varexp = "a*b/c**2"
                tmpvar = mktemp("tmpvarout")
                tmpvar = tmpvar//".fits"
                dqexp = "im1 || im2"
                imexpr (varexp, tmpvar, tempim2//"["//l_var_ext//"]",
                    l_flatim//"["//l_var_ext//"]",
                    l_flatim//"["//l_sci_ext//"]",
                    dims="auto", outtype="auto", refim="auto", bwidth="0",
                    btype="nearest", bpixval=0., rangecheck=yes, verbose=no,
                    exprdb="none")
                fxinsert (tmpvar, outimg//"[2]", groups="", verbose-,
                    >& "dev$null")
                gemhedit (outimg//"[3]", "EXTNAME", l_var_ext, "Extension name")
                gemhedit (outimg//"[3]", "EXTVER", 1, "Extension version")
                imdelete (tmpvar, verif-, >& "dev$null")
                mask1 = mktemp("tmpmask1")
                mask1 = mask1//".fits"
                mask2 = mktemp("tmpmask2")
                mask2 = mask2//".fits"
                mask3 = mktemp("tmpmask3")
                mask3 = mask3//".fits"
                imcopy (tempim2//"["//l_dq_ext//"]", mask1, verb-)
                imcopy (l_flatim//"["//l_dq_ext//"]", mask2, verb-)
                addmasks (mask1//","//mask2, mask3, expr=dqexp, flags="")
                fxinsert (mask3, outimg//"[3]", groups="", verb-, >& "dev$null")
                gemhedit (outimg//"[4]", "EXTNAME", l_dq_ext, "Extension name")
                gemhedit (outimg//"[4]", "EXTVER", 1, "Extension version")
            }

            if (!imaccess(outimg)) {
                printlog ("ERROR - GSREDUCE: Output image not written for \
                    some reason. Check logfile.", l_logfile, yes)
                nbad = nbad+1 
                goto error
            }     

            if (tempim1 != img)
                imdelete (tempim1//","//tempim2, verify=no, >& "dev$null")
            if (l_fl_vardq)
                imdelete (mask1//","//mask2//","//mask3, ver-, >& "dev$null")

            gemhedit (outimg//"[0]", "GSFLATIM", l_flatim,
                "Flatfield Image used by gsreduce")
            #end of l_fl_flat 
        } else {
            imrename (tempim2, outimg, verbose+)
        }

        #cut spectra into individual spectra if this is MOS
        if (l_fl_cut) {
#            imgets (tempim2//"[0]", "OBSMODE", > "dev$null")
            imgets (outimg//"[0]", "OBSMODE", > "dev$null")
            if (imgets.value == "MOS")
                mosmode = yes
            else { 
                mosmode = no
                printlog ("GSREDUCE: Image "//img// " is "//imgets.value//". \
                    Will not cut spectra into separate science extensions",
                    l_logfile, l_verbose)
            }  

            if (mosmode) { 
                printlog ("GSREDUCE: Calling GSCUT to cut spectra into \
                    separate science extensions", l_logfile, l_verbose)
                tmpout = mktemp("tmpoutimage")//".fits"
                imrename (outimg, tmpout, verb-) 
                gscut (inimage=tmpout, outimage=outimg, secfile="",
                    logfile=l_logfile, fl_vardq=l_fl_vardq,
                    gratingdb=l_gratingdb, filterdb=l_filterdb,
                    xoffset=l_xoffset, yoffset=l_yoffset, yadd=l_yadd,
                    verbose=l_verbose, sci_ext=l_sci_ext, var_ext=l_var_ext,
                    dq_ext=l_dq_ext, gradimage=l_gradimage, refimage=l_refimage)
                imdelete (tmpout, verify-, >& "dev$null")
                
				if (gscut.status != 0) {
					printlog ("ERROR - GSREDUCE: GSCUT returned with error", \
					    l_logfile, l_verbose)
					nbad = nbad + 1
				    goto error
				}
                # Add titles, use MDFROW in header to match extension with MDF
                if (l_fl_title) {
                    printlog ("GSREDUCE: Updating titles of MOS spectra",
                        l_logfile, l_verbose)
                    imgets (outimg//"[0]", "NSCIEXT", >& "dev$null")
                    nobjects = int(imgets.value)
                    if (nobjects==0) {
                        printlog ("ERROR - GSREDUCE: NSCIEXT==0, no objects \
                            in gscut output.", l_logfile, yes)
                        goto error
                    }
                    #               tinfo(outimg//"[MDF]",ttout-) 
                    #               nobjects=tinfo.nrows
                    for (n_i=1; n_i<=nobjects; n_i+=1) {
                        imgets (outimg//"["//l_sci_ext//","//str(n_i)//"]",
                            "MDFROW", >& "dev$null")
                        if (imgets.value=="0") {
                            printlog ("ERROR - GSREDUCE: MDFROW not found.",
                                l_logfile, yes)
                            goto error
                        }
                        tprint (outimg//"[MDF]", prparam-, prdata+, showrow-,
                            showhdr-, showunits-, col="ID",
                            rows=int(imgets.value)) | scan(objectid)
                        hedit (outimg//"["//l_sci_ext//","//str(n_i)//"]",
                            "i_title", objectid, add+, addonly-, delete-, show-,
                            verify-, update+)
                    }
                }
            }
        }
        # line 602-603 was moved from after gsappwave

        # ----  gsappwave ----------
        if (l_fl_gsappwave) {
            printlog ("GSREDUCE: Calling GSAPPWAVE", l_logfile, l_verbose)
            gsappwave (outimg, logfile=l_logfile, gratingdb = l_gratingdb,
                filterdb="gmos$data/GMOSfilters.dat", key_dispaxis="DISPAXIS",
                dispaxis=1, verbose=l_verbose)
        }

        #Update the PHU
        gemdate ()
        printlog ("GSREDUCE: Updating PHU and cleaning up",
            l_logfile, l_verbose)
        gemhedit (outimg//"[0]", "GSREDUCE", gemdate.outdate,
            "UT Time stamp for GSREDUCE")
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI")

        #End of While loop over all images 
    }

    if (l_giflat)
        hedit (l_flatim//"[0]", "GIFLAT", "produced by gsflat", add-, addonly-,
            delete+, show-, verify-, update+, >& "dev$null")
    scanfile1 = ""

    goto clean

error:
    status = 1
    print (" ")
    printlog ("ERROR - GSREDUCE: Program execution failed with "//str(nbad)//\
        " errors", l_logfile, yes)
    goto clean

clean: 
    scanfile1 = ""
    scanfile2 = ""
    printlog (" ", l_logfile, l_verbose)
    printlog ("GSREDUCE done", l_logfile, l_verbose)
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile, l_verbose)
    delete (temp1//","//inlist//","//outlist, verify=no, >& "dev$null")

end
