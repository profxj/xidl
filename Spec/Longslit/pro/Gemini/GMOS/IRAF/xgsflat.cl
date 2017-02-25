# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure gsflat (inflats,specflat,fl_answer)

# Create a GMOS normalized spectral flatfield (LONGSLIT or MOS) from GCAL flats
#
# Version    Feb 28, 2002  ML,BM,IJ  v1.3 release
#            Aug 15, 2002  BM  add x/yoffset parameters for new gscut
#            Sept 20, 2002    v1.4 release
#            Oct 14, 2002  IJ don't use "i" inside a script
#	     Aug 20, 2003  MB add yadd parameter for new gscut
#                          allow detector by detector fitting
#            Sep 11, 2003  KL IRAF2.12 - new/modified parameters
#                             hedit - addonly
#                             imcombine - headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#            Sep 14, 2004  BM,IJ edge finding of slits, fl_vardq support, cleaning

string  inflats     {"",prompt="Input flatfields"} 
string  specflat    {"",prompt="Output normalized flat (MEF)"}
bool    fl_slitcorr {no,prompt="Correct output for Illumination/Slit-Function"}
string  slitfunc    {"",prompt="Slit Function (MEF output of gsslitfunc)"} 
bool    fl_keep     {no,prompt="Keep imcombined flat?"} 
string  combflat    {"",prompt="Filename for imcombined flat"}  
 
bool    fl_over     {no,prompt="Subtract overscan level"} 
bool    fl_trim     {yes,prompt="Trim off overscan region"}
bool    fl_bias     {yes,prompt="Subtract bias image"}
bool    fl_dark     {no,prompt="Subtract (scaled) dark image"}
bool    fl_fixpix   {yes,prompt="Interpolate across chip gaps"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames"}
string  bias        {"",prompt="Bias image"}
string  dark        {"",prompt="Dark image"}
string  key_exptime {"EXPTIME",prompt="Exposure time header keyword"}
string  key_biassec {"BIASSEC",prompt="Header keyword for overscan strip image section"} 
string  key_datasec {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}
string  rawpath     {"",prompt="GPREPARE: Path for input raw images"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string  key_mdf     {"MASKNAME",prompt="Header keyword for the MDF"}
string  mdffile     {"",prompt="MDF to use if keyword not found"}
string  mdfdir      {"gmos$data/", prompt="MDF database directory"}
string  bpm         {"",prompt="Name of bad pixel mask file or image"}
string  gaindb      {"default",prompt="Database with gain data"}
string  gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
string  bpmfile     {"gmos$data/chipgaps.dat",prompt="Info on location of chip gaps"}
string  refimage    {"",prompt="Reference image for slit positions"}
real    sat         {65000.,prompt="Saturation level in raw images"}
real    xoffset     {INDEF,prompt="X offset in wavelength [nm]"}
real    yoffset     {INDEF,prompt="Y offset in unbinned pixels"}
real    yadd        {0,prompt="Additional pixels to add to each end of MOS slitlet lengths"}	
bool    fl_usegrad  {no,prompt="Use gradient method to find MOS slits"}
bool    fl_emis     {no,prompt="mask emission lines from lamp (affected pixels set to 1. in output)"} 

#Response/Profile fitting 
bool    fl_inter    {no,prompt="Fit response interactively?"}
bool    fl_answer   {yes,prompt="Continue interactive fitting?"}
bool    fl_detec    {no, prompt="Fit response detector by detector rather than slit by slit?"} 
string  function    {"spline3",min="spline3|legendre|chebyshev|spline1", prompt="Fitting function for response"}
string  order       {"15",prompt="Order of fitting function, minimum value=1"} 
real    low_reject  {3,prompt="Low rejection in sigma of response fit"}
real    high_reject {3,prompt="High rejection in sigma of response fit"}
int     niterate    {2,prompt="Number of rejection iterations in response fit"}

#Gemcombine 
string  combine     {"average",prompt="Combination operation"}
string  reject      {"avsigclip",prompt="Rejection algorithm"}
string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
string  scale       {"mean",prompt="Image scaling"}
string  zero        {"none",prompt="Image zeropoint offset"}
string  weight      {"none",prompt="Image weights"}
string  statsec     {"",prompt="Statistics section"}
real    lthreshold  {INDEF,prompt="Lower threshold"}
real    hthreshold  {INDEF,prompt="Upper threshold"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {0,min=0,prompt="Minimum to keep or maximum to reject"}  
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,prompt="Upper sigma clipping factor"}
string  key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
string  key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron         {3.5,min=0.,prompt="Readout noise rms in electrons"}
real    gain        {2.2,min=0.00001,prompt="Gain in e-/ADU"}
string  snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,prompt="Radius (pixels) for neighbor rejection"}

#Colbias
bool    ovs_flinter {no,prompt="Interactive overscan fitting?"}
bool    ovs_med     {no,prompt="Use median instead of average in column bias?"} 
string  ovs_func    {"chebyshev",min="spline3|legendre|chebyshev|spline1", prompt="Overscan fitting function"}
int     ovs_order   {1,prompt="Order of overscan fitting function"}
real    ovs_lowr    {3.,prompt="Low sigma rejection factor"}
real    ovs_highr   {3.,prompt="High sigma rejection factor"}
int     ovs_niter   {2,prompt="Number of rejection iterations"}

#General 
string  logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"} 


begin

    #local variables of parameters 
    string  l_inflats, l_specflat, l_logfile, l_slitfunc
    string  l_bias, l_dark, l_bpm, l_key_biassec, l_key_datasec, l_ovs_func 
    string  l_combflat, l_bpmfile, l_refimage
    string  l_key_mdf, l_mdffile, l_combine, l_reject, l_scale, l_zero
    string  l_weight, l_statsec, l_key_exptime
    string  l_key_ron, l_key_gain 
    string  l_snoise, l_masktype
    string  l_rfunction, l_sci_ext, l_var_ext, l_dq_ext
    string  l_mdfdir, l_gaindb, l_rawpath, l_gratingdb, l_filterdb
    string  l_order
    real    l_lthreshold, l_hthreshold, l_lsigma, l_hsigma,l_ron,l_gain
    real    l_sigscale, l_maskvalue
    real    l_pclip, l_grow, l_sat, l_ovs_lowr, l_ovs_highr
    real    l_rlowrej, l_rhighrej
    real    l_xoffset, l_yoffset, l_yadd
    int     l_nlow, l_nhigh, l_nkeep, l_rorder, l_rdorder[6], l_ovs_order
    int     l_ovs_niter, l_niterate
    bool    l_slitcorr, l_mclip, l_rinteractive, l_fl_vardq, l_verbose
    bool    l_redkeep
    bool    l_ovs_flinter,l_ovs_med,l_fl_over,l_fl_bias,l_fl_dark,l_fl_trim 
    bool    l_fl_fixpix
    bool    l_fl_answer
    bool    l_fl_detec
    bool    l_fl_emis, l_fl_usegrad
    bool    l_fl_gmosaic

    #other parameters used within this task 

    file    temp1, temp2
    bool    reducing, gpflag, giflag, gsflag, gsfl1, mosflag, mosfl1, florder
    bool    specmode 
    string  filt1, filtn, grat1, gratn, cutfile, l_gradimage
    string  mdf, inlist, filelist, scilist, img, specsec, imgraw
    string  response, suf, combsig
    string  infile[20], headtest
    string  scisec, sciflat, flatfit, inscoo, combflatsciext
    string  gprep, gired, specsecfile, tmpflat
    real    firstexp, l_sgain, l_sron, gaineff, roneff, expdiff
    int     colpos1, colpos2, compos, seclen, x1, y1, nextens, next1
    int     n_i, n_j, n_k, l, n, idum, nbad, len, nslits, nim, n_ccd, n_row
    bool    fl_delgrad
    struct  sdate

    char    l_sample
    int     x11, x12, x21, x22, Xbin, Xmax
    int     norder
    char    keyfound
    int     junk


    #Make parameter assignments  
    l_inflats=inflats ; l_specflat=specflat
    l_combflat=combflat ; l_redkeep=fl_keep
    l_fl_over=fl_over ; l_fl_bias=fl_bias ; l_fl_dark=fl_dark
    l_fl_fixpix=fl_fixpix
    l_bias=bias ; l_dark=dark ; l_bpm=bpm   ; l_bpmfile=bpmfile
    l_key_mdf=key_mdf ; l_mdffile=mdffile 
    l_key_biassec=key_biassec ; l_key_datasec=key_datasec 
    l_ovs_flinter=ovs_flinter ; l_ovs_med=ovs_med 
    l_ovs_func=ovs_func ; l_ovs_order=ovs_order  
    l_slitcorr=fl_slitcorr ; l_slitfunc=slitfunc
    l_combine=combine ; l_reject=reject
    l_scale=scale ; l_zero=zero ; l_weight=weight
    l_statsec=statsec ; l_key_exptime=key_exptime
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma
    l_verbose=verbose ; l_mclip=mclip
    junk = fscan (key_gain, l_key_gain)
    junk = fscan (key_ron, l_key_ron)
    l_gain=gain ; l_ron=ron ; l_snoise=snoise 
    l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow
    l_sci_ext = sci_ext ; l_dq_ext = dq_ext ; l_var_ext=var_ext
    l_ovs_lowr=ovs_lowr ; l_ovs_highr=ovs_highr ; l_ovs_niter=ovs_niter
    l_fl_vardq=fl_vardq
    l_rinteractive=fl_inter ;  l_rfunction=function;  
    l_order=order	
    l_rlowrej=low_reject ; l_rhighrej=high_reject; l_niterate=niterate
    l_logfile=logfile; l_sat=sat
    l_fl_trim=fl_trim
    l_mdfdir=mdfdir ; l_gaindb=gaindb ; l_rawpath=rawpath
    l_gratingdb=gratingdb ; l_filterdb=filterdb
    l_xoffset=xoffset ; l_yoffset=yoffset; l_yadd=yadd
    l_fl_emis=fl_emis
    l_fl_detec=fl_detec
    l_refimage=refimage ; l_fl_usegrad=fl_usegrad
    l_masktype = masktype
    l_maskvalue = maskvalue

    status = 0
    if (l_fl_detec) {
        l_fl_gmosaic = no
        l_fl_fixpix = no
    } else
        l_fl_gmosaic = yes

    #Emission Line Masking not enabled yet
    if (l_fl_emis) {
        printlog ("WARNING - GSFLAT: Emission line masking is not yet \
            implemented", l_logfile, yes)
        printlog ("                  Setting fl_emis = no", l_logfile, yes)
        l_fl_emis = no
    }

    # Keep some parameters from changing by outside world
    cache ("imgets", "hedit", "gmos", "tinfo", "fparse", "gimverify", "gemdate")

    #Define temporary files
    specsecfile = mktemp("tmpsecfile")
    scilist = mktemp("tmpscilist")
    tmpflat = mktemp("tmpcombflat")
    scisec = mktemp("tmpscisec")
    if (!l_redkeep) 
        l_combflat = mktemp("tmpcombflat")
    if (l_fl_vardq)
        combsig = mktemp("tmpsig")
    response = mktemp("tmpresponse")
    temp1 = mktemp("tmpfilelist")
    temp2 = mktemp("tmpfilelist")
    mdf = mktemp("tmpmdffile")//".fits" 

    #Check logfile ... 
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    } else if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSFLAT: both gsflat.logfile and gmos.logfile \
                are empty.", l_logfile, l_verbose)
            printlog ("                  Using default file gmos.log.",
                l_logfile, l_verbose)
        }
    }
    date | scan(sdate)
    printlog ("",l_logfile,l_verbose)
    printlog ("------------------------------------------------------------\
        --------------------", l_logfile, l_verbose)
    printlog ("GSFLAT -- "//sdate,l_logfile,l_verbose)
    printlog ("",l_logfile,l_verbose)
    printlog ("Input images or list              = "//l_inflats,
        l_logfile,l_verbose)
    printlog ("Output spectral flat              = "//l_specflat,
        l_logfile,l_verbose)
    printlog ("Correct for slit function         = "//l_slitcorr,
        l_logfile,l_verbose) 
    if (l_slitcorr) 
        printlog ("Slitfunc image from GSSLITFUNCTION= "//l_slitfunc,
            l_logfile,l_verbose)
    printlog ("Keep combined flat combflat       = "//l_redkeep,
        l_logfile,l_verbose) 
    if (l_redkeep)
        printlog ("Combined flat (not normalized)    = "//l_combflat,
            l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("Fitting parameters for Spectral Flatfield: ",l_logfile,l_verbose)
    printlog("  interactive = "//l_rinteractive,l_logfile,l_verbose)
    printlog("  function    = "//l_rfunction,l_logfile,l_verbose)
    printlog("  order       = "//l_order,l_logfile,l_verbose)
    printlog("  low_reject  = "//l_rlowrej,l_logfile,l_verbose)
    printlog("  high_reject = "//l_rhighrej,l_logfile,l_verbose)
    printlog("  niterate    = "//l_niterate,l_logfile,l_verbose)
    printlog("  Fit detectors separately fl_detec = "//l_fl_detec,
        l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("Overscan Subtraction fl_over   = "//l_fl_over,l_logfile,l_verbose)
    printlog("Trim image           fl_trim   = "//l_fl_trim,l_logfile,l_verbose)
    printlog("Bias Subtraction     fl_bias   = "//l_fl_bias,l_logfile,l_verbose)
    printlog("Dark Subtraction     fl_dark   = "//l_fl_dark,l_logfile,l_verbose)
    printlog("VAR & DQ planes      fl_vardq  = "//l_fl_vardq,l_logfile,l_verbose) 
    printlog("Fixpix chip gaps     fl_fixpix = "//l_fl_fixpix,l_logfile,l_verbose)
    printlog("Fit detec. by detec. fl_detec  = "//l_fl_detec,l_logfile,l_verbose)
    printlog("Use gradient method  fl_usegrad= "//l_fl_usegrad,l_logfile,l_verbose)
    printlog("Mask emission lines  fl_emis   = "//l_fl_emis,l_logfile,l_verbose)
    printlog("Reference image      refimage  = "//l_refimage,l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("bias     = "//l_bias,l_logfile,l_verbose)
    printlog("dark     = "//l_dark,l_logfile,l_verbose)
    printlog("sci_ext  = "//l_sci_ext,l_logfile,l_verbose)
    printlog("var_ext  = "//l_var_ext,l_logfile,l_verbose)
    printlog("dq_ext   = "//l_dq_ext,l_logfile,l_verbose)
    printlog("bpm      = "//l_bpm,l_logfile,l_verbose)
    printlog("key_mdf  = "//l_key_mdf,l_logfile,l_verbose)
    printlog("mdffile  = "//l_mdffile,l_logfile,l_verbose)
    printlog("mdfdir   = "//l_mdfdir,l_logfile,l_verbose)
    printlog("Chip gaps= "//l_bpmfile,l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)

    nbad = 0

    #The checks on things begin ... 

    if (l_inflats == "" || l_inflats==" ") {
        printlog ("ERROR - GSFLAT: Input Flats not specified", l_logfile, yes)
        nbad = nbad+1
    }

    # check existence of list files
    if (substr(l_inflats,1,1) == "@") {
        inlist = substr (l_inflats, 2, strlen(l_inflats))
        if (!access(inlist)) {
            printlog ("ERROR - GSFLAT: Input list "//inlist//" not found",
                l_logfile, verbose+)
            nbad = nbad+1 
        }
    }

    if (l_specflat == "" || l_specflat==" ") {
        printlog ("ERROR - GSFLAT: output specflat not specified",
            l_logfile, yes)
        nbad = nbad+1
    }
    if (l_redkeep && (l_combflat == "" || l_combflat==" ")) {
        printlog ("ERROR - GSFLAT: output combflat not specified",
            l_logfile, yes)
        nbad = nbad+1
    }

    # Check that the output file does not already exist. If so, exit.
    gimverify(l_specflat) ; l_specflat=gimverify.outname
    if (gimverify.status != 1) {  
    printlog("ERROR - GSFLAT: Output file "//l_specflat//" already exists.",
    l_logfile,yes)
    nbad=nbad+1
    }

    # Check that the output combined flatfield does not already exist. If so, 
    # exit.
    gimverify (l_combflat) ; l_combflat=gimverify.outname
    if (gimverify.status  != 1) {  
        printlog ("ERROR - GSFLAT: Output file "//l_combflat//" already \
            exists.", l_logfile,yes)
        nbad = nbad+1
    }

    if (l_sci_ext =="" || l_sci_ext ==" ") {
        printlog ("ERROR - GSFLAT: extension name sci_ext is missing",
            l_logfile,yes)
        nbad = nbad+1
    }

    #If dq/var propogation is requested, make sure the names are given 
    if(l_fl_vardq) { 
        if (l_dq_ext=="" || l_dq_ext ==" ") {
            printlog ("ERROR - GSFLAT: extension name dq_ext is missing",
                l_logfile, yes)
            nbad = nbad+1
        } else if (l_var_ext=="" || l_var_ext ==" ") { 
            printlog ("ERROR - GSFLAT: extension name var_ext is missing",
                l_logfile, yes)
            nbad = nbad+1
        } 
    }

    #check no commas in sci_ext, var_ext and dq_ext
    if (stridx(",",l_sci_ext)>0 || stridx(",",l_var_ext)>0 || 
        stridx(",",l_dq_ext)>0 ) {
        
        printlog ("ERROR - GSFLAT: sci_ext, var_ext or dq_ext contains commas, \
            give root name only", l_logfile, yes)
        nbad = nbad+1
    }

    #Now check for slit-function file if slitcorr=yes is specified in the 
    #parameters
    if (l_slitcorr) {
        gimverify (l_slitfunc) ; l_slitfunc = gimverify.outname
        if (gimverify.status != 0) { 
            printlog ("ERROR - GSFLAT:  Input slit-function file : "//\
                l_slitfunc//" either does not exist", l_logfile, yes)
            printlog ("                 or is not a valid MEF.", l_logfile, yes)
            nbad = nbad+1
        } else { 
            imgets (l_slitfunc//"[0]", "GSSLITFU") 
            if (imgets.value == "0" || imgets.value == " ") {
                printlog ("ERROR - GSFLAT: Input file "//l_slitfunc//" was \
                    not produced by GSSLITFUNCTION.", l_logfile, yes)
                nbad = nbad+1
            }
        }
    }
    # slit correction not yet implemented for detector by detector fitting
    if (l_slitcorr && l_fl_detec) {
        printlog ("ERROR - GSFLAT: Slit Function correction not implemented \
            for detector by detector fitting.", l_logfile, yes)
        nbad = nbad+1
    }

    #check on existence of grating and filter databases, and chipgaps
    if (!access(l_gratingdb)) { 
        printlog ("ERROR - GSFLAT: gratings database file : "//l_gratingdb//\
            " does not exist.", l_logfile, yes)
        nbad = nbad+1
    }
    if (!access(l_filterdb)) { 
        printlog ("ERROR - GSFLAT: filters database file : "//l_filterdb//\
            " does not exist.", l_logfile, yes)
        nbad = nbad+1
    }

    if ((l_gaindb != "") && (l_gaindb != "default")) {
        if (access(l_gaindb) == no) {
            printlog ("ERROR - GSFLAT: gain database file : "//l_gaindb//\
                " does not exist.", l_logfile, yes)
            nbad = nbad+1
        }
    }
    
    if (l_bpmfile != "" && !access(l_bpmfile)) { 
        printlog ("ERROR - GSFLAT: Chip gap definition file : "//l_bpmfile//\
            " does not exist.", l_logfile, yes)
        nbad = nbad+1
    } else if (l_bpmfile != "") {
        fields (l_bpmfile,"1,2",lines="2") | scan(x11,x12)
        fields (l_bpmfile,"1,2",lines="3") | scan(x21,x22)
        l_sample = "1:"//str(x11-1)//","//str(x12+1)//":"//str(x21-1)//","//str(x22+1)//":6218"
    } else
        l_sample = "*"

    #check that the order for fitting the normalization is rational for 
    #requested reduction technique (slit vs. detector at a time)
    
    norder = 0
    for (n_i=1; n_i<=6; n_i+=1)
        l_rdorder[n_i] = 0
        
    print (l_order) | tokens("STDIN", newl-) | match(",","STDIN",stop+) | \
        count("STDIN") | scan(norder)
    print (l_order) | tokens("STDIN", newl-) | match(",","STDIN",stop+) | \
        fields("STDIN",1,lines="1") | scan(l_rorder)
    if (l_rorder < 1) {
        printlog ("ERROR - GSFLAT: order of fitting function must be greater \
            than zero", l_logfile, yes)
        nbad = nbad+1
    }
    if (norder!=1 && norder!=3 && norder!=6) {
        printlog ("ERROR - GSFLAT: order of fitting function must contain \
            either 1, 3, or 6 elements.", l_logfile, yes)
        nbad = nbad+1
    } else if((norder==3 || norder==6) && !l_fl_detec && (l_rorder > 0)) {
        printlog ("WARNING - GSFLAT: Fitting normalization slit by slit; \
            will use only first", l_logfile, yes)
        printlog ("                  element, ["//str(l_rorder)//"], of \
            fitting order array, ["//l_order//"]", l_logfile,yes)
    } else if (norder==1 && l_fl_detec) {
        for (n_i=1; n_i<=6; n_i+=1)
            l_rdorder[n_i] = l_rorder
    } else if ((norder==3 || norder==6) && (l_rorder > 0)) {
        florder = yes
        for (n_i=1; n_i<=norder; n_i+=1) {
            print (l_order) | tokens("STDIN", newl-) | \
                match(",","STDIN",stop+) | fields("STDIN",1,lines=str(n_i)) | \
                scan(l_rdorder[n_i])
            if (l_rdorder[n_i] < 1)
                florder=no
        }
        if (florder==no) {
            printlog ("ERROR - GSFLAT: order of fitting function must be \
                greater than zero", l_logfile, yes)
            nbad = nbad+1
        }
    } 

    #Assuming we're ok so far, check that all images in the list exist, and 
    # are of the correct filetype

    if (nbad == 0) { 

        #Create list of input flat images 
        files (l_inflats, > temp1)
        scanfile = temp1 

        nim = 0
        while (fscan(scanfile,img) != EOF) {
            nim = nim+1 
            gimverify (img) ; img = gimverify.outname//".fits"
            imgraw = gimverify.outname//".fits"
            if (gimverify.status>=1 && l_rawpath!="") {
                gimverify (l_rawpath//img)
                imgraw = gimverify.outname//".fits"
            }
            if (gimverify.status != 0) { 
                if (gimverify.status == 1) { 
                    printlog ("ERROR - GSFLAT: Input image : "//img//" does \
                        not exist.", l_logfile, yes)
                    nbad = nbad+1 
                } else if (gimverify.status == 2) { 
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is in \
                        OIF format, not MEF.", l_logfile, yes)
                    nbad = nbad+1
                } else if (gimverify.status == 3) { 
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is a \
                        GEIS file, not MEF.", l_logfile, yes)
                    nbad = nbad+1
                } else if (gimverify.status == 4) { 
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is a \
                        simple FITS file, not MEF.", l_logfile, yes)
                    nbad = nbad+1
                }
            }

        }
        # Exit if problems found
        if (nbad > 0) {
            printlog ("ERROR - GSFLAT: "//nbad//" fatal errors found.",
                l_logfile, yes)
            goto error
        }

        scanfile = temp1 

        nim = 0
        while (fscan(scanfile,img) != EOF) {
            nim = nim+1 
            gimverify (img) ; img = gimverify.outname//".fits"
            imgraw = gimverify.outname//".fits"
            if (gimverify.status>=1 && l_rawpath!="") {                
                gimverify(l_rawpath//img)
                imgraw=gimverify.outname//".fits"
            }

            # Check for existing GSREDUCEd image 
            if (imaccess("gs"//img)) {
                imgets ("gs"//img//"[0]","GSREDUCE")
                if (imgets.value!="" && imgets.value!=" " && imgets.value!="0") {
                    printlog ("GSFLAT: Using image gs"//img//" previously \
                        processed with GSREDUCE", l_logfile, l_verbose)
                    img = "gs"//img
                    infile[nim] = img
                } else {
                    printlog ("ERROR - GSFLAT: Found output image gs"//img//\
                        " which is not GSREDUCEd", l_logfile, yes)
                    nbad = nbad+1
                }
            } else {
                infile[nim] = img
                img = imgraw
            }

            nextens = 0
            hselect (img//"[0]","NEXTEND",yes) | scan(nextens)
            if (nextens==0) {
                fxhead (img, format_file="", long_header-, count_lines-) | \
                    match ("STDIN", "IMAGE", stop-) | count | scan(nextens)
            }

            if (nim == 1) 
                next1 = nextens
            if (nextens != next1) {
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same number of extensions", l_logfile, yes)
                nbad = nbad+1
            } 

            imgets (img//"[0]", "MASKTYP", >& "dev$null") 
            if (imgets.value == "1" )
                specmode = yes
            else
                specmode = no
            if (!specmode) { 
                printlog ("ERROR - GSFLAT: Input file "//img//" not of type \
                    LONGSLIT or MOS.", l_logfile, yes)
                nbad = nbad+1 
            }

            # -------------
            # check if GSREDUCE has been run if reductions are not 
            # requested here
            
            imgets (img//"[0]", "GSREDUCE", >& "dev$null")
            if (imgets.value == " " || imgets.value == "" || \
                imgets.value == "0") {
                
                gsflag = no
                if (nim==1)
                    gsfl1=no
                else if (gsflag!=gsfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GSREDUCE'd images", l_logfile, yes)
                    nbad = nbad+1
                }
            } else  {
                gsflag = yes
                if (nim == 1)
                    gsfl1 = yes
                else if (gsflag!=gsfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GSREDUCE'd images", l_logfile, yes)
                    nbad = nbad+1
                }
            }

            if (!gsflag && (!l_fl_over && !l_fl_bias &&!l_fl_dark &&!l_fl_trim)) {
                printlog ("ERROR - GSFLAT: Flat Image - "//img//" has not \
                    been GSREDUCED", l_logfile, yes)
                printlog ("                and no reduction flags are \
                    requested here.", l_logfile, yes)
                nbad = nbad+1
            }

            if (gsflag)
                if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim) { 
                    printlog ("WARNING - GSFLAT: Image "//img//" has been \
                        processed with GSREDUCE", l_logfile, yes)
                    printlog ("                  Ignoring the processing flags",
                        l_logfile,yes)
                    l_fl_over=no
                    l_fl_bias=no
                    l_fl_dark=no
                    l_fl_trim=no
                }

            imgets (img//"[0]", "GMOSAIC", >& "dev$null")
            if (imgets.value != " " && imgets.value != "" && imgets.value!="0") {
                printlog ("WARNING - GSFLAT: Image "//img//" has been \
                    processed with GMOSAIC", l_logfile, yes)
                mosflag = yes
                if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim) { 
                    printlog ("                  Ignoring the processing flags",
                        l_logfile,yes)
                    l_fl_over=no
                    l_fl_bias=no
                    l_fl_dark=no
                    l_fl_trim=no
                } else {
                    printlog ("                  No mosaicing done",
                        l_logfile,yes)
                    l_fl_gmosaic = no
                }
                if (nim==1)
                    mosfl1 = yes
                else if (mosflag!=mosfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GMOSAIC'd images", l_logfile, yes)
                    nbad = nbad+1
                }
                if (l_fl_detec) {
                    printlog ("ERROR - GSFLAT: fl_detec set for fitting \
                        detector by detector, but input image is GMOSAIC'd",
                        l_logfile, yes)
                    nbad = nbad+1
                }
            } else {
                mosflag = no
                if (nim==1)
                    mosfl1 = no
                else if (mosflag!=mosfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GMOSAIC'd images", l_logfile, yes)
                    nbad = nbad+1
                }
            }

            # -------------
            # check exposure times
            imgets (img//"[0]", l_key_exptime, >& "dev$null") 
            if (imgets.value == "0" || imgets.value == " ") {
                printlog ("ERROR - GSFLAT: Image header parameter "//\
                    l_key_exptime//" not found", l_logfile, yes)
                nbad = nbad+1
            }
            if (nim == 1)
                firstexp = real(imgets.value)
            else {
                expdiff = (real(imgets.value)-firstexp)/firstexp
                if (expdiff > 0.1)
                    printlog ("WARNING - GSFLAT: "//img//" - exposure time \
                        different by more than 10%", l_logfile, yes)
            }


            #Check that images in list all have the same FILTER1 and 
            #GRATING ID's and the same central wavelength. 
            hselect (img//"[0]", "FILTER2,GRATING", yes) | scan (filtn, gratn)
            if (nim == 1) {
                filt1 = filtn
                grat1 = gratn
            }
            if (filtn != filt1) {
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same filter ID.", l_logfile, yes)
                nbad = nbad+1
            }
            if (gratn != grat1){
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same grating ID.", l_logfile, yes)
                nbad = nbad+1
            }


        }  #End of while loop over images
        scanfile = ""

    } #End of if-loop for "nbad=0 so far"


    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - GSFLAT: "//nbad//" fatal errors found.",
            l_logfile, yes)
        goto error
    }

    #END OF BASIC CHECKS ...
    #Make sure the input names has .fits on them
    
    delete (temp1, verif-, >& "dev$null")
    for (n_i=1; n_i<=nim; n_i+=1)
        files (infile[n_i], >> temp1)
    scanfile = temp1

    #NOW, call GSREDUCE to overscan, bias-subtract, and dark-subtract if 
    #requested (most likely)
    
    reducing = no
    if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim) {
        reducing = yes
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Calling GSREDUCE.", l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)
        gsreduce ("@"//temp1, outimages="", outpref="gs", logfile=l_logfile,
            verbose=l_verbose, fl_gmosaic=l_fl_gmosaic, fl_fixpix=l_fl_fixpix,
            fl_over=l_fl_over, fl_bias=l_fl_bias, fl_dark=l_fl_dark, fl_flat-,
            fl_cut-, bias=l_bias, dark=l_dark, flatim="", fl_title=no,
            bpmfile=l_bpmfile, key_exptime=l_key_exptime,
            key_biassec=l_key_biassec, key_datasec=l_key_datasec,
            fl_vardq=l_fl_vardq, sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext, key_mdf=l_key_mdf, mdffile=l_mdffile,
            mdfdir=l_mdfdir, bpm=l_bpm, key_ron=l_key_ron, key_gain=l_key_gain,
            ron=l_ron, gain=l_gain, sat=l_sat, ovs_flinter=l_ovs_flinter,
            ovs_med=l_ovs_med, ovs_func=l_ovs_func, ovs_order=l_ovs_order,
            ovs_lowr=l_ovs_lowr, ovs_highr=l_ovs_highr, ovs_niter=l_ovs_niter,
            fl_trim=l_fl_trim, gaindb=l_gaindb, rawpath=l_rawpath,
            gratingdb=l_gratingdb, filterdb=l_filterdb, gradimage="",
            refimage="")

        if (gsreduce.status != 0) { 
            printlog ("ERROR - there was a problem running GSREDUCE.  \
                Stopping now.", l_logfile, yes)
            nbad = nbad+1
            goto error
        } 

        printlog ("GSFLAT: Returned from GSREDUCE.", l_logfile, l_verbose)
        for (n_i=1; n_i<=nim; n_i+=1)
            files("gs"//infile[n_i], >> temp2)
        scanfile = temp2
    }

    #Start Process of combining flats 
    #i loop is over the # of images

    n_i = 0
    while (fscan(scanfile,img) !=EOF) { 

        n_i = n_i+1

        printlog ("GSFLAT: Image #"//n_i//" - Working on : "//img,
            l_logfile, l_verbose)

        #check gain and readnoise values - input parameters vs. header 
        #keywords/values.

        if ( l_key_gain != "" ) {
            keyfound = ""
            if ( l_fl_detec )
                hselect (img//"["//l_sci_ext//",1]", l_key_gain, yes) | \
                    scan (keyfound)
            else
                hselect (img//"[0]", l_key_gain, yes) | scan (keyfound)

            if (keyfound == "") {
                printlog ("WARNING - GSFLAT: keyword "//l_key_gain//" not \
                    found in "//img, l_logfile, l_verbose)
                printlog ("                  Using GAIN = "//str(l_gain),
                    l_logfile, yes)
                l_sgain = l_gain
            } else
                l_sgain = real(keyfound)
        } else
            l_sgain = l_gain

        if ( l_key_ron != "" ) {
            keyfound = ""
            if ( l_fl_detec )
                hselect (img//"["//l_sci_ext//",1]", l_key_ron, yes) | \
                    scan (keyfound)
            else
                hselect (img//"[0]", l_key_ron, yes) | scan (keyfound)

            if (keyfound == "") {
                printlog ("WARNING - GSFLAT: keyword "//l_key_ron//" not found \
                    in "//img, l_logfile, l_verbose)
                printlog ("                  Using RON = "//str(l_ron),
                    l_logfile, yes)
                l_sron = l_ron
            } else
                l_sron = real(keyfound)
        } else
            l_sron = l_ron

        #Make list for imcombine of science extensions 
        print (img, >> scilist)

    } #end of while loop

    #Copy MDF file of the first image to disk for later use
    #Note that the MDFs of the input images should all be identical
    #First check that the MDF exists
    fields (scilist, 1, lines="1") | scan(img)
    tinfo (img//"[MDF]", ttout-, >& "dev$null")
    if (tinfo.tbltype!="fits") {
        printlog ("ERROR - GSFLAT: MDF file does not exist", l_logfile, yes)
        nbad = nbad+1
        goto error
    }
    tcopy (img//"[MDF]", mdf, verbose=no, >& "dev$null")

    #Now, combine images
    #special cases for low numbers of images

    printlog (" ", l_logfile, l_verbose)

    if (nim == 1) {
        printlog ("WARNING - GSFLAT: only one image.", l_logfile, yes)
        copy (img, l_combflat//".fits", verbose-)
        if (l_fl_vardq)  
            combsig = infile[1]//"["//l_var_ext//"]"

    } else {
        if (nim <= 5){
            if (nim == 2){
                printlog ("WARNING - GSFLAT: only combining two images, \
                 turning off rejection", l_logfile, yes)
                printlog ("                  and setting combine type to \
                    average.", l_logfile, yes)
                l_reject = "none"
                l_combine = "average" 
            } else
                printlog ("GSFLAT: combining five or less images.",
                    l_logfile, yes)
        }

        #gemcombine takes care of fl_vardq setting
#JXP Kludge
#            reject=l_reject, offsets="none", masktype=l_masktype,
#            maskvalue=l_maskvalue, bpmfile="", scale=l_scale,
        gemcombine ("@"//scilist, l_combflat, combine=l_combine,
            reject=l_reject, offsets="none", 
            bpmfile="", scale=l_scale,
            zero=l_zero, nrejfile="", weight=l_weight, lsigma=l_lsigma,
            hsigma=l_hsigma, statsec=l_statsec, expname=l_key_exptime,
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, nlow=l_nlow,
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, key_ron=l_key_ron,
            key_gain=l_key_gain, ron=l_sron, gain=l_sgain, snoise=l_snoise,
            sigscale=l_sigscale, pclip=l_pclip, grow=l_grow, logfile=l_logfile,
            fl_vardq=l_fl_vardq)
        if (l_fl_vardq)
            combsig = l_combflat//"["//l_var_ext//"]"

        printlog ("GSFLAT: GEMCOMBINE complete.  Output in file : "//\
            l_combflat, l_logfile, l_verbose)
    }  # End of more than one image

    delete (scilist, verify-, >& "dev$null")

    printlog (" ", l_logfile, l_verbose)


    ## Branch point for detector by detector or slit by slit fitting
    if (l_fl_detec) {	## fitting detector by detector
    
        ##l_combflat is a MEF file with 3 or 6 science extensions
        ##At this point, we aren't going to do anything with potential 
        ##VAR/DQ planes
        
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Begin fitting response functions for detectors",
            l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)

        copy (l_combflat//".fits", tmpflat//".fits")

        imgets (tmpflat//".fits[0]", "NSCIEXT")

        if (int(imgets.value)==6 && (l_rdorder[4] < 1 || l_rdorder[5] < 1 || l_rdorder[6] < 1)) {
            printlog ("WARNING- GSFLAT: order of fitting function must be \
                greater than zero", l_logfile, yes)
            printlog ("WARNING- GSFLAT: setting order of fitting function \
                to "//str(l_rdorder[1])//" for all amps.", l_logfile, yes)
            l_rdorder[2] = l_rdorder[1]
            l_rdorder[3] = l_rdorder[1]
            l_rdorder[4] = l_rdorder[1]
            l_rdorder[5] = l_rdorder[1]
            l_rdorder[6] = l_rdorder[1]
        }

        for (n_ccd=1; n_ccd<=int(imgets.value); n_ccd+=1) { 
            printlog ("GSFLAT: fitting response functions for CCD"//str(n_ccd),
                l_logfile, l_verbose)
            fit1d (tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                l_combflat//"["//l_sci_ext//","//str(n_ccd)//",overwrite]",
                type="ratio", axis=1, interactive=l_rinteractive, naverage=1,
                function=l_rfunction, order=l_rdorder[n_ccd],
                low_reject=l_rlowrej, high_reject=l_rhighrej,
                niterate=l_niterate, grow=1)
            if (l_fl_vardq) {
                imexpr ("a*b**2/c**2",l_combflat//"["//l_var_ext//","//str(n_ccd)//",overwrite]",
                    tmpflat//"["//l_var_ext//","//str(n_ccd)//"]",
                    l_combflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                    tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]",verbose-)
                imcopy (tmpflat//"["//l_dq_ext//","//str(n_ccd)//"]",
                    l_combflat//"["//l_dq_ext//","//str(n_ccd)//",overwrite]",
                    verbose-)
            }
        }

        suf = substr (l_specflat, strlen(l_specflat)-3, strlen(l_specflat))
        if (suf!="fits")
            l_specflat = l_specflat//".fits"
        rename (l_combflat//".fits", l_specflat, field="all")
        rename (tmpflat//".fits", l_combflat//".fits", field="all")

    } else {		## fitting slit by slit
    
        ##l_combflat is a MEF file with a single mosaiced science extension.

        #VAR/DQ made outside the loop
        #The combined flat is l_combflat. 
        #For LONGSLIT mode the MDF will have 3 lines (2 bridges, three spectra 
        #sections) We want to keep these in the same SCI extensions however 
        #even though we will do some of the calibrations separately for each 
        #section.  GSCUT will now do this rather than this task, but some 
        #things we still need to take care of here.

        tinfo (mdf, tbltype="fits", subtype="binary", ttout=no)
        nslits = tinfo.nrows
        printlog ("GSFLAT : Found "//str(nslits)//" Slits in MDF.",
            l_logfile, l_verbose) 

        # make dummy images for inserting response sections for both science 
        # and variance extensions (if needed)
        
        printlog (" ", l_logfile, l_verbose)
        combflatsciext = l_combflat//"["//l_sci_ext//",1]"
        imarith (combflatsciext,"/",combflatsciext, response, verbose-)

        #Call gscut to get the images sections for each spectrum in the MDF 
        fl_delgrad = no
        l_gradimage = ""
        if (l_fl_usegrad)
            l_gradimage = l_combflat

        printlog ("GSFLAT: Calling GSCUT to determine image sections for \
            spectra.", l_logfile, l_verbose)
        scisec = mktemp("tmpscisec")
        gscut (l_combflat, outimag=scisec, fl_update+, secfile=specsecfile,
            logfile=l_logfile, fl_vard=l_fl_vardq, gratingdb=l_gratingdb,
            filterdb=l_filterdb, verbose=l_verbose, xoffset=l_xoffset,
            yoffset=l_yoffset, yadd=l_yadd, refimage=l_refimage,
            gradimage=l_gradimage) 

        if (!access(specsecfile) ) { 
            printlog ("ERROR - GSFLAT: GSCUT failed in some way. ",
                l_logfile, yes)
            nbad = nbad+1
            goto error
        } 
        printlog("GSFLAT: Returned from GSCUT.", l_logfile, l_verbose)

        # Set the sample to avoid the chip gaps, if possible
        if (l_sample!="*") {
            hselect (l_combflat//"["//l_sci_ext//",1]", "CCDSUM", yes) | \
                translit ("STDIN", '"', "", delete+) | scan(Xbin)
            hselect (l_combflat//"["//l_sci_ext//",1]", "i_naxis1", yes) | \
                scan(Xmax)
            x11 = int((x11-1)/Xbin+0.5)
            x12 = int((x12-1)/Xbin+0.5)
            x21 = int((x21-1)/Xbin+0.5)
            x22 = int((x22-1)/Xbin+0.5)
            l_sample = "1:"//str(x11-1)//","//str(x12+1)//":"//str(x21-1)//","//str(x22+1)//":"//Xmax
        }
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Begin fitting response functions for slit(s)",
            l_logfile, l_verbose)
        printlog ("GSFLAT: Sample in fitting "//l_sample,
            l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)

        #Now we process each slitlet 
        count (specsecfile) | scan (nslits)

        scanfile = specsecfile
        n_i = 1
        while(fscan(scanfile,specsec) !=EOF) { 
            # Temporary FITS files re-used within this loop must have
            # a different name at each iteration (FITS Kernel cache is 
            # based on file names)
            
            flatfit = mktemp("tmpfit")
            sciflat = mktemp("tmpflat")
            inscoo = mktemp("tmpinscoo")

            # Get x1,y1 from MDF, MDFROW from gscut image
            hselect (scisec//"["//l_sci_ext//","//str(n_i)//"]",
                "MDFROW", yes) | scan(n_row)
            tprint (scisec//".fits[mdf]", columns="SECX1,SECY1", row=str(n_row),
                prparam-, prdata+, pwidth=80, plength=0, showr-, orig_row+,
                showhdr-, showunits-, option="plain", align+, sp_col="", 
                lgroup=0) | scan(x1,y1)
            print (x1," ",y1, > inscoo) 

            if (l_verbose)
                printlog ("GSFLAT: Fitting response for slit#"//str(n_i)//\
                    " - SpecSec = "//specsec, l_logfile, l_verbose)

            fit1d (scisec//"["//l_sci_ext//","//str(n_i)//"]", sciflat,
                type="ratio", axis=1, interactive=l_rinteractive,
                sample=l_sample, naverage=1, function=l_rfunction,
                order=l_rorder, low_reject=l_rlowrej, high_reject=l_rhighrej,
                niterate=l_niterate, grow=0.)
            if (l_slitcorr)
                imarith (sciflat,"*", 
                    l_slitfunc//"["//l_sci_ext//"]"//specsec,
                    sciflat, verbose-)
            # insert sections into larger images
            iminsert (response, sciflat, response, "replace", coordfile=inscoo,
                offset1=0, offset2=0, xcol="c1", ycol="c2")
            imdelete (sciflat, verify-, >& "dev$null")
            delete (inscoo, verify-)
            n_i = n_i+1 

            # Give the user a chance to fit the rest non-interactively
            if (l_rinteractive && nslits>1) {
                l_fl_answer = fl_answer
                l_rinteractive = l_fl_answer
            }


        } #end of while loop over slits 
        scanfile = ""
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Finished Loop over slits. Now packing output MEF...",
            l_logfile, l_verbose)

        # update headers
        hedit (response, "BPM", "", add-, addonly-, delete+, verify-,
            show-, update+)

        #Copy the PHU from the first of the "reduced" flatfields 
        suf = substr (l_specflat, strlen(l_specflat)-3, strlen(l_specflat))
        if (suf!="fits") 
            l_specflat = l_specflat//".fits"

        wmef (input=mdf, output=l_specflat, extnames="MDF",
            phu=l_combflat//".fits", verbose-, >& "dev$null")
        if (wmef.status != 0) { 
            printlog ("ERROR - GSFLAT: problem writing output MEF.",
                l_logfile, yes)
            nbad = nbad+1
            goto error 
        }

        #Now the science extension/flatfield
        fxinsert (response//".fits", l_specflat//"[1]", groups="", verbose=yes,
            >& "dev$null")

        gemhedit (l_specflat//"[2]", "EXTNAME", l_sci_ext, "Extension name")
        gemhedit (l_specflat//"[2]", "EXTVER", 1, "Extension number")
        hedit (l_specflat//"[2]", "BPM", value="", add-, addonly-, delete+,
            verify-, show-, update+)

        #The VAR/DQ planes if requested ...
        if (l_fl_vardq) {
            imexpr ("a*b**2/c**2",l_specflat//"["//l_var_ext//",1,append]",
                l_combflat//"["//l_var_ext//",1]",
                l_specflat//"["//l_sci_ext//",1]",
                l_combflat//"["//l_sci_ext//",1]", verbose-)
            imcopy (l_comflat//"["//l_dq_ext//",1]",
                l_specflat//"["//l_dq_ext//",1,append]", verbose-)
        }

    }  ## end of slit by slit fitting section

    #Update the PHU
    gemdate ()
    hedit (l_specflat//"[0]", "BPM", "", add-, addonly-, delete+, verify-,
        show-, update+)
    gemhedit (l_specflat//"[0]", "NCOMBINE", nim,
        "Number of images in IMCOMBINE")
    gemhedit (l_specflat//"[0]", "GSFLAT", gemdate.outdate, 
        "UT Time stamp for GSFLAT")
    gemhedit (l_specflat//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Time of last modification with GEMINI")
    if( l_slitcorr)
        gemhedit (l_specflat//"[0]", "GSLITCOR", gemdate.outdate,
            "UT Time stamp for Slit Function Correction")
    if (l_fl_detec)
        gemhedit (l_specflat//"[0]", "GSFMODE", "DETECTOR",
            "GSFLAT fitting mode")
    else
        gemhedit (l_specflat//"[0]", "GSFMODE", "SLIT", "GSFLAT fitting mode")
        
    gemhedit (l_specflat//"[0]", "GSFFCT", l_rfunction,
        "GSFLAT fitting function")
    gemhedit (l_specflat//"[0]", "GSFORDER", l_order, "GSFLAT fitting order")

    gemhedit (l_specflat//"[0]", "GAINORIG", l_sgain, "Original Gain value")
    gemhedit (l_specflat//"[0]", "RONORIG", l_sron, "Original Read-noise value")
    gaineff = real(nim)*l_sgain
    roneff = l_sron/sqrt(real(nim))
    if (l_key_gain == "")
        l_key_gain="GAIN"
    if (l_key_ron == "")
        l_key_ron="RON"
    gemhedit (l_specflat//"[0]", l_key_gain, gaineff, "Rescaled Gain value")
    gemhedit (l_specflat//"[0]", l_key_ron, roneff, "Rescaled Read-noise value")

    # clean up
    goto clean

error:
    status = 1
    printlog (" ", l_logfile, yes) 
    printlog ("ERROR - GSFLAT: Program execution failed with "//\
        str(nbad)//" errors.", l_logfile, verbose=yes) 

clean:
    scanfile = ""
    if (status==0) { 
        printlog (" ", l_logfile, yes) 
        date | scan(sdate)
        printlog ("GSFLAT done  "//sdate, l_logfile, yes)
    }
    if (!l_redkeep)
        delete (l_combflat//".fits", ver-, >& "dev$null")
    imdelete (scisec, verif-, >& "dev$null")
    imdelete (response, verif-, >& "dev$null")
    delete (mdf//","//temp1//","//temp2//","//specsecfile, ver-, >& "dev$null")

    # close log file
    printlog ("-------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)

end
