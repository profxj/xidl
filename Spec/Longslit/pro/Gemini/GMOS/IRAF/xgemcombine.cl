# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure gemcombine (input, output) 

# Combine MEF files
# 
# Version Feb 28, 2002  BM  v1.3 release
#         Aug 10, 2002  BM,IJ  v1.4 release
#         Nov 20, 2002  IJ allow rejection with only 2 input images
#         Aug 12, 2003  KL IRAF2.12 - hedit: addonly- ; imcombine: headers,
#                           bpmasks,expmasks,outlimits ="", rejmask->rejmasks,
#                           plfile->nrejmasks
#	  Nov 07, 2003  KL change default nkeep from 0 to 1

string input        {"", prompt = "Input MEF images"}
string output       {"", prompt = "Output MEF image"}

string title        {"", prompt = "Title for output SCI plane"}
string combine      {"average", enum = "average|median", prompt = "Combination operation"}
string reject       {"avsigclip", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip", prompt = "Rejection algorithm"}
string offsets      {"none", prompt = "Input image offsets"}
string  masktype    {"none", enum="none|goodvalue", prompt = "Mask type"}
real    maskvalue   {0., prompt = "Mask value"}
string scale        {"none", prompt = "Image scaling"}
string zero         {"none", prompt = "Image zeropoint offset"}
string weight       {"none", prompt = "Image weights"}
string statsec      {"[*,*]", prompt = "Statistics section"}
string expname      {"EXPTIME", prompt = "Exposure time header keyword"}
real   lthreshold   {INDEF, prompt = "Lower threshold"}
real   hthreshold   {INDEF, prompt = "Upper threshold"}
int    nlow         {1, min = 0, prompt = "minmax: Number of low pixels to reject"}
int    nhigh        {1, min = 0, prompt = "minmax: Number of high pixels to reject"}
int    nkeep        {1, prompt = "Minimum to keep or maximum to reject"}  
bool   mclip        {yes, prompt = "Use median in sigma clipping algorithms?"}
real   lsigma       {3., min = 0, prompt = "Lower sigma clipping factor"}
real   hsigma       {3., min = 0, prompt = "Upper sigma clipping factor"}
string key_ron      {"RDNOISE", prompt = "Keyword for readout noise in e-"}
string key_gain     {"GAIN", prompt = "Keyword for gain in electrons/ADU"}
#JXP -- Added mask parameters
#string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
#real    maskvalue   {0., prompt="Mask value"}
real   ron          {0.0, min = 0., prompt = "Readout noise rms in electrons"}
real   gain         {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}
string snoise       {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}
real   sigscale     {0.1, min = 0., prompt = "Tolerance for sigma clipping scaling correction"}
real   pclip        {-0.5, prompt = "pclip: Percentile clipping parameter"}
real   grow         {0.0, min = 0., prompt = "Radius (pixels) for neighbor rejection"}
string bpmfile      {"", prompt = "Name of bad pixel mask file or image."}
string nrejfile     {"", prompt = "Name of rejected pixel count image."}
string sci_ext      {"SCI", prompt = "Name(s) or number(s) of science extension"}
string var_ext      {"VAR", prompt = "Name(s) or number(s) of variance extension"}
string dq_ext       {"DQ", prompt = "Name(s) or number(s) of data quality extension"}
bool   fl_vardq     {no, prompt = "Make variance and data quality planes?"}
string logfile      {"gemcombine.log", prompt = "Log file"}
bool   fl_dqprop    {no, prompt = "Propagate all DQ values?"}
bool   verbose      {yes, prompt = "Verbose output?"}
int    status       {0, prompt = "Exit status (0=good)"}
struct  *flist      {"", prompt = "Internal use only"} 

begin
        
        string l_input, l_output, l_statsec, l_expname, l_bpmfile
        string l_logfile, l_offsets, l_masktype
        string l_combine, l_reject, l_scale, l_zero, l_weight
        string l_key_ron, l_key_gain, l_snoise
        string l_sgain, l_sron
        real   l_ron, l_gain, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
        real   l_grow, l_sigscale, l_pclip, l_maskvalue
        int    l_nlow, l_nhigh, l_nkeep
        bool   l_verbose, s_verbose, l_mclip, l_fl_dqprop, l_fl_vardq, debug
        string l_sci_ext, l_var_ext, l_dq_ext, l_title, tsci, tvar, tdq
        string sextlist, vextlist, dextlist, extlist, sci[600], var[600], dq[600]
        string filelist, combout, combdq, combsig, img, dqsum, dqfits, dqsumold, suf
        string tmpdq, scilist, dqlist, plimg, tmplog, sect, infile[500], mdf
        string exphu, phuimg
        real mean, epadu, rdnoise, gaineff, roneff
        int i, j, k, l, n, idum, nbad, len, nsci, nvar, ndq, junk
        int nimg, nextnd, secpos1, secpos2, nsciext
        bool useextver, usemdf, l_cleanbpm
        struct sdate, line
        char	l_nrejfile = ""
        
        # Query parameters
        l_input = input; l_output = output
        l_title = title
        l_logfile = logfile
        l_combine = combine; l_reject = reject; l_offset = offset
#JXP kludge (June 2009)
#        l_masktype= masktype
#        l_maskvalue = maskvalue
        l_masktype = "none"
        l_maskvalue = 0.
        l_scale = scale; l_zero = zero; l_weight = weight
        l_statsec = statsec; l_expname = expname
        l_lthreshold = lthreshold; l_hthreshold = hthreshold
        l_nlow = nlow; l_nhigh = nhigh; l_nkeep = nkeep
        l_lsigma = lsigma; l_hsigma = hsigma
        l_verbose = verbose; l_mclip = mclip
        l_key_gain = key_gain; l_key_ron = key_ron
        l_gain = gain; l_ron = ron; l_snoise = snoise 
        l_sgain = l_key_gain; l_sron = l_key_ron
        l_sigscale = sigscale; l_pclip = pclip; l_grow = grow
        l_sci_ext = sci_ext; l_dq_ext = dq_ext; l_var_ext = var_ext
        l_bpmfile = bpmfile; l_fl_dqprop = fl_dqprop; l_fl_vardq = fl_vardq
        junk = fscan (nrejfile, l_nrejfile)
        
        status = 0
        useextver = no
        usemdf = no
        debug = no

        if (debug) print ("in gemcombine")
        
        # Keep imgets parameters from changing by outside world
        cache ("imgets", "hedit", "keypar", "gimverify", "gemextn", "gemdate") 
        
        s_verbose = hedit.show
        hedit.show = l_verbose
        
        # Start logging to file
        if (l_logfile == "STDOUT") {
            l_logfile = ""
            l_verbose = yes
        }
        
        date | scan (sdate) 
        printlog ("--------------------------------------------------------------------------------", l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE -- " // sdate, l_logfile, verbose = l_verbose) 
        printlog (" ", l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: input = " // l_input, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: output = " // l_output, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: title = " // l_title, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: bpmfile = " // l_bpmfile, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: sci_ext = " // l_sci_ext, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: var_ext = " // l_var_ext, l_logfile, verbose = l_verbose) 
        printlog ("GEMCOMBINE: dq_ext = " // l_dq_ext, l_logfile, verbose = l_verbose) 
        printlog (" ", l_logfile, verbose = l_verbose) 
        
        # Verify that the extension names are not empty, otherwise exit gracefully
        if (l_sci_ext == "" || l_sci_ext == " ") {
            printlog ("ERROR - GEMCOMBINE: Extension name sci_ext is missing", l_logfile, yes) 
            goto error
        }
        if (l_dq_ext == "" || l_dq_ext == " ") {
            printlog ("ERROR - GEMCOMBINE: Extension name dq_ext is missing", l_logfile, yes) 
            goto error
        } 
        
        # Define temporary files
        filelist = mktemp ("tmpfilelist") 
        scilist = mktemp ("tmpscilist") 
        dqlist = mktemp ("tmpdqlist")
        dqsum = mktemp ("tmpdqsum") 
        tmpdq = mktemp ("tmpdq") 
        tmplog = mktemp ("tmplog") 
        sextlist = mktemp ("tmpsextlist") 
        vextlist = mktemp ("tmpvextlist") 
        dextlist = mktemp ("tmpdextlist") 
        extlist = mktemp ("tmpextlist") 
        mdf = mktemp ("tmpmdf") 
	
        # temporary files wannabes (values assigned later)
        combout = ""
        combdq = ""
        combsig = ""
        dqfits = ""
        

        if ("" != l_nrejfile) {
            if (no == l_fl_vardq) {
                printlog ("WARNING - GEMCOMBINE: nrejfile given, but \
                    fl_vardq-, so no output.", l_logfile, verbose+)
            } else {
                gemextn (l_nrejfile, check="absent", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="", replace="", outfile="dev$null", logfile="",
                    glogpars="", verbose=l_verbose)
                if (gemextn.fail_count > 0 || 1 != gemextn.count) {
                    printlog ("ERROR - GEMCOMBINE: problem with nrejfile.", \
                        l_logfile, verbose+)
                    goto error
                }
                gemextn (l_nrejfile, check="absent", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="exten", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) | scan (line)
                l_nrejfile = line
            }
        }
	

        #check that there are input files
        if (l_input == "" || l_input == " ") {
            printlog ("ERROR - GEMCOMBINE: input files not specified", l_logfile, yes) 
            goto error
        }
        
        # check existence of list file
        if (substr (l_input, 1, 1) == "@") {
            len = strlen (l_input) 
            if (no==access (substr (l_input, 2, len) ) ) {
                printlog ("ERROR - GEMCOMBINE: " // substr (l_input, 2, len) // " does not exist.", l_logfile, yes) 
                goto error
            }
        }
        
        # Check for an image section
        sect = ""
        secpos1 = stridx ("[", l_input) 
        secpos2 = stridx ("]", l_input) 
        if (secpos1 > 0 && secpos2 > 0) {
            sect = substr (l_input, secpos1, secpos2) 
            l_input = substr (l_input, 1, secpos1-1) 
        } else if (secpos1 > 0 || secpos2 > 0) {
            printlog ("ERROR - GEMCOMBINE: mismatched brackets in " // l_input, l_logfile, yes) 
            goto error
        }
        files (l_input, sort-, > filelist) 
        
        #check that an output file is given
        if (l_output == "" || l_output == " ") {
            printlog ("ERROR - GEMCOMBINE: output files not specified", l_logfile, yes) 
            goto error
        }
        
        # Check that the output file does not already exist. If so, exit.
        if (imaccess (l_output) ) {  
            printlog ("ERROR - GEMCOMBINE: Output file " // l_output // " already exists.", l_logfile, yes) 
            goto error
        }
        
        #put suffix on output image
        len = strlen (l_output) 
        suf = substr (l_output, (len-4), len) 
        if (suf != ".fits") {
            l_output = l_output // ".fits"
        }
        
        # Make list of sci/var/dq extensions
        files (l_sci_ext, sort-, > sextlist) 
        count (sextlist) | scan (nsci) 
        files (l_var_ext, sort-, > vextlist) 
        count (vextlist) | scan (nvar) 
        files (l_dq_ext, sort-, > dextlist) 
        count (dextlist) | scan (ndq) 
        if ( (nsci != nvar) || (nsci != ndq) ) {
            printlog ("ERROR - GEMCOMBINE: different numbers of SCI, VAR, or DQ extensions.", l_logfile, yes) 
            goto error
        }
        
        if (nsci > 600) {
            printlog ("ERROR - GEMCOMBINE: the number of science extensions is larger than sci, var, and dq array lengths", l_logfile, yes) 
            goto error
        }
        
        join (sextlist // "," // vextlist // "," // dextlist, "", output = extlist) 
        delete (sextlist // "," // vextlist // "," // dextlist, verify-, >& "dev$null") 
        flist = extlist
        i = 0
        while (fscan (flist, tsci, tvar, tdq) != EOF) {
            i = i+1
            sci[i] = tsci
            var[i] = tvar
            dq[i] = tdq
        }
        flist = ""
        delete (extlist, verify-, >& "dev$null") 
        
        # Check that all images in the list exist and are MEF
        nbad = 0
        k = 0
        flist = filelist
        while (fscan (flist, img) != EOF) {
            secpos1 = stridx ("[", img) 
            secpos2 = stridx ("]", img) 
            if (secpos1 > 0 && secpos2 > 0) {
                sect = substr (img, secpos1, secpos2) 
                img = substr (img, 1, secpos1-1) 
            } else if (secpos1 > 0 || secpos2 > 0) {
                printlog ("ERROR - GEMCOMBINE: mismatched brackets in " // img, l_logfile, yes) 
                goto error
            }
            gimverify (img) 
            if (gimverify.status > 0) {
                nbad = nbad+1
            } else {
                k = k+1
                # name w/o suffix
                infile[k] = gimverify.outname // ".fits"
                if (k == 1) {
                    imgets (infile[k] // "[0]", "NSCIEXT", >& "dev$null") 
                    if (debug) print (infile[k] // "[0]" // " nsciext: " // imgets.value) 
                    nsciext = int (imgets.value) 
                    if (nsciext > 600) {
                        printlog ("ERROR - GEMCOMBINE: the number of science extensions is larger than sci, var, and dq array lengths", l_logfile, yes) 
                        goto error
                    }
                    # check for MDF
                    imgets (infile[k] // "[0]", "NEXTEND", >& "dev$null") 
                    if (debug) print (imgets.value) 
                    idum = int (imgets.value) 
                    for (l = 1; l <= idum; l += 1) {
                        keypar (infile[k] // "[" // l // "]", "EXTNAME", silent+) 
                        if (keypar.value == "MDF" && keypar.found) {
                            tcopy (infile[k] // "[MDF]", mdf // ".fits", verbose-) 
                            usemdf = yes
                        }
                    }
                    # check if EXTVER used
                    if (nsci == 1 && nsciext > 0) {
                        if (imaccess (infile[k] // "[" // sci[k] // "," // nsciext // "]") ) {
                            tsci = sci[1]
                            tvar = var[1]
                            tdq = dq[1]
                            for (i = 1; i <= nsciext; i += 1) {
                                sci[i] = tsci // "," // i
                                var[i] = tvar // "," // i
                                dq[i] = tdq // "," // i
                            }
                            nsci = nsciext
                            useextver = yes
                        } else {
                            if (debug) print ("strange test failed")
                            if (debug) print (infile[k] // "[" // sci[k] // "," // nsciext // "]")
                        }
                    }
                }
            }    
        }
        #end of while loop
        flist = ""
        nimg = k
        
        # Exit if problems found
        if (nbad > 0) {
            printlog ("ERROR - GEMCOMBINE: " // nbad // " image(s) either do not exist or are not MEF files.", l_logfile, yes) 
            goto error
        }
        
        # Check offset file
        if (l_offset != "none" && !access (l_offset) ) {
            printlog ("ERROR - GEMCOMBINE: " // l_offset // " Offset file not found", l_logfile, yes) 
            goto error
        }
        
        
        # images ok, so continue
        l_cleanbpm = no # clean headers of BPM after combining
        if (debug) print ("nsci: " // nsci)
        for (i = 1; i <= nsci; i += 1) {
            # temporary files used within this loop
            combout = mktemp("tmpcombout")
            combdq = mktemp("tmpcombdq")
            combsig = mktemp("tmpcombsig")
            dqfits = mktemp("tmpdqfits")
	    
            n = 0
            gaineff = 0.0
            roneff = 0.0
            flist = filelist
            while (fscan (flist, img) != EOF) {
                # Check if an image section is present
                secpos1 = stridx ("[", img) 
                secpos2 = stridx ("]", img) 
                if (secpos1 > 0) {
                    sect = substr (img, secpos1, secpos2) 
                    img = substr (img, 1, secpos1-1) 
                    len = secpos1-1
                } else {
                    len = strlen (img) 
                }
                suf = substr (img, len-4, len) 
                # Strip the suffix of the input file:
                if (suf != ".fits" && imaccess (img // ".fits") ) {
                    img = img // ".fits"
                }
                if (no==imaccess (img // "[" // sci[i] // "]") ) {
                    printlog ("ERROR - GEMCOMBINE: Could not access " // img // "[" // sci[i] // "]", l_logfile, yes) 
                    goto error
                }
                n = n+1
                # if n=1, save the name for the output phu
                if (n == 1) {
                    phuimg = img
                }
                # check gain and readnoise values 
                if (l_key_gain != "" && l_key_gain != " ") {
                    if (l_sgain != str (l_gain) ) {
                        imgets (img // "[" // sci[i] // "]", l_key_gain, >& "dev$null") 
                        if (imgets.value == "0") {
                            imgets (img // "[0]", l_key_gain, >& "dev$null") 
                        }
                        if (imgets.value == "0") {
                            # only warn if it's going to be used
                            if (l_reject == "ccdclip" || l_reject == "crreject") {
                                printlog ("WARNING - GEMCOMBINE: keyword " // l_key_gain // " not found in " // img, l_logfile, yes) 
                                printlog ("Using gain = " // str (l_gain), l_logfile, l_verbose) 
                            }
                            l_sgain = str (l_gain) 
                        } else { 
                            l_sgain = imgets.value
                        }
                    } 
                } else {
                    l_sgain = str (l_gain) 
                }
                if (l_key_ron != "" && l_key_ron != " ") {
                    if (l_sron != str (l_ron) ) {
                        imgets (img // "[" // sci[i] // "]", l_key_ron, >& "dev$null") 
                        if (imgets.value == "0") {
                            imgets (img // "[0]", l_key_ron, >& "dev$null") 
                        }
                        if (imgets.value == "0") {
                            # only warn if it's going to be used
                            if (l_reject == "ccdclip" || l_reject == "crreject") {
                                printlog ("WARNING - GEMCOMBINE: keyword " // l_key_ron // " not found in " // img, l_logfile, yes) 
                                printlog ("Using ron = " // l_ron, l_logfile, l_verbose) 
                            }
                            l_sron = str (l_ron) 
                        } else { 
                            l_sron = imgets.value
                        }
                    } 
                } else {
                    l_sron = str (l_ron) 
                }
                gaineff = gaineff+real (l_sgain) 
                roneff = roneff+real (l_sron) **2
                
                # Make sure exposure time is in the science header
                keypar (img // "[" // sci[i] // "]", l_expname, silent+) 
                if (keypar.found == no) {
                    # check the phu
                    keypar (img // "[0]", l_expname, silent+) 
                    if (keypar.found == no) {
                        printlog ("ERROR - GEMCOMBINE: " // l_expname // " not found", l_logfile, l_verbose) 
                        goto error
                    } else {
                        hedit (img // "[" // sci[i] // "]", l_expname, keypar.value, add+, 
                            addonly-, del-, verify-, show-, update+) 
                    }
                }
                # science extension
                print (img // "[" // sci[i] // "]" // sect, >> scilist) 
                # DQ extension
                if (imaccess (img // "[" // dq[i] // "]") ) {
                    imcopy (img // "[" // dq[i] // "]", tmpdq // "_" // n // ".pl", verbose-) 
                    if (l_fl_dqprop) {
                        print (tmpdq // "_" // n // ".pl", >> dqlist) 
                    }
                    hedit (img // "[" // sci[i] // "]", "BPM", tmpdq // "_" // n // ".pl", add+, 
                        addonly-, del-, verify-, show-, update+) 
                    l_cleanbpm = yes
                } else if (imaccess (l_bpmfile) ) {
                    hedit (img // "[" // sci[i] // "]", "BPM", l_bpmfile, add+, addonly-, del-, 
                        verify-, show-, update+) 
                    l_cleanbpm = yes
                }
            }
            flist = ""
            
            #combine images
            #special cases for low numbers of images
            if (n == 1) {
                printlog ("ERROR - GEMCOMBINE: only one image to combine.", l_logfile, yes) 
                goto error
            }
            if (n <= 5) {
                printlog ("WARNING - GEMCOMBINE: five or less images to combine.", l_logfile, yes) 
            }
            if (l_reject == "minmax" && n <= nhigh + nlow) {
                printlog ("ERROR - GEMCOMBINE: Too few images for minmax \
                    parameters (nhigh, nlow).", l_logfile, yes) 
                goto error
            }
            if (l_fl_vardq == no) {
                combdq = ""
                combsig = ""
            }
            if (l_verbose) {
                if (debug) print ("calling imcombine")
                if (debug) print (l_scale)
                if (debug) print (l_zero)
                if (debug) print (l_weight)
                imcombine ("@" // scilist, combout, headers = "", bpmasks = "", rejmasks = "", 
                    nrejmasks = combdq, expmasks = "", sigmas = combsig, logfile = tmplog, 
                    combine = l_combine, reject = l_reject, project = no, outtype = "real", 
                    outlimits = "", offsets = l_offset, masktype = l_masktype, maskvalue = l_maskvalue, 
                    blank = 0, scale = l_scale, zero = l_zero, weight = l_weight, statsec = l_statsec, 
                    expname = l_expname, lthreshold = l_lthreshold, hthreshold = l_hthreshold, 
                    nlow = l_nlow, nhigh = l_nhigh, nkeep = l_nkeep, mclip = l_mclip, lsigma = l_lsigma, 
                    hsigma = l_hsigma, rdnoise = l_sron, gain = l_sgain, snoise = l_snoise, 
                    sigscale = l_sigscale, pclip = l_pclip, grow = l_grow) 
            } else {
                if (debug) print ("calling imcombine")
                imcombine ("@" // scilist, combout, headers = "", bpmasks = "", rejmasks = "", 
                    nrejmasks = combdq, expmasks = "", sigmas = combsig, logfile = tmplog, 
                    combine = l_combine, reject = l_reject, project = no, outtype = "real", 
                    outlimits = "", offsets = l_offset, masktype = l_masktype, maskvalue = l_maskvalue, 
                    blank = 0, scale = l_scale, zero = l_zero, weight = l_weight, statsec = l_statsec, 
                    expname = l_expname, lthreshold = l_lthreshold, hthreshold = l_hthreshold, 
                    nlow = l_nlow, nhigh = l_nhigh, nkeep = l_nkeep, mclip = l_mclip, lsigma = l_lsigma, 
                    hsigma = l_hsigma, rdnoise = l_sron, gain = l_sgain, snoise = l_snoise, 
                    sigscale = l_sigscale, pclip = l_pclip, grow = l_grow, >& "dev$null") 
            }
            if (debug) print ("called imcombine")
            if (access (tmplog) ) {
                if (l_logfile != "" && l_logfile != " ") {
                    type (tmplog, >> l_logfile) 
                }
                if (l_verbose) {
                    type (tmplog) 
                }
                delete (tmplog, verify-, >& "dev$null") 
            }
            
            # New effective gain and readnoise
            if (l_combine == "average") {
                roneff = sqrt (roneff) 
            } else {
                gaineff = 2.*gaineff/3.
                roneff = sqrt (2.*roneff/3.) 
            }
            
            # clean headers
            if (l_cleanbpm) 
                hedit ("@" // scilist, "BPM", "", add-, addonly-, delete+, verify-, show-, update+) 
            # update gain and readnoise
            imgets (combout, "GAINORIG", >& "dev$null") 
            if (imgets.value == "0") {
                gemhedit (combout, "GAINORIG", (real (l_sgain) ), "Input gain") 
            }
            imgets (combout, "RONORIG", >& "dev$null") 
            if (imgets.value == "0") {
                gemhedit (combout, "RONORIG", (real (l_sron) ), "Input read-noise") 
            }
            hedit (combout, l_key_gain, gaineff, add+, addonly-, del-, verify-, show-, update+) 
            hedit (combout, l_key_ron, roneff, add+, addonly-, del-, verify-, show-, update+) 
            
            # make variance image by squaring combsig
            #imarith(combsig,"*",combsig,combsig,verbose-)
            
            # Use an "inverse" bad pixel mask (values are the number of pixels used)
            # to get the variance right.
            if (l_fl_vardq) {
                imcopy (combdq // ".pl", combdq // ".fits", verbose-) 
                if ("" != l_nrejfile) {
                    imcopy (combdq // ".pl", \
                    l_nrejfile // "[" // i // ",append].fits", verbose-)
                }
                delete (combdq // ".pl", verify-)
                imcalc (combsig // ".fits," // combdq // ".fits", combsig // ".fits", "(im1**2/(" // n // "-im2))", pixtyp = "real", verb-) 
                
                imgets (combsig, "GAINORIG", >& "dev$null") 
                if (imgets.value == "0") {
                    gemhedit (combsig, "GAINORIG", (real (l_sgain) **2), "Input gain") 
                }
                imgets (combsig, "RONORIG", >& "dev$null") 
                if (imgets.value == "0") {
                    gemhedit (combsig, "RONORIG", (real (l_sron) **2), "Input read-noise") 
                }
                hedit (combsig, l_key_gain, (gaineff**2), add+, addonly-, del-, verify-, 
                    show-, update+) 
                hedit (combsig, l_key_ron, (roneff**2), add+, addonly-, del-, verify-, 
                    show-, update+) 
                
                # DQ file
                imcalc (combdq, combdq // ".fits", "if (im1 == " // n // ") then 1 else 0", verbose-) 
                # propagate DQ values if requested
                #
                # WARNING: something is wrong in the following block.
                #    The image dqsum is used but it is never create!  KL
                if (access (dqlist) && l_fl_dqprop) {
                    printlog (" ", l_logfile, l_verbose) 
                    printlog ("Propagating DQ values", l_logfile, l_verbose) 
                    flist = dqlist
                    for (j = 1; j <= n; j += 1) {
                        # tmp file used within this loop
                        dqsumold = mktemp("tmpdqsumold")
			
                        idum = fscan (flist, plimg) 
                        imarith (plimg // sect, "*", combdq, plimg, verbose-) 
                        if (j == 1)
                            imcopy (combdq, dqsumold // ".fits", verbose-) 
                        else
                            imrename (dqsum, dqsumold, verbose-) 
                        addmasks (dqsumold // "," // plimg, dqsum // ".fits", "im1 || im2", flags = " ") 
                        imdelete (dqsumold, verify-, >& "dev$null")
                    }
                    flist = ""
                    imdelete (combdq // ".fits", verify-, >& "dev$null") 
                    imrename (dqsum, combdq // ".fits", verbose-) 
                    delete (dqlist, verify-) 
                }
                # Headers
                hedit (combsig // "[0]," // combdq, "BPM", "", add-, addonly-, delete+, verify-, 
                    show-, update+) 
            }
            # update headers
            hedit (combout, "BPM", "", add-, addonly-, delete+, verify-, show-, update+) 
            if (l_title != "" && l_title != " ") {
                hedit (combout, "i_title", l_title, add+, addonly-, del-, verify-, 
                    show-, update+) 
            }
            
            #make output MEF file
            files (sci[i], sort-) | scan (tsci) 
            files (var[i], sort-) | scan (tvar) 
            files (dq[i], sort-) | scan (tdq) 
            if (i == 1) {
                if (l_fl_vardq) {
                    wmef (combout // "," // combsig // "," // combdq, l_output, 
                        extnames = tsci // "," // tvar // "," // tdq, phu = phuimg, 
                        verbose-, >& "dev$null") 
                    if (useextver) {
                        hedit (l_output // "[1]," // l_output // "[2]," // l_output // "[3]", 
                            "EXTVER", i, add+, addonly-, del-, verify-, show-, update+) 
                    }
                } else {
                    if (debug) {
                        print ("calling wmef")
                        wmef (combout, l_output, extnames = tsci, phu = phuimg, 
                            verbose+)
                        print ("wmef completed")
                    } else {
                        wmef (combout, l_output, extnames = tsci, phu = phuimg, 
                            verbose-, >& "dev$null")
                    }
                    if (useextver) {
                        hedit (l_output // "[1]", "EXTVER", i, add+, addonly-, del-, 
                            verify-, show-, update+) 
                    }
                }
                if (wmef.status > 0) {
                    goto error
                }
                imgets (phuimg // "[0]", l_expname, >& "dev$null") 
                exphu = imgets.value
                imgets (combout, l_expname, >& "dev$null") 
                if (real (exphu) != real (imgets.value) ) {
                    hedit (l_output // "[0]", l_expname, imgets.value, add+, addonly-, del-, 
                        verify-, show-, update+) 
                }
                # update gain and readnoise in PHU?
                imgets (phuimg // "[0]", "GAIN", >& "dev$null") 
                if (imgets.value != "0") {
                    imgets (phuimg // "[0]", "GAINORIG", >& "dev$null") 
                    if (imgets.value == "0") {
                        gemhedit (l_output // "[0]", "GAINORIG", (real (l_sgain) ), "Input gain") 
                    }
                    imgets (phuimg // "[0]", "RONORIG", >& "dev$null") 
                    if (imgets.value == "0") {
                        gemhedit (l_output // "[0]", "RONORIG", (real (l_sron) ), "Input read-noise") 
                    }
                    hedit (l_output // "[0]", l_key_gain, gaineff, add-, addonly-, del-, 
                        verify-, show-, update+) 
                    hedit (l_output // "[0]", l_key_ron, roneff, add-, addonly-, del-, 
                        verify-, show-, update+) 
                }
            } else {
                if (l_fl_vardq) {
                    imcopy (combdq, dqfits // ".fits", verbose-) 
                    fxinsert (combout // ".fits," // combsig // ".fits," // dqfits // ".fits", l_output // "[" // ( (i-1) *3) // "]", "0", verbose-, >& "dev$null") 
                    imgets (l_output // "[0]", "NEXTEND") 
                    if (debug) print (imgets.value) 
                    nextnd = int (imgets.value) +3
                    hedit (l_output // "[0]", "NEXTEND", nextnd, add-, addonly-, del-, 
                        verify-, show-, update+) 
                    hedit (l_output // "[" // ( (i-1) *3+1) // "]", "EXTNAME", tsci, add+, 
                        addonly-, del-, verify-, show-, update+) 
                    hedit (l_output // "[" // ( (i-1) *3+2) // "]", "EXTNAME", tvar, add+, 
                        addonly-, del-, verify-, show-, update+) 
                    hedit (l_output // "[" // ( (i-1) *3+3) // "]", "EXTNAME", tdq, add+, 
                        addonly-, del-, verify-, show-, update+) 
                    if (useextver) {
                        hedit (l_output // "[" // ( (i-1) *3+1) // "]", "EXTVER", i, add+, 
                            addonly-, del-, verify-, show-, update+) 
                        hedit (l_output // "[" // ( (i-1) *3+2) // "]", "EXTVER", i, add+, 
                            addonly-, del-, verify-, show-, update+) 
                        hedit (l_output // "[" // ( (i-1) *3+3) // "]", "EXTVER", i, add+, 
                            addonly-, del-, verify-, show-, update+) 
                    }
                    imdelete (dqfits // ".fits", verify-) 
                } else {
                    if (debug) print ("calling fxinsert")
                    fxinsert (combout // ".fits", l_output // "[" // (i-1) // "]", "0", verbose-, >& "dev$null") 
                    if (debug) print ("fxinsert completed")
                    imgets (l_output // "[0]", "NEXTEND") 
                    if (debug) print (imgets.value) 
                    nextnd = int (imgets.value) +1
                    hedit (l_output // "[0]", "NEXTEND", nextnd, add-, addonly-, del-, 
                        verify-, show-, update+) 
                    hedit (l_output // "[" // i // "]", "EXTNAME", tsci, add+, addonly-, del-, 
                        verify-, show-, update+) 
                    if (useextver) {
                        hedit (l_output // "[" // i // "]", "EXTVER", i, add+, addonly-, del-, 
                            verify-, show-, update+) 
                    }
                }
            }
            imdelete (combout // "," // combsig // "," // combdq // "," // dqfits, verify-, >& "dev$null") 
            imdelete (tmpdq // "*.pl", verify-, >& "dev$null") 
            delete (scilist, verify-, >& "dev$null") 
        }
        # attach MDF
        if (usemdf) {
            if (l_verbose) {
                printlog ("GEMCOMBINE: Attaching the MDF", l_logfile, l_verbose) 
            }
            fxinsert (mdf // ".fits", l_output // "[" // (i-1) // "]", "1", verbose-, >& "dev$null") 
            imgets (l_output // "[0]", "NEXTEND") 
            if (debug) print (imgets.value) 
            nextnd = int (imgets.value) +1
            hedit (l_output // "[0]", "NEXTEND", nextnd, add-, addonly-, del-, verify-, 
                show-, update+) 
        }
        
        # Update PHU
        gemdate ()
        nhedit (l_output // "[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", comfile="NULL", after="",
            before="", update+, add+, addonly-, delete-, verify-, show-)
        nhedit (l_output // "[0]", "GEMCOMB", gemdate.outdate, 
            "UT Time stamp for GEMCOMBINE", comfile="NULL", after="",
            before="", update+, add+, addonly-, delete-, verify-, show-) 
        # clean up
        goto clean
        
error:
        {
            status = 1
            goto clean
        }
        
clean:
        {
            delete (filelist // "," // scilist, verify-, >& "dev$null") 
            imdelete (combout // "," // combsig // "," // combdq // ".fits", verify-, >& "dev$null") 
            imdelete (dqsum // "," // mdf // ".fits", verify-, >& "dev$null") 
            imdelete (tmpdq // "*.pl," // combdq // ".pl", verify-, >& "dev$null") 
            delete (sextlist // "," // vextlist // "," // dextlist // "," // extlist, verify-, >& "dev$null") 
            hedit.show = s_verbose
            # close log file
            printlog (" ", l_logfile, l_verbose) 
            printlog ("GEMCOMBINE done", l_logfile, l_verbose) 
            printlog ("--------------------------------------------------------------------------------", l_logfile, l_verbose) 
            flpr
        }
        
end
