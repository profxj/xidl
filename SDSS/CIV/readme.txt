;;;;;;;;;
;; readme.txt
;; 3 Jun 2011
;; 
;; This file describes use of SDSS CIV/SiIV/CaII/... search recipe for
;; $XIDL_DIR/SDSS/CIV codes.  You can run all of the steps on a
;; partial or full sample with sdss_runpipeline.
;;
;; History:
;;   03 Jun 11  by KLC


$XIDL_DIR/SDSS/CIV Contents:
        Misc.
        --dr7qso_srt.fit.gz: Schneider et al. (2010) DR7 QSO catalog
        of 105783 objects, sorted by plate "string"
        (four character, with 0's prepended to reach length).
	--dr7qso_srt.list: appropriately formatted spectrum list file in same
	order as dr7qso_srt.fit
	--dr7qso_noBAL.fit.gz/.list: 99569 subset of Schneider et al.
	    structure above and its list.
	--dr7qso_BAL_cat.fit.gz: 6214 BAL QSO subset from Shen et al. (2011)
            value-added DR7 QSO catalog.
	--dr7qso_BAL.fit.gz/.list: 6214 subset of Schneider et al.
	    structure and its list.
	--sdss_balqso_map.sav: Simple file that maps standard QSO
        structure to their flags.
	--dr7qso_splfail.list: spectrum list where inverse-variance
	weighted spline continuum fit could not be accomplised due to
	low S/N.  
	--todo.txt: list of changes to be made to code
	--sdss_civ.lst: list of transitions to plot in sdss_chkciv and
	sdss_fixciv. 
	--sdss_runparallelsrch.sh: shell script to run sdss_fndlin or 
            sdss_civsearch on multiple processors.

	Codes.
        --sdss_chkciv.pro: GUI to rate candidates.
        --sdss_civsearch.pro: Finds the candidates from generic
        absorption ilne list; instantiates sdsscivstrct.
        --sdss_ewciv.pro: Measures EW of candidates.
        --sdss_fixciv.pro: GUI to adjust bounds and EW measures (even
        continuum).
        --sdss_fndciv.pro: Called by sdss_civsearch to do the
        candidate search; fills qalcharstrct.
        --sdss_fndlin.pro: Finds the generic absorption lines redward
        of Lya (and also fits continua); instantiates sdsscontistrct.
        --sdss_functions.pro: Generic gaggle of little codes used by
        everything else. 
        --sdss_runpipeline.pro: Runs the full suite of codes on random
        sample. 
        --sdsscivstrct__define.pro: Doublet structure for candidates
        and for real lines.
        --sdsscontistrct__define.pro: Continuum structure that also
        holds automatically detected absorption lines. 
        --Extras: readcol w/ count= keyword


1. Set SDSSPATH environment variable (e.g. $HOME/SDSS) and SDSSDR (e.g. DR7_QSO).
    Everything will be read/written with respect to $SDSSPATH/$SDSSDR/
    Many files following the naming convention:
    1d_26/pppp/1d/spSpec-jjjjj-pppp-fff*.fit
    pppp is plate, jjjjj is spectroscopic MJD, and fff is the fiber number.
    All may be gzipped.
    Spectra are in $SDSSPATH/$SDSSDR/spectro/ and have nothing extra before .fit.
    Continuum fits are in $SDSSPATH/$SDSSDR/conti/ and -eigconti.fit are from	
    PCA fits and -splconti.fit are from bspline fits.
    Automatically detected lines are in $SDSSPATH/$SDSSDR/abslin/ and have
    -abslin.fit.

    I also recommend adding the following lines to your .idlstartup file:
        !quiet = 1 ; suppress messages
        .r sdss_functions
        !quiet = 0
    This enables all the functions in sdss_functions to be immediately callable. 


2. Retrieve DR7_QSO/, which includes:
        a)  spectro/ (8.0Gb and 105783 files)
             Contains quasar spectra and the header redshifts are not trusted.
	     You can retrieve individual files with wget from:
	     http://das.sdss.org/spectro/1d_26/pppp/1d/spSpec*.fit
        b) conti/ (9.6Gb and 2x105783 files)
            Continuum fit for each spectrum based on PCA 
	    eigenspectra (for more information, contact M Kao) and 
	    based on bspline that can be generated with sdss_fitspline (in
            sdss_functions). 
 	    If spline files exist, can short cut sdss_fndlin with /usesplfil.
	    Can make spline files with:
idl>sdss_fndlin_fitallspl,getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list',getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit',pixlim=7,/clobber


3. Create or retrieve abslin/ directory (13Gb and 105783 files).
     These files contain the sdsscontistrct structures that have all
     information of the continuum fit and the automatically detected
     centroids. 
     To create, run sdss_fndlin like:
idl>sdss_fndlin,getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list',getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit',/eig,pixlim=7,/tcpu,/clobber
      The spectro files do not have the right redshifts written to their
      headers, must defer to the QSO structure.


4. Search for the desired doublet with sdss_civsearch. Use dblt_name=
    to define desired doublet (passed to dblt_retrieve()).  
idl>sdss_civsearch,getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list',getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit','sdss_civcand_full.fits',/tcpu
    The output file has the candidates in first extension. (Previously
    there were two extensions.)


[Note: use $XIDL_DIR/SDSS/CIV/sdss_runparallelsrch.sh to cull out the
jobs onto multiple processors by starting new IDL sessions:
To run sdss_fndlin on 4 processors:
xterm>$XIDL_DIR/SDSS/CIV/sdss_runparallelsrch.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list
$XIDL_DIR/SDSS/CIV/dr7qso_srt.fit "/tcpu,/clobber,pixlim=7" 
To run sdss_civsearch on 2 processors, for SiIV:
xterm>$XIDL_DIR/SDSS/CIV/sdss_runparallelsrch.sh 2 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list
$XIDL_DIR/SDSS/CIV/dr7qso_srt.fit $SDSSPATH/$SDSSDR/sdss_siivcand_full.fits "/tcpu,dblt_name='SiIV'"
Woot!]


5. Visually evaluate candidate doublets with sdss_chkciv GUI.
idl>sdss_chkciv,'sdss_civcand_full.fits','sdss_civcand_rated.fits'
    The output will be a single-extension file with the sdsscivstrct,
    now with the rating[0] tag filled in (see sdss_getrating() in
    sdss_functions). 
    If sdss_civcand_rated.fits already exists, have the option to
    restart from that, so long it is a copy of the input.


6. Optional: sort the output from Step 5 before plugging it into the
    next GUI.  For example, don't bother any more with the rejected
    ones. 


7. Modify the candidates with sdss_fixciv. May not need this step...

