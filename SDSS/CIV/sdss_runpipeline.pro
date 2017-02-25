;+ 
; NAME:
; sdss_runpipeline
;    Version 2.0
;
; PURPOSE:
;    Run full pipeline
;
; CALLING SEQUENCE:
;   searchForCIV2 [dblt_name=, nlist=, zlim=, absdir=, list_fil=,
;            and_fil=, rated_fil=, bad_fil=, /eig, /plot, /usesplfil, /clobber,
;            debug, /help] 
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   list_fil= -- SDSS DR7 QSO catalo file list [default: dr7qso_pipelinetest.list]
;   [sdssstrct_fil] -- SDSS DR7 QSO catalog structure in order of
;                      list_fil; same name as list_fil with .fit [so
;                      default: dr7qso_pipelinetest.fit] 
;   cand_fil= -- FITS file with 2 extensions from sdss_civsearch
;   rated_fil= -- FITS file of rated CIV candidates from sdss_chkciv
;   fixed_fil= -- FITS file of fixed CIV doublets from sdss_fixciv
;   bad_fil= -- FITS file of bad CIV candidates from sdss_chkciv
;
; OPTIONAL KEYWORDS:
;   dblt_name= -- ion name (default: 'CIV')
;   nlist= -- number of QSOs to search (default: 10L)
;   zlim= -- 2-element array of low, high redshift to search (default:
;            none unless dblt_name is 'MgII' or 'CaII', then set to
;            small range where those detectable in SDSS)
;   absdir= -- name of output directory, will be created if DNE
;              (default: abslint/) 
;   list_fil=, cand_fil=, rated_fil=, bad_fil=
;   /eig -- eigen-continuum is default
;   /plot -- show useful plots
;   /usesplfil -- use SPLCONTI spline continuum (faster)
;   /clobber -- overwrite everything
;   /debug -- print lots of information (hopefully useful)
;   /help -- print out syntax
;   _extra= -- passed to sdss_fndlin and sdss_civsearch; see code
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_runpipeline,/eig,/debug,/usesplfil,/clobber
;
; PROCEDURES/FUNCTIONS CALLED:
;   dblt_retrieve()
;   sdss_functions (sdss_getrandqsosmpl, directly)
;   sdss_fndlin
;   sdss_ewciv
;   sdss_civsearch
;   sdss_chkciv
;   sdss_fixciv
;
; REVISION HISTORY:
;   02-Jun-2011 Major revamp to be generic and run with new
;               structures, KLC
;-
;------------------------------------------------------------------------------

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_runpipeline, dblt_name=dblt_name, keeplist=keeplist, nlist=nlist, $
                      zlim=zlim,absdir=absdir,list_fil=list_fil,cand_fil=cand_fil, $
                      rated_fil=rated_fil, fixed_fil=fix_fil, bad_fil=bad_fil, $
                      eig=eig, plot=plot, usesplfil=usesplfil,clobber=clobber, $
                      debug=debug, help=help, _extra=extra

  if keyword_set(help) then begin
     print,'Syntax - searchForCIV2 [dblt_name=, nlist=, /keeplist, zlim=, absdir=, '
     print,'                  list_fil=, cand_fil=, rated_fil=, bad_fil=, /eig, '
     print,'                  /plot, /usesplfil, /clobber, /debug, /help, _extra=]'
     return
  endif 

  ;; Figure out where all the time goes!
  profiler,/reset
  profiler
  profiler,/system

  ;; set defaults
  sdssdir = sdss_getsdssdir()
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  dblt = dblt_retrieve(dblt_name)
  if not keyword_set(nlist) then nlist = 10L ; number to test
  if not keyword_set(zlim) then begin
     ;; be picky about range for high-rwave doublets
     wvlim = sdss_getspecwave()
     case strlowcase(dblt.ion) of 
        'caii': zlim = [0.,wvlim[1]/dblt.wvI-1] ; 1.34
        'mgii': zlim = [0.,wvlim[1]/dblt.wvI-1] ; 2.29
        else: zlim = 0                  ; search range, array of two
     endcase 
  endif else if n_elements(zlim) ne 2 then $
     stop,'sdds_runpipeline: zlim must be [low,high] array'
  if not keyword_set(absdir) then absdir = 'abslint/' ; don't mess with real stuff
  ;; Does directory exist?
  test = file_search(sdssdir+absdir,count=ntest)
  if ntest eq 0 then spawn,'mkdir -p '+sdssdir+absdir

  ;; Y/N flags (1: yes, 0: no) 
  if not keyword_set(eig) then eig = 0 
  if not keyword_set(clobber) then clobber = 0 
  if not keyword_set(debug) then debug = 0
  if not keyword_set(plot) then plot = 0           ; show stuff
  if not keyword_set(usesplfil) then usesplfil = 0 ; cheat the spline; faster

  ;; File names
  if not keyword_set(list_fil) then list_fil = 'dr7qso_pipelinetest.list'
  sdssstrct_fil = strmid(list_fil,0,strpos(list_fil,'.',/reverse_search))+'.fit'
  if not keyword_set(cand_fil) then $
     cand_fil = strlowcase(dblt.ion)+'cand_pipelinetest.fit'
  if not keyword_set(rated_fil) then $
     rated_fil = strlowcase(dblt.ion)+'rated_pipelinetest.fit'
  if not keyword_set(fixed_fil) then $
     fixed_fil = strlowcase(dblt.ion)+'fixed_pipelinetest.fit'
  if not keyword_set(bad_fil) then $
     bad_fil = strlowcase(dblt.ion)+'bad_pipelinetest.fit'

  ;; Clear some space
  print,''
  print,''

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; STEP ZERO
  ;; Generate the sample, which can be a random draw or the full
  ;; sample (in scrambled order) and/or for a specific redshift range
  ;; (for effective searching).
  ;; _extra= includes /bal, rmaglim=
  test = file_search(sdssstrct_fil,count=ntest)
  if (ntest eq 0 or keyword_set(clobber)) and not keyword_set(keeplist) then begin
     sdss_getrandqsosmpl,nlist,sdssstrct_fil,list=list_fil,zlim=zlim,$
                         absdir=absdir,_extra=extra
     print,''
     print,''
  endif 
  
  
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; STEP ONE
  ;; First find wavelengths of interest (absorption) using
  ;; sdss_qsolin.pro. Toggle eigenspectra fits on and off
  ;; using /EIG
  ;; Note: pixlim= for pca_qsotilt and lsnr= for sdss_fndlin_srch()
  ;; /tcpu prints out time information
  ;; /redo makes sure to do the actual finding even if clobber not set
  ;; _extra= includes istrt=, /nocompress
  sdss_fndlin,list_fil, sdssstrct_fil, /tcpu,/redo,clobber=clobber,eig=eig, $
              pixlim=7, debug=debug, cchk=plot, usesplfil=usesplfil,_extra=extra
  print,''
  print,''
  
  
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; STEP TWO
  ;; Now use sdss_civsearch.pro to find candidate doublets
  ;; _extra= includes zmin=, zmax=, nfil=, rmax=, inifil=
  ;;         and lsnr2=, dvtol=, dvgal=, dvqso=, dvbal=, dz=, gapsz=
  ;;         (to sdss_fndciv) 
  sdss_civsearch, list_fil, sdssstrct_fil, cand_fil, /tcpu, debug=debug, $
                  dblt_name=dblt_name, _extra=extra
  ;; Read back in
;  qstr = xmrdfits(cand_fil,1,/silent)
  civcand = xmrdfits(cand_fil,1,/silent)
  
  
  if size(civcand,/type) eq 8 then begin 
     print,'' 
     print,'' 
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Step THREE
     ;; Now use sdss_chkciv.pro to give grades.
     ;; _extra includes altdblt_name= (default: MgII)
     sdss_chkciv,cand_fil,rated_fil,dblt_name=dblt_name,debug=debug,_extra=extra
     civstr = xmrdfits(rated_fil,1,/silent)
     test = where(civstr.rating[0] eq sdss_getrating(/unrated),ntest)
     if ntest ne 0 then $
        print,'sdss_runpipeline: left some candidates unrated ',ntest
     gd = where(civstr.rating[0] ne sdss_getrating(/unrated) and $
                civstr.rating[0] ne sdss_getrating(/bad),ngd)

     if ngd ne 0 then begin
        mwrfits,civstr[gd],fixed_fil,/create,/silent
        print,'sdss_runpipeline: created file of CIV doublets to fix ',fixed_fil
  
        ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; STEP FOUR
        ;; Now use sdss_fixciv.pro to check that system are CIV
        ;; absorption and to adjust fits.  To move boundaries of
        ;; fits, click with left mouse button to set left boundary and
        ;; right button to set right boundary.  When the plot stops
        ;; changing, you're at the end of your list and you
        ;; should just exit out by pressing "Done"
        print,'' 
        print,'' 
        print,'sdss_runpipeline: last step (sdss_fixciv) not working'
        ;; sdss_fixciv may over write
;     sdss_fixciv, fixed_fil
        endif 
  endif else if not keyword_set(debug) then $
     print,'sdss_runpipeline: no candidates found; nothing to plot'
  

  ;; Finish and report
  profiler,/report,filename='sdss_runpipeline.profile'
  print,'sdss_runpipeline: IDL profiler report saved to sdss_runpipeline.profile'

  if keyword_set(debug) then stop,'sdss_runpipeline: stop before exiting'

end

