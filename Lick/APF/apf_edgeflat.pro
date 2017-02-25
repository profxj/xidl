;+ 
; NAME:
; apf_edgeflat
;     Version 1.1
;
; PURPOSE:
;  This is mainly a frontend to x_edgeflat.  See that routine for
;  extensive (and the most recent) details
;
; CALLING SEQUENCE:
;  apf_edgeflat, apf, setup, [chip], /CHK, 
;
; INPUTS:
;   apf     -  HIRES structure
;   setup    -  Setup identifier 
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to apf_fittflat to create a 2D
;  solution.  Typical name:  'Flats/TStr_B_01.fits'
;
; OPTIONAL KEYWORDS:
;  /CHK -- Check the order edges interactively (similar to INTER)
;  /INTER -- Determine the order edges interactively
;  SMSHROW -- Row where order edges are determined (default: 1/2 way
;             up the image)
;  THRESH  -- Threshold for an order edge on the red side
;              (default: 100.)
;  /CLOBBER -- Overwrite the output structure
;  P_NSIG  --  Number of sigma significance that an order edge should
;             have for the red side (default: 50.)
;  NSIG  --  Number of sigma significance that an order edge should
;             have fo the blue side (default: 5.0)
;  MINFLAT -- Mininum counts in flat on blue side to use orders (100)
;  IEXTRAP -- Number of orders to extract in extrapolation off each
;             edge, e.g. [1,1]  to do one off each edge of the chip
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_edgeflat, apf, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;   apf_edgeflat_clean
;   apf_getfil
;   x_edgeflat
;
; REVISION HISTORY:
;   Feb-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro apf_edgeflat, apf, setup, CHK=chk, SMSHROW=smshrow, $
                    CLOBBER=clobber, P_NSIG=p_nsig, THRESH=thresh, $
                    DEBUG = debug, $
                    INTER=inter, NSIG=nsig, MINFLAT=minflat, NOBASIS=nobasis, $
                    MNSEP=mnsep, $
                    NOSTOP=nostop, NOSORDR=nosordr, IEXTRAP=iextrap, _EXTRA=EXTRA

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_edgeflat, apf, setup, [chip], P_NSIG=, THRESH=, /CLOBBER,'
      print, '           /CHK, /INTER, /DEBUG, NSIG=, MINFLAT= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( P_NSIG) then p_nsig = 50.
  if not keyword_set( NSIG) then nsig = 5.0
  if not keyword_set( THRESH) then thresh = 100.
  if not keyword_set( BGUESS) then BGUESS = 75L
  if not keyword_set( NCOEFF) then ncoeff = 6L
  if not keyword_set( MINFLAT) then minflat = 100.
  if not keyword_set( NOOVERLAP ) then nooverlap = 4.
  if NOOVERLAP LT 0 THEN NOOVERLAP = 0

  x_psclose
  !p.multi=[0,1,1]

; Loop on chip
  ;apf_rdxlog, 'flat_log', 'apf_edgeflat', /ophdr

  if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [0L,0L]
  ;; Log
  ;apf_rdxlog, 'flat_log', mssg
  ;mssg = 'apf_edgeflat: Extrap = '+strjoin(strtrim(extrap,2))
  ;apf_rdxlog, 'flat_log', mssg

  ;; Check output file
  trc_fil = apf_getfil('tflat_str', setup, /name, CHKFIL=chkf)
  if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
     print, 'apf_edgeflat: Trace exists, moving on..'
     return
  endif

  ;; Order structure
  ordr_fil = apf_getfil('ordr_str', setup, CHKFIL=chkf, /name)
  if CHKF NE 0 AND $
     not keyword_set( CLOBBER ) then begin
     print, 'apf_edgeflat: Order structure exists. Moving on!', ordr_fil
     return
  endif


  ;; Read in flat 
  flat = apf_getfil('qtz_fil', setup, CHIP=qq, FIL_NM=flat_fil)
  flativar = xmrdfits(flat_fil, 1, /silent)
  sz = size(flat, /dimensions)

  ;; Set binning
  cbin = round(2048. / sz[0])
  rbin = round(4096. / sz[1])
  
  ;; Output files
  qafil = apf_getfil('qa_trcflat', setup, CHIP=qq)
  logfil = apf_getfil('flat_log')
  
  x_edgeflat, flat, cbin, trc_fil, ordr_fil, TNSIG=25., NOOVERLAP=nooverlap, $
              CHK=chk, QAFIL=qafil, FLATIVAR=flativar, SORDR=sordr, MNSEP=mnsep, $
              EXTRAP=extrap, _EXTRA=EXTRA, LOGFIL=logfil
  
  ;; Write 
  print, 'apf_edgeflat: Check the QA file == ', qafil

  ;; 
  print, 'apf_edgeflat: You may now wish to check the results ' + $
    'with apf_chktrcflat'
  print, 'apf_edgeflat: All done!'

  return
end
