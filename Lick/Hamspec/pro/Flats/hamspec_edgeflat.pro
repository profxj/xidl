;+ 
; NAME:
; hamspec_edgeflat
;     Version 1.1
;
; PURPOSE:
;  This is mainly a frontend to x_edgeflat.  See that routine for
;  extensive (and the most recent) details
;
; CALLING SEQUENCE:
;  hamspec_edgeflat, hamspec, setup /CHK, 
;
; INPUTS:
;   hamspec     -  HIRES structure
;   setup    -  Setup identifier 
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to hamspec_fittflat to create a 2D
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
;   hamspec_edgeflat, hamspec, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;   hamspec_edgeflat_clean
;   hamspec_getfil
;   x_edgeflat
;
; REVISION HISTORY:
;   Feb-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro hamspec_edgeflat, hamspec, setup, CHK=chk, SMSHROW=smshrow, $
                    CLOBBER=clobber, P_NSIG=p_nsig, THRESH=thresh, $
                    DEBUG = debug, $
                    INTER=inter, NSIG=nsig, MINFLAT=minflat, NOBASIS=nobasis, $
                    MNSEP=mnsep, $
                    NOSTOP=nostop, NOSORDR=nosordr, IEXTRAP=iextrap, _EXTRA=EXTRA

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_edgeflat, hamspec, setup,  P_NSIG=, THRESH=, /CLOBBER,'
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
;  hamspec_rdxlog, 'flat_log', 'hamspec_edgeflat', /ophdr

  if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,1L]

  ;; Log
;  hamspec_rdxlog, 'flat_log', mssg
;  mssg = 'hamspec_edgeflat: Extrap = '+strjoin(strtrim(extrap,2))
;  hamspec_rdxlog, 'flat_log', mssg

  ;; Check output file
  trc_fil = hamspec_getfil('tflat_str', setup, /name, CHKFIL=chkf)
  if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
     print, 'hamspec_edgeflat: Trace exists, moving on..'
     return
  endif

  ;; Order structure
  ordr_fil = hamspec_getfil('ordr_str', setup, CHKFIL=chkf, /name)
  if CHKF NE 0 AND $
     not keyword_set( CLOBBER ) then begin
     print, 'hamspec_edgeflat: Order structure exists. Moving on!', ordr_fil
     return
  endif


  ;; Read in flat 
  flat = hamspec_getfil('qtz_fil', setup, FIL_NM=flat_fil)
  flativar = xmrdfits(flat_fil, 1, /silent)
  sz = size(flat, /dimensions)

  ;; Set binning
  ise = where(hamspec.setup EQ setup and hamspec.flg_anly EQ 1)
  cbin = hamspec[ise[0]].colbin

  ;; Output files
  qafil = hamspec_getfil('qa_trcflat', setup)
  logfil = hamspec_getfil('flat_log')

  ;; Using INST_STATS for Hamspec.  Odds are this will fail for some 
  ;; users!!  INST_STATS = [RMS, MEAN_POS, MEAN_NEG]  see
  ;; x_trcflat.pro
  WIDEN = 2. ;; Ok for Dewar 6 and 2.5" slit
  if strtrim(hamspec[0].ccd,2) EQ 'e2vCCD203-824kx4kthin' then begin
     ncoeff = 5L
     poly_n = 3
     n_xpfit = 4
     eordr = 3
  endif
  x_edgeflat, flat, cbin, trc_fil, ordr_fil, NOOVERLAP=nooverlap, N_XPFIT=n_xpfit,$
              WIDEN=widen, ncoeff=ncoeff, poly_n=poly_n, $;TNSIG=10., $
              CHK=chk, QAFIL=qafil, FLATIVAR=flativar, SORDR=sordr, MNSEP=mnsep, $
              EORDR=eordr, $
              EXTRAP=extrap, _EXTRA=EXTRA, LOGFIL=logfil, INST_STATS=[0.001, 0., 0.]

      
  ;; Write 
  print, 'hamspec_edgeflat: Check the QA file == ', qafil

  ;; 
  print, 'hamspec_edgeflat: You may now wish to check the results ' + $
    'with hamspec_chktrcflat'
  print, 'hamspec_edgeflat: All done!'

  return
end
