;+ 
; NAME:
; hires_edgeflat
;     Version 1.1
;
; PURPOSE:
;  This is mainly a frontend to x_edgeflat.  See that routine for
;  extensive (and the most recent) details
;
; CALLING SEQUENCE:
;  hires_edgeflat, hires, setup, [chip], /CHK, 
;
; INPUTS:
;   hires     -  HIRES structure
;   setup    -  Setup identifier 
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to hires_fittflat to create a 2D
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
;   hires_edgeflat, hires, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;   hires_edgeflat_clean
;   hires_getfil
;   x_edgeflat
;
; REVISION HISTORY:
;   Feb-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro hires_edgeflat, hires, setup, chip, CHK=chk, SMSHROW=smshrow, $
                    CLOBBER=clobber, P_NSIG=p_nsig, THRESH=thresh, $
                    DEBUG = debug, $
                    INTER=inter, NSIG=nsig, MINFLAT=minflat, NOBASIS=nobasis, $
                    MNSEP=mnsep, $
                    NOSTOP=nostop, NOSORDR=nosordr, IEXTRAP=iextrap, _EXTRA=EXTRA

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_edgeflat, hires, setup, [chip], P_NSIG=, THRESH=, /CLOBBER,'
      print, '           /CHK, /INTER, /DEBUG, NSIG=, MINFLAT= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if n_elements(CHIP) EQ 0 then chip = [1L,2L,3L]
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
  hires_rdxlog, 'flat_log', 'hires_edgeflat', /ophdr

  for ii=0L,n_elements(chip)-1 do begin
      qq = chip[ii]
      ;; Chip
      case qq of 
          -1: begin
              mssg= 'hires_edgeflat: Tracing Single Chip'
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [0L,0L]
          end
          1: begin
              mssg= 'hires_edgeflat: Tracing BLUE CCD' 
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,0L]
          end
          2: begin
              mssg= 'hires_edgeflat: Tracing GREEN CCD' 
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,1L]
          end
          3: begin
              mssg= 'hires_edgeflat: Tracing RED CCD' 
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,1L]
              gd = where(hires.setup EQ setup AND hires.chip EQ qq $
                         AND hires.flg_anly NE 0 AND hires.type EQ 'TFLT')
              TNSIG=100.
              ;; 2nd order?
              if hires[gd[0]].cross EQ 'RED' $
                AND not keyword_set(NOSORDR) then begin
                  mssg2= 'hires_edgeflat:  Attempting to block 2nd order light'
                  hires_rdxlog, 'flat_log', mssg2
                  print, 'hires_edgeflat:  Set /NOSORDR to stop me!'
                  ;; Use green if it exists
                  gostr = hires_getfil('ordr_str',setup, CHIP=2)
                  np = n_elements(gostr[0].lhedg)
                  MNSEP = gostr[1].lhedg[np/2]-gostr[0].lhedg[np/2] 
                  SORDR=1
              endif
          end
          else: stop
      endcase

      ;; Log
      hires_rdxlog, 'flat_log', mssg
      mssg = 'hires_edgeflat: Extrap = '+strjoin(strtrim(extrap,2))
      hires_rdxlog, 'flat_log', mssg

      ;; Check output file
      trc_fil = hires_getfil('tflat_str', setup, CHIP=qq, /name, CHKFIL=chkf)
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'hires_edgeflat: Trace exists, moving on..'
          continue
      endif

      ;; Order structure
      ordr_fil = hires_getfil('ordr_str', setup, CHIP=qq, CHKFIL=chkf, /name)
      if CHKF NE 0 AND $
        not keyword_set( CLOBBER ) then begin
          print, 'hires_edgeflat: Order structure exists. Moving on!', ordr_fil
          continue
      endif


      ;; Read in flat 
      flat = hires_getfil('qtz_fil', setup, CHIP=qq, FIL_NM=flat_fil)
      flativar = xmrdfits(flat_fil, 1, /silent)
      sz = size(flat, /dimensions)

      ;; Set binning
      cbin = round(2048. / sz[0])
      if qq LT 0 then rbin = round(2048. / sz[1]) else rbin = round(4096. / sz[1])

      ;; Output files
      qafil = hires_getfil('qa_trcflat', setup, CHIP=qq)
      logfil = hires_getfil('flat_log')

      if qq EQ -1 then begin
         smshrow = round(1100./rbin) ;; Avoid the inkspot
         N_x0fit = 5
      endif

      x_edgeflat, flat, cbin, trc_fil, ordr_fil, NOOVERLAP=nooverlap, $
                  CHK=chk, QAFIL=qafil, FLATIVAR=flativar, SORDR=sordr, MNSEP=mnsep, $
                  EXTRAP=extrap, _EXTRA=EXTRA, LOGFIL=logfil, SMSHROW=smshrow, $ 
                  N_x0fit=n_x0fit
      
      ;; Write 
;      mwrfits, trc_str, trc_fil, /create
      print, 'hires_edgeflat: Check the QA file == ', qafil

  endfor

  ;; 
  print, 'hires_edgeflat: You may now wish to check the results ' + $
    'with hires_chktrcflat'
  print, 'hires_edgeflat: All done!'

  return
end
