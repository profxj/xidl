;+ 
; NAME:
; uves_edgeflat
;     Version 1.1
;
; PURPOSE:
;  This is mainly a frontend to x_edgeflat.  See that routine for
;  extensive (and the most recent) details
;
;
; CALLING SEQUENCE:
;   
;  uves_edgeflat, uves, setup, [side], /CHK, 
;
; INPUTS:
;   uves     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to uves_fittflat to create a 2D
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
;             have fo the blue side (default: 2.0)
;  MINFLAT -- Mininum counts in flat on blue side to use orders (100)
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_edgeflat, uves, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;   uves_edgeflat_clean
;   uves_getfil
;   x_edgeflat
;
; REVISION HISTORY:
;   Feb-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro uves_edgeflat, uves, setup, side, CHK=chk, SMSHROW=smshrow, $
                    CLOBBER=clobber, P_NSIG=p_nsig, THRESH=thresh, DEBUG=debug, $
                    INTER=inter, NSIG=nsig, MINFLAT=minflat, NOBASIS=nobasis, $
                    MNSEP=mnsep, $
                    NOSTOP=nostop, NOSORDR=nosordr, IEXTRAP=iextrap, _EXTRA=EXTRA

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_edgeflat, uves, setup, [side], P_NSIG=, THRESH=, /CLOBBER,'
      print, '           /CHK, /INTER, /DEBUG, NSIG=, MINFLAT= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( SIDE ) then side = [1L]
  if not keyword_set( P_NSIG) then p_nsig = 50.
  if not keyword_set( NSIG) then nsig = 5.0
  if not keyword_set( THRESH) then thresh = 100.
  if not keyword_set( BGUESS) then BGUESS = 75L
  if not keyword_set( NCOEFF) then ncoeff = 6L
  if not keyword_set( MINFLAT) then minflat = 100.
  if not keyword_set( NOOVERLAP ) then nooverlap = 4.

  x_psclose
  !p.multi=[0,1,1]

; Loop on chip
;  uves_rdxlog, 'flat_log', 'uves_edgeflat', /ophdr

  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]

      ;; Log
;      uves_rdxlog, 'flat_log', mssg
;      mssg = 'uves_edgeflat: Extrap = '+strjoin(strtrim(extrap,2))
;      uves_rdxlog, 'flat_log', mssg

      ;; Check output file
      idx = where(uves.setup EQ setup and uves.side EQ qq)
      trc_fil = uves_getfil('tflat_str', setup, WCEN=uves[idx[0]].xdangl, $
                            /name, CHKFIL=chkf)
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'uves_edgeflat: Trace exists, moving on..'
          continue
      endif

      ;; Side
      case qq of 
          1: begin
              mssg= 'uves_edgeflat: Tracing BLUE CCD' 
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,0L]
          end
          2: begin
              mssg= 'uves_edgeflat: Tracing RED CCD' 
              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,1L]
              mnsep = 50. / uves[idx[0]].colbin
              tnsig = 1000.
          end
;          3: begin
;              mssg= 'uves_edgeflat: Tracing RED CCD' 
;              if keyword_set(IEXTRAP) then extrap = iextrap else extrap = [1L,1L]
;              gd = where(uves.setup EQ setup AND uves.chip EQ qq $
;                         AND uves.flg_anly NE 0 AND uves.type EQ 'TFLT')
;              TNSIG=100.
;              ;; 2nd order?
;              if uves[gd[0]].cross EQ 'RED' $
;                AND not keyword_set(NOSORDR) then begin
;                  mssg2= 'uves_edgeflat:  Attempting to block 2nd order light'
;                  uves_rdxlog, 'flat_log', mssg2
;                  print, 'uves_edgeflat:  Set /NOSORDR to stop me!'
;                  ;; Use green if it exists
;                  gostr = uves_getfil('ordr_str',setup, CHIP=2)
;                  np = n_elements(gostr[0].lhedg)
;                  MNSEP = gostr[1].lhedg[np/2]-gostr[0].lhedg[np/2] 
;                  SORDR=1
;              endif
;          end
          else: stop
      endcase

      ;; Order structure
      ordr_fil = uves_getfil('ordr_str', setup, WCEN=uves[idx[0]].xdangl, $
                             CHKFIL=chkf, /name)
      if CHKF NE 0 AND $
        not keyword_set( CLOBBER ) then begin
          print, 'uves_edgeflat: Order structure exists. Moving on!', ordr_fil
          continue
      endif


      ;; Read in flat 
      flat = uves_getfil('qtz_fil', setup, WCEN=uves[idx[0]].xdangl, FIL_NM=flat_fil)
      flativar = xmrdfits(flat_fil, 1, /silent)
      sz = size(flat, /dimensions)

      ;; Set binning
      cbin = uves[idx[0]].colbin
      rbin = uves[idx[0]].rowbin

      ;; Output files
      qafil = uves_getfil('qa_trcflat', setup, WCEN=uves[idx[0]].xdangl)
;      logfil = uves_getfil('flat_log')

      x_edgeflat, flat, cbin, trc_fil, ordr_fil, NOOVERLAP=nooverlap, $
        CHK=chk, QAFIL=qafil, FLATIVAR=flativar, SORDR=sordr, MNSEP=mnsep, $
        EXTRAP=extrap, _EXTRA=EXTRA, LOGFIL=logfil, TNSIG=tnsig
      
      ;; Write 
;      mwrfits, trc_str, trc_fil, /create
      if keyword_set(QAFIL) then print, 'uves_edgeflat: Check the QA file == ', qafil

  endfor

  ;; 
  print, 'uves_edgeflat: You may now wish to check the results ' + $
    'with uves_chktrcflat'
  print, 'uves_edgeflat: All done!'

  return
end
