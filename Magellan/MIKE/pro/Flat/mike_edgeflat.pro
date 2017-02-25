;+ 
; NAME:
; mike_edgeflat
;     Version 1.1
;
; PURPOSE:
;  To trace the order edges of the individual orders.  And the results are
;  fit (which replaces mike_fittflat).  The fit is a low order fit to the slit
;  width, and then a higher order fit to the order centers.  See
;  x_edgeflat for additional keywords (e.g. EXTRAP, NOOVERLAP, etc.) and discussion.
;  1.  Identify the order edges (interactive is recommended).
;  2.  Performs an order by order tracing of the order edges using
;      trace_crude.
;  3.  Perform (iteratively) a PCA analysis on the coefficients of the
;      individual traces.
;  4.  Create and write to disk a structure summarizing the fits.
;  5.  Do a 2-part fit that replaces the fittflat call.
;  See x_edgeflat for a full set of KEYWORDS
;
;
; CALLING SEQUENCE:
;   
;  mike_edgeflat, mike, setup, [side], /CHK, 
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to mike_fittflat to create a 2D
;  solution.  Typical name:  'Flats/TStr_B_01.fits'
;
; OPTIONAL KEYWORDS:
;  /CHK -- Check the order edges interactively (similar to INTER)
;  /INTER -- Determine the order edges interactively
;  /CLOBBER -- Overwrite the output structure
;  NCOEFF  -- Number of Legendre coefficients for tracing [default: 6]
;  NOOVERLAP -- If orders are closer than NOOVERLAP, then they are not
;               included in the PCA fit but are extrapolated from the
;               PCA analysis.
;  KEEP_FRAC -- This value sets the minimum fractional amount that an
;               order can be on the CCD to be included in the PCA
;               [Default: 0.9]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_edgeflat, mike, 1, 1, /clobber
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_edgeflat
;   mike_getfil
;   x_trcflat_edges
;
; REVISION HISTORY:
;   1-Nov-2004  Written by SB, replaces mike_trcflat and mike_fittflat
;                QA file is still in mike_trcflat
;   1-May-2005  Consumed by x_edgeflat (JXP)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro mike_edgeflat, mike, setup, side, CHK=chk, NOOVERLAP=nooverlap, $
                   CLOBBER=clobber, _EXTRA=EXTRA

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_edgeflat, mike, setup, [side], /CLOBBER,'
      print, '           /CHK, /DEBUG, MINFLAT=, _EXTRA [v2.0]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
;  if not keyword_set( P_NSIG) then p_nsig = 50.
;  if not keyword_set( NSIG) then nsig = 5.0
;  if not keyword_set( THRESH) then thresh = 100.
;  if not keyword_set( BGUESS) then BGUESS = 75L
;  if not keyword_set( NCOEFF) then ncoeff = 6L
;  if not keyword_set( MINFLAT) then minflat = 100.

  x_psclose
  !p.multi=[0,1,1]

; Loop on side

  for iside=0L,n_elements(side)-1 do begin
      qq = side[iside]
      ;; SIDE
      if qq EQ 1 then print, 'mike_edgeflat: Tracing BLUE trace flat' $
      else print, 'mike_edgeflat: Tracing RED trace flat'


      ;; Check output file
      trc_fil = mike_getfil('tflat_str', setup, SIDE=qq, /name, CHKFIL=chkf)
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'mike_edgeflat: Trace exists, moving on..'
          continue
      endif

      ;; OVERLAP
      idx = where(mike.setup EQ setup AND mike.side EQ qq AND $
                  mike.type EQ 'TFLT')
      cbin = mike[idx[0]].colbin
      if not keyword_set( NOOVERLAP) then nooverlap = 2.1 * 3 / float(cbin)

      ;; Read in flat 
      flat_file = mike_getfil('tflat_fil', setup, SIDE=qq, /name)
      flat = xmrdfits(flat_file, /silent)
      flativar = xmrdfits(flat_file,1, /silent)

      ;; Set binning
      sz = size(flat, /dimensions)
      cbin = round(2048. / sz[0])
      rbin = round(4096. / sz[1])

      ;; Files
      qafil = mike_getfil('qa_trcflat', setup, SIDE=qq)
      trc_fil = mike_getfil('tflat_str', setup, SIDE=qq, /name, CHKFIL=chkf)
      ordr_fil = mike_getfil('ordr_str', setup, SIDE=qq, CHKFIL=chkf, /name)
      if CHKF NE 0 AND $
        not keyword_set( CLOBBER ) then begin
          print, 'mike_edgeflat: Order structure exists. Moving on!', ordr_fil
          continue
      endif
      
;      start_guess_order = qq EQ 2 ? 80 : BGUESS

      ;; Call
      x_edgeflat, flat, cbin, trc_fil, ordr_fil, ZERO_OUT=(qq EQ 2), $
        QAFIL=qafil, _EXTRA=EXTRA, FLATIVAR=flativar, NOOVERLAP=NOOVERLAP
      
  endfor
  ;; 
  print, 'mike_edgeflat: You may now wish to check the results with mike_chktrcflat'
  print, 'mike_edgeflat: All done!'

  return
end
