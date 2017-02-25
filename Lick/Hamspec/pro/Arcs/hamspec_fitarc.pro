;+ 
; NAME:
; hamspec_fitarc   
;     Version 1.2
;
; PURPOSE:
;  To identify and centroid arc lines in each order.  There is
;  actually no fitting done at this stage other than to reject bad
;  lines.   The main program calls x_fitarc as its driver.
;  The algorithm x_fitarc does the following:
;
;    1) Input the arc image from hamspec_mkarc 
;    2) Input an archived arc solution appropriate to the chip
;    3) Extract 1D (boxcar) spectra down each order :: extract_arc
;    4) Cross correlate (FFT) against the archived 1D arc spectrum,
;    this gives the order number and the pixel offset
;    5) Automatically identify a set of lines (x_templarc)
;    6) Perform a low order fit to these lines
;    7) Reidentify all lines with this fit and refit 
;    8) Write arc solutions (one per order) to a fits file 
;    9) If the orders extend beyond the archived solution, attempt to
;    extrapolate to the remaining orders.  The idea is to use the
;    known wavelengths from the good orders to extrapolate a guess at
;    the solution.  Then use this guess to peak up on the arc lines in
;    these additional orders.
;
; CALLING SEQUENCE:
;   
;  hamspec_fitarc, hamspec, setup, [obj_id], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   hamspec     -  HIRES structure
;   setup    -  Integer defining setup
;   [obj_id]   -  Object identifier
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit for pre-identified lines
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/hamspec_thar.lst
;   /CHK      - Manually check steps along the way
;   /DEBUG    - Debugging
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   IORDR     - Initial order for analysis
;   /CLOBBER  - Overwrite previous fits
;   SHFTPRM=  - Fit structure for shifting the orders of the arc
;   /SEXTRAP  - Extrapolate wavelength solution in x_fitarc
;   GUESSARC  - Archived file to use in guessing the wavelength
;               solution
;   /NOLOG    - Fit without logarithmic wavelengths
;   IFSIG     - Sigma significance for lines in final fit
;   IPSIG     - Sigma significance for lines in initil fit
;   TWOGUESS= - Archived solution for guess for 'additional' orders
;               not covered by first archived solution
;   /NOFGAUSS - Do NOT use Gaussian fitting for the arc lines
;   /NOSECND  - Do NOT try to to recover other orders
;   STRETCH   - Stretch value, in pixels
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_fitarc, hamspec, 1, 1
;   hamspec_fitarc, hamspec, 1, /clobb, arcfil='Arcs/Arc_B0171.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;
; REVISION HISTORY:
;   20-Feb-2005 Created by JXP (based on mike_fitarc)
;   03-Jan-2008 Made FGAUSS the default
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hamspec_fitarc, hamspec, setup, obj_id, $
                  INTER=inter, LINLIST=linlist, ARCFIL=arcfil, EXTRAP=extrap, $
                  CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                  IORDR=iordr, PINTER=pinter, NOLOG=nolog, NOSECND=nosecnd, $
                  SPINTER=spinter, SINTER=sinter, TWOGUESS=twoguess, $
                  GUESSARC=guessarc, SEXTRAP=sextrap, FWEIGHT=fweight, $
                  NOFGAUSS=nofgauss, STRETCH=stretch, _EXTRA=extra



;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_fitarc, hamspec, setup, obj_id, /INTER, IORDR=, /DEBUG,'
      print, '     /CHK, SIGREJ=, /CLOBBER, /SPINTER, LINLIST=, ARCFIL= '
      print, '     GUESSARC=, STRETCH= [v1.2]'
      return
  endif 

  
;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/old_hires_thar.lst' 
  print, 'hamspec_fitarc: Using line list -- ', linlist

  ;; Fitting procedure -- Gauss is the default
  if not keyword_set(NOFGAUSS) then FGAUSS = 1
  if not keyword_set(STRETCH) then stretch = 0
      
  ;if not keyword_set(ARCFIL) then begin
  ;   ;; Grab all obj indices
  idx = where(hamspec.flg_anly NE 0 AND $
               hamspec.setup EQ setup) 


  ;; Order structure
  ordr_str = hamspec_getfil('ordr_str', setup, fil_nm=ordr_fil)
  ordr_str.flg_anly = 1
  nordr = n_elements(ordr_str)

  ;; Trim Dewar4
  if strtrim(hamspec[idx[0]].ccd,2) EQ 'e2vCCD203-824kx4kthin' then begin
     ordr_str[0:20].flg_anly = 0
     ordr_str[nordr-20:*].flg_anly = 0  ;; This didn't do anything yet..
     TRIM_FFT = 1
     linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ham_dewar4_thar.lst' 
    ; linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/old_hires_thar.lst'

     FEW_LINES=1
  endif
  ;; KLUDGE TO GET ROLLING  29 Aug 2012
  ;ordr_str = ordr_str[15:*]  ;; #35 = Order 98
  ;ordr_str = ordr_str[120:*]  ;; #35 = Order 98
  ;ordr_str = ordr_str[5:*]  ;; #35 = Order 98
      ;;  Would need to modify x_edgtrcflat to enable the following
  ;;  lines of code:  JXP  Oct 2011
  ;gdo = where(ordr_str.flg_anly EQ 1)
  ;ordr_str = ordr_str[gdo]

  ;; LOOP
  
  ;; ROW_SHIFT
  if keyword_set( ROW_SHIFT ) then delvarx, row_shift
     
     ;; ARCS
  arcfil = hamspec_getfil('arc_fil', setup, /name)
     
  ;; IFSIG
  IFSIG = [ [0,1,10], [1,150,10.]] 
  IPSIG = [ [0,1,20], [1,150,20.]] 

  ;; Check for arcfil
  if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
     print, 'hamspec_fitarc: Use hamspec_procarc to create arc file first!!', $
            arcfil
     return
  endif
     
  ;; Outfil
  pos = strpos(arcfil, '.fits')
  out_fil = hamspec_getfil('arc_fit', setup, $
                           /name, CHKFIL=chkfil,  FRAME=frame)
  
  ;; QAFIL
  qafil = hamspec_getfil('qa_arcfit', setup, $
                         /name, CHKFIL=chkfil,  FRAME=frame)
  
  ;; Guessarc
  if not keyword_set( INTER ) then begin
     if not keyword_set( GUESSARC ) then begin
        case strtrim(hamspec[idx[0]].ccd,2) of 
           'e2vCCD203-824kx4kthin': if hamspec[idx[0]].date GT 2456385.4 $
             THEN guessarc = getenv('XIDL_DIR')+'/Lick/Hamspec/calibs/hamspec_dewar4_thar4.idl' $
             ELSE guessarc = getenv('XIDL_DIR')+'/Lick/Hamspec/calibs/hamspec_dewar4_thar1.idl' 
           else: guessarc = getenv('XIDL_DIR')+'/Lick/Hamspec/calibs/hamspec_dewar6.idl'
        endcase
                                ;guessarc = hamspec_arctempl(hamspec[idx])
        print, 'hamspec_fitarc:  Using fil ', guessarc, $
               ' for our guess.'
     endif else begin 
        ;; Add on the path
        stop
        ;guessarc = getenv('HIRES_CALIBS')+'/ARCS/'+guessarc
     endelse
  endif 
     
  ;; Peak width
  head = xheadfits(arcfil)
  sz1 = sxpar(head,'NAXIS2')
                                ;rbin = round(3990./sz1)
  pkwdth = 4L                   ;(6L / rbin) > 2L
     
  ;; First iteration
  print, 'hamspec_fitarc: Using arc file ', arcfil
  rslt = x_fitarc( arcfil, ordr_str, out_fil, $
                   LINLIST=linlist, TRIM_FFT=trim_fft, $ 
                   CHK=chk, CLOBBER=clobber, $
                   INTER=INTER, $
                   PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                   IORDR=iordr,  GUESSARC=guessarc, $
                   OROW_SHIFT=row_shift, IFSIG=ifsig, $
                   IPSIG=ipsig, $
                   MEDWID=2L, $
                   QAFIL=qafil, NOLOG=nolog, $
                   SATUR=45000.*hamspec[idx[0]].gain, $
                   PKWDTH=pkwdth,$
                   STRETCH=stretch, $
                   MXORDER=180L, $
                   FWEIGHT=fweight, $
                   FGAUSS=fgauss, $
                   CONSTR_POS=0.1, $
                   FIT_MSK=fit_msk, $
                   FEW_LINES=few_lines, $
                   NOEXTRAP=(keyword_set(EXTRAP) EQ 0), $
                   _EXTRA=extra)
  if rslt NE -1 then begin  ;; Otherwise file exists

     ;; Trim further (Dewar 4)
     if strtrim(hamspec[idx[0]].ccd,2) EQ 'e2vCCD203-824kx4kthin' then begin
        print, 'hamspec_fitarc: Trimming dewar4 to orders 65-146'
        gdo = where(ordr_str.order LT 65 OR ordr_str.order GT 146)
        ordr_str[gdo].flg_anly = 0
     endif
     ;; Write order structure
     mwrfits, ordr_str, ordr_fil, /create, /silent
     
  endif

; All done
  print, 'hamspec_fitarc: All done!'

  return
end
