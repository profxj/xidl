;+ 
; NAME:
; apf_fitarc   
;     Version 1.1
;
; PURPOSE:
;  To identify and centroid arc lines in each order.  There is
;  actually no fitting done at this stage other than to reject bad
;  lines.   The main program calls x_fitarc as its driver.
;  The algorithm x_fitarc does the following:
;
;    1) Input the arc image from apf_mkarc 
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
;  apf_fitarc, apf, setup, [obj_id], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   apf     -  HIRES structure
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
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/apf_thar.lst
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_fitarc, apf, 1, 1
;   apf_fitarc, apf, 1, /clobb, arcfil='Arcs/Arc_B0171.fits'
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
pro apf_fitarc, apf, setup, obj_id, $
                  INTER=inter, LINLIST=linlist, ARCFIL=arcfil, EXTRAP=extrap, $
                  CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                  IORDR=iordr, PINTER=pinter, NOLOG=nolog, NOSECND=nosecnd, $
                  SPINTER=spinter, SINTER=sinter, TWOGUESS=twoguess, $
                  GUESSARC=guessarc, SEXTRAP=sextrap, FWEIGHT=fweight, $
                  NOFGAUSS=nofgauss, _EXTRA=extra



;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_fitarc, apf, setup, obj_id, /INTER, IORDR=, /DEBUG,'
      print, '     /CHK, SIGREJ=, /CLOBBER, /SPINTER, LINLIST=, ARCFIL= '
      print, '     GUESSARC= [v1.1]'
      return
  endif 

  
;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/old_hires_thar.lst' 
  print, 'apf_fitarc: Using line list -- ', linlist

  ;; Fitting procedure -- Gauss is the default
  if not keyword_set(NOFGAUSS) then FGAUSS = 1
      
  if not keyword_set(ARCFIL) then begin
     ;; Grab all obj indices
     indx = where(apf.flg_anly NE 0 AND $
                  apf.obj_id EQ obj_id AND apf.setup EQ setup AND $
                  (strtrim(apf.type,2) EQ 'OBJ' OR $
                   strtrim(apf.type,2) EQ 'STD'), nindx)
     if nindx EQ 0 then begin
        print, 'apf_fitarc: No Obj found!  Skipping..' 
        return
     endif
     
     ;; Unique
     arcfil_all= apf[indx].arc_fil
     asort = sort(arcfil_all)
     auniq = uniq(arcfil_all[asort])
  endif else begin ;; Single frame
     arcfil_all = [arcfil]
     auniq = [0L]
     asort = [0L]
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-4,4))
     indx = where(apf.frame EQ frame)
  endelse


  ;; Order structure
  ordr_str = apf_getfil('ordr_str', setup, fil_nm=ordr_fil)
  ordr_str.flg_anly = 1
  ;; KLUDGE TO GET ROLLING  29 Aug 2012
  ordr_str = ordr_str[15:17]
      ;;  Would need to modify x_edgtrcflat to enable the following
  ;;  lines of code:  JXP  Oct 2011
;      gdo = where(ordr_str.flg_anly EQ 1)
;      ordr_str = ordr_str[gdo]

  ;; LOOP
  for mm=0L,n_elements(auniq)-1 do begin
     
     ;; ROW_SHIFT
     if keyword_set( ROW_SHIFT ) then delvarx, row_shift
     
     ;; ARCS
     arcfil = arcfil_all[asort[auniq[mm]]]
     idx = indx[asort[auniq[mm]]]
     
     ;; IFSIG
     IFSIG = [ [0,1,10], [1,150,10.]] 
     IPSIG = [ [0,1,20], [1,150,20.]] 

     ;; Check for arcfil
     if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
        print, 'apf_fitarc: Use apf_procarc to create arc file first!!', $
               arcfil
        continue
     endif
     
     ;; Outfil
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-4,4))
     out_fil = apf_getfil('arc_fit', setup, $
                          /name, CHKFIL=chkfil,  FRAME=frame)
     
          ;; QAFIL
     qafil = apf_getfil('qa_arcfit', setup, $
                        /name, CHKFIL=chkfil,  FRAME=frame)
     
     ;; Guessarc
     if not keyword_set( INTER ) then begin
        if not keyword_set( GUESSARC ) then begin
           guessarc = getenv('HIRES_CALIBS')+'/ARCS/hires_tmplR2x1G2.idl'
           ;guessarc = apf_arctempl(apf[idx])
           print, 'apf_fitarc:  Using fil ', guessarc, $
                  ' for our guess.'
        endif else begin 
           ;; Add on the path
           guessarc = getenv('HIRES_CALIBS')+'/ARCS/'+guessarc
        endelse
     endif 
     
     ;; Peak width
     head = xheadfits(arcfil)
     sz1 = sxpar(head,'NAXIS2')
     rbin = round(3990./sz1)
     pkwdth = (6L / rbin) > 2L
     
     ;; First iteration
     rslt = x_fitarc( arcfil, ordr_str, out_fil, $
                      LINLIST=linlist,$ 
                      CHK=chk, CLOBBER=clobber, $
                      /INTER, $
                      PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                      IORDR=iordr,$; GUESSARC=guessarc, $
                      OROW_SHIFT=row_shift, IFSIG=ifsig, $
                      IPSIG=ipsig, $
                      MEDWID=2L, $
                      QAFIL=qafil, NOLOG=nolog, $
                      SATUR=45000.*apf[idx].gain, $
                      PKWDTH=pkwdth,$
                      FWEIGHT=fweight, $
                      FGAUSS=fgauss, $
                      CONSTR_POS=0.15, $
                      FIT_MSK=fit_msk, $
                      NOEXTRAP=(keyword_set(EXTRAP) EQ 0), $
                      _EXTRA=extra)
     if rslt NE -1 then begin  ;; Otherwise file exists
        ;; Write order structure
        mwrfits, ordr_str, ordr_fil, /create, /silent
        
     endif
  endfor

; All done
  print, 'apf_fitarc: All done!'

  return
end
