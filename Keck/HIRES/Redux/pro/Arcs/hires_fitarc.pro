;+ 
; NAME:
; hires_fitarc   
;     Version 1.1
;
; PURPOSE:
;  To identify and centroid arc lines in each order.  There is
;  actually no fitting done at this stage other than to reject bad
;  lines.   The main program calls x_fitarc as its driver.
;  The algorithm x_fitarc does the following:
;
;    1) Input the arc image from hires_mkarc 
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
;  hires_fitarc, hires, setup, [obj_id, chip], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   hires     -  HIRES structure
;   setup    -  Integer defining setup
;   [obj_id]   -  Object identifier
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit for pre-identified lines
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/hires_thar.lst
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
;   hires_fitarc, hires, 1, 1
;   hires_fitarc, hires, 1, /clobb, arcfil='Arcs/Arc_B0171.fits'
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
pro hires_fitarc, hires, setup, obj_id, chip, $
                  INTER=inter, LINLIST=linlist, ARCFIL=arcfil, EXTRAP=extrap, $
                  CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                  IORDR=iordr, PINTER=pinter, NOLOG=nolog, NOSECND=nosecnd, $
                  SPINTER=spinter, SINTER=sinter, TWOGUESS=twoguess, $
                  GUESSARC=guessarc, SEXTRAP=sextrap, FWEIGHT=fweight, $
                  NOFGAUSS=nofgauss, _EXTRA=extra



;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_fitarc, hires, setup, obj_id, [chip], /INTER, IORDR=, /DEBUG,'
      print, '     /CHK, SIGREJ=, /CLOBBER, /SPINTER, LINLIST=, ARCFIL= '
      print, '     GUESSARC= [v1.1]'
      return
  endif 

  
;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 
  print, 'hires_fitarc: Using line list -- ', linlist

  if keyword_set(ARCFIL) then begin
      ;; Grab Chip from name
      case strmid(arcfil,9,1) of
          'S': chip = -1L
          'B': chip = 1L
          'G': chip = 2L
          'R': chip = 3L
          else: stop
      endcase
  endif

  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  
  ;; Fitting procedure -- Gauss is the default
  if not keyword_set(NOFGAUSS) then FGAUSS = 1
      
  ;; Loop on chip
  for kk=0L,n_elements(chip)-1 do begin
      qq = chip[kk]
      ;; CHIP
      case qq of 
         -1: begin
            print, 'hires_fitarc: Fitting Single arc'
            ccdsz = [2048L,2048]
         end
         1: print, 'hires_fitarc: Fitting BLUE arc' 
         2: print, 'hires_fitarc: Fitting GREEN arc' 
         3: print, 'hires_fitarc: Fitting RED arc' 
         else: stop
      endcase

      if not keyword_set(ARCFIL) then begin
          ;; Grab all obj indices
          indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                       hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                       (strtrim(hires.type,2) EQ 'OBJ' OR $
                        strtrim(hires.type,2) EQ 'STD'), nindx)
          if nindx EQ 0 then begin
              print, 'hires_fitarc: No Obj found!  Skipping..' 
              continue
          endif
          
          ;; Unique
          arcfil_all= hires[indx].arc_fil
          asort = sort(arcfil_all)
          auniq = uniq(arcfil_all[asort])
      endif else begin ;; Single frame
          arcfil_all = [arcfil]
          auniq = [0L]
          asort = [0L]
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))
          indx = where(hires.frame EQ frame AND hires.chip EQ chip)
      endelse


      ;; Order structure
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip, fil_nm=ordr_fil)
      ordr_str.flg_anly = 1
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
          if hires[idx].cross EQ 'UV' then begin
              IFSIG = [ [0,1,10], [1,150,10.]] 
              IPSIG = [ [0,1,20], [1,150,20.]] 
          endif else begin
              IFSIG = [ [0,1,20], [1,46,5.], [46, 150, 20.]] 
              IPSIG = [ [0,1,40], [1,46,5.], [46, 150, 40.]] 
          endelse

          if hires[idx].lampfil EQ 'ug5' then begin
              IFSIG = [ [0,1,10], [1,150,10.]] 
              IPSIG = [ [0,1,20], [1,150,20.]] 
          endif

          ;; Check for arcfil
          if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
              print, 'hires_fitarc: Use hires_procarc to create arc file first!!', $
                arcfil
              continue
          endif

          ;; Outfil
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))
          out_fil = hires_getfil('arc_fit', setup, CHIP=chip, $
                                 /name, CHKFIL=chkfil,  FRAME=frame)

          ;; QAFIL
          qafil = hires_getfil('qa_arcfit', setup, CHIP=chip, $
                                 /name, CHKFIL=chkfil,  FRAME=frame)

          ;; Guessarc
          if not keyword_set( INTER ) then begin
              if not keyword_set( GUESSARC ) then begin
                  guessarc = hires_arctempl(hires[idx], SINGLE=(CHIP EQ -1))
                  print, 'hires_fitarc:  Using fil ', guessarc, $
                         ' for our guess.'
              endif else begin 
                  ;; Add on the path
                  guessarc = getenv('HIRES_CALIBS')+'/ARCS/'+guessarc
              endelse
          endif 

          ;; Peak width
          head = xheadfits(arcfil)
          sz1 = sxpar(head,'NAXIS2')
          if qq NE -1 then rbin = round(3990./sz1) else rbin = round(2048./sz1)
          pkwdth = (6L / rbin) > 2L
          ;; First iteration
          rslt = x_fitarc( arcfil, ordr_str, out_fil, $
                           INTER=inter, LINLIST=linlist,$ 
                           CHK=chk, CLOBBER=clobber, /HIRES, $
                           PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                           IORDR=iordr, GUESSARC=guessarc, $
                           OROW_SHIFT=row_shift, IFSIG=ifsig, $
                           IPSIG=ipsig, $
                           QAFIL=qafil, NOLOG=nolog, $
                           SATUR=45000.*hires[idx].gain, $
                           PKWDTH=pkwdth,$
                           CCDSZ=ccdsz, $
                           FWEIGHT=fweight, $
                           FGAUSS=fgauss, $
                           CONSTR_POS=0.15, $
                           FIT_MSK=fit_msk, $
                           NOEXTRAP=(keyword_set(EXTRAP) EQ 0), $
                           _EXTRA=extra)
          if rslt NE -1 then begin  ;; Otherwise file exists
              ;; Write order structure
              mwrfits, ordr_str, ordr_fil, /create, /silent
              
              ;; Second iteration
              nofit = where(FIT_MSK EQ 0B, nno)
              if nno NE 0 then begin
                  if keyword_set(NOSECND) then begin
                      ordr_str[nofit].flg_anly = -1
                  endif else begin
                      ;; Deal with 'extra' orders
                      if nno GT 1 then begin
                          ;; Look for a gap between masked orders
                          od = ordr_str[nofit].order
                          sep = od - shift(od,-1)
                          sep = sep[0:nno-2]
                          isep = where(sep NE (-1), nsep)
                          if nsep NE 0 then nloop = 2 else begin
                              nloop = 1
                              isep  = nno-1
                          endelse
                      endif else begin
                          nloop = 1
                          isep = 0
                      endelse
                      ;; Loop on the sets of 'extra' orders
                      for ll=0,nloop-1 do begin
;                         if ll EQ 1 then stop
                          if ll EQ 0 then inofit = nofit[0:isep] $
                          else inofit = nofit[isep+1:*]
                          if not keyword_set(TWOGUESS) then begin
                              guessarc = hires_arctempl(hires[idx], $
                                                        ordr_str[inofit].order, $
                                                       SINGLE=(CHIP EQ -1)) 
                          endif else guessarc = TWOGUESS
                          print, 'hires_fitarc:  Now using ', guessarc
;                          print, 'hires_fitarc:  Shifting by ', row_shift
                          iordr = min(ordr_str[inofit].order, MAX=fordr)
                          rslt = x_fitarc( arcfil, ordr_str, out_fil, $
                                           INTER=sinter, LINLIST=linlist,$ 
                                           CHK=chk, CLOBBER=clobber, /HIRES, $
                                           PINTER=spinter, DEBUG=debug, $
                                           SIGREJ=sigrej, $;IROW_SHIFT=row_shift, $
                                           IORDR=iordr, FORDR=fordr, $
                                           CONSTR_POS=0.15, $
                                           GUESSARC=guessarc, PKWDTH=pkwdth, $
                                           QAFIL=qafil, NOLOG=nolog, $
                                           IFSIG=ifsig, IPSIG=ipsig, $
                                           SATUR=45000., CLIP=0.3, $
                                           NOEXTRAP=(keyword_set(SEXTRAP) EQ 0) $
                                           ,_EXTRA=extra)
;                                           IROW_SHIFT=1, MXSHFT=150L, $
                      endfor
                  endelse
                  ;; Write order structure
                  mwrfits, ordr_str, ordr_fil, /create, /silent
              endif
;                          
              if size(rslt,/tname) NE 'STRING' then stop
          endif
      endfor
  endfor

; All done
  print, 'hires_fitarc: All done!'

  return
end
