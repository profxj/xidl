;+ 
; NAME:
; uves_fitarc   
;     Version 1.1
;
; PURPOSE:
;  To identify and centroid arc lines in each order.  There is
;  actually no fitting done at this stage other than to reject bad
;  lines.   The main program calls x_fitarc as its driver.
;  The algorithm x_fitarc does the following:
;
;    1) Input the arc image from uves_mkarc 
;    2) Input an archived arc solution appropriate to the side
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
;  uves_fitarc, uves, setup, obj_id, [side], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   uves     -  HIRES structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit for pre-identified lines
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/uves_thar.lst
;   /CHK      - Manually check steps along the way
;   /DEBUG    - Debugging
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   IORDR     - Initial order for analysis
;   /CLOBBER  - Overwrite previous fits
;   SHFTPRM=  - Fit structure for shifting the orders of the arc
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_fitarc, uves, 1, 1
;   uves_fitarc, uves, 1, /clobb, arcfil='Arcs/Arc_B0171.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;
; REVISION HISTORY:
;   20-Feb-2005 Created by JXP (based on mike_fitarc)
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uves_fitarc, uves, setup, side, $
                  INTER=inter, LINLIST=linlist, ARCFIL=arcfil, EXTRAP=extrap, $
                  CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                  IORDR=iordr, PINTER=pinter, NOLOG=nolog, NOSECND=nosecnd, $
                  SPINTER=spinter, SINTER=sinter, TWOGUESS=twoguess, $
                  GUESSARC=guessarc, SEXTRAP=sextrap, _EXTRA=extra



;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_fitarc, uves, setup, [side], /INTER, IORDR=, /DEBUG,'
      print, '     /CHK, SIGREJ=, /CLOBBER, /SPINTER, LINLIST=, ARCFIL= '
      print, '     GUESSARC= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 
  if keyword_set(ARCFIL) then begin
      if not keyword_set(SETUP) then begin
          prs = strsplit(arcfil, '_')
          ;; Grab setup from name
          setup = long(prs[2])
      endif
  endif

  if not keyword_set( side ) then side = [1L]
      
  ;; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]
      ;; SIDE
      case qq of 
          1: print, 'uves_fitarc: Fitting BLUE arc' 
          2: print, 'uves_fitarc: Fitting RED arc' 
          else: stop
      endcase

      if not keyword_set(ARCFIL) then begin
          ;; Grab all obj indices
          indx = where(uves.flg_anly NE 0 AND uves.side EQ qq AND $
                       uves.setup EQ setup AND $
                       strtrim(uves.type,2) EQ 'ARC', nindx)
          if nindx EQ 0 then begin
              print, 'uves_fitarc: No Arc found!  Skipping..' 
              continue
          endif
          
          wcen = uves[indx[0]].xdangl

          ;; Unique
          arcfil_all= [uves_getfil('arc_fil', setup, WCEN=wcen, /name)]
          asort = sort(arcfil_all)
          auniq = uniq(arcfil_all[asort])
      endif else begin ;; Single frame
          arcfil_all = [arcfil]
          auniq = [0L]
          asort = [0L]
          wcen = uves_getwcen(arcfil,FRAME=frame, /AFIL)
          indx = where(uves.setup EQ setup AND uves.side EQ qq AND $
                      uves.type EQ 'ARC')
      endelse

      ;; Order structure
      ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen, fil_nm=ordr_fil)
      ordr_str.flg_anly = 1

      ;; DEBUGGING
;      ordr_str = ordr_str[0:31]

      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin
 
          ;; ROW_SHIFT
          if keyword_set( ROW_SHIFT ) then delvarx, row_shift

          ;; ARCS
          arcfil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; IFSIG
;          if uves[idx].cross EQ 'UV' then begin
              IFSIG = [ [0,1,10], [1,200,5.]] 
              IPSIG = [ [0,1,20], [1,200,20.]] 
;          endif else begin
;              IFSIG = [ [0,1,20], [1,46,5.], [46, 150, 20.]] 
;              IPSIG = [ [0,1,40], [1,46,15.], [46, 150, 40.]] 
;          endelse

;          if uves[idx].lampfil EQ 'ug5' then begin
;              IFSIG = [ [0,1,10], [1,150,10.]] 
;              IPSIG = [ [0,1,20], [1,150,20.]] 
;          endif

          ;; Check for arcfil
          if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
              print, 'uves_fitarc: Use uves_procarc to create arc file first!!', $
                arcfil
              continue
          endif

          ;; Outfil
          out_fil = uves_getfil('arc_fit', setup, WCEN=wcen, FRAME=frame, $
                                 /name, CHKFIL=chkfil)

          ;; QAFIL
          qafil = uves_getfil('qa_arcfit', setup, WCEN=wcen,  FRAME=frame, /name)

          ;; Guessarc
          if not keyword_set( INTER ) then begin
              if not keyword_set( GUESSARC ) then begin
                  case uves[idx].side of
                      1: cside = 'blue'
                      2: cside = 'red'
                      else: stop
                  endcase
                  guessarc = getenv('XIDL_DIR')+ $
                             '/VLT/UVES/pro/Arcs/Templates/uves_arc_'+$
                             cside+'_'+strtrim(round(uves[idx].xdangl),2)+'.idl'
              endif
              print, 'uves_fitarc:  Using fil ', guessarc, $
                     ' for our guess.'
          endif 

          ;; Peak width
          head = xheadfits(arcfil)
          sz1 = sxpar(head,'NAXIS2')
          if side EQ 1 then begin
              rbin = round(2940./sz1)
              ccdsz = [2048L, sz1]
              ostep = -1
          endif else begin
              rbin = round(4050./sz1)
              ccdsz = [2048L, sz1]
              ostep = -1
          endelse
          pkwdth = 6L / rbin
          
          ;; First iteration
          rslt = x_fitarc( arcfil, ordr_str, out_fil, $
                           INTER=inter, LINLIST=linlist,$ 
                           CHK=chk, CLOBBER=clobber, CCDSZ=ccdsz, $
                           PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                           IORDR=iordr, GUESSARC=guessarc, $
                           OROW_SHIFT=row_shift, IFSIG=ifsig, $
                           IPSIG=ipsig, OSTEP=ostep, $
                           QAFIL=qafil, NOLOG=nolog, $
                           SATUR=30000.*uves[idx].gain, $
                           PKWDTH=pkwdth,$
                           CONSTR_POS=0.15, $
                           FIT_MSK=fit_msk, $
                           NOEXTRAP=(keyword_set(EXTRAP) EQ 0), $
                           _EXTRA=extra)

          ;; 
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
                      for ll=0,-1 do begin
;                      for ll=0,nloop-1 do begin
                          if ll EQ 0 then inofit = nofit[0:isep] $
                          else inofit = nofit[isep+1:*]
                          if not keyword_set(TWOGUESS) then begin
                              guessarc = uves_arctempl(uves[idx], $
                                                        ordr_str[inofit].order) 
                          endif else guessarc = TWOGUESS
                          print, 'uves_fitarc:  Now using ', guessarc
;                          print, 'uves_fitarc:  Shifting by ', row_shift
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
                                           NOEXTRAP=(keyword_set(SEXTRAP) EQ 0))
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
  print, 'uves_fitarc: All done!'

  return
end
