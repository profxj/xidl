;+ 
; NAME:
; esi_echtrcstd   
;     Version 1.1
;
; PURPOSE:
;    Trace a standard star in each ordrer
;
; CALLING SEQUENCE:
;   
;  esi_echtrcstd, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  Image with trcstdered light removed
;
; OPTIONAL KEYWORDS:
;  /NOFND   - Do not repeat step to find object
;  /NOSKY   - Do not repeat sky subtraction
;  /CHK     - Show the final trace
;  MXERR=   - Maximum error in trace centering to include 
;                 (default: 0.3 pix)
;  OFF=     - Offset used in find object routine
;  /CUAR    - Data has only CuAr lamps
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Only setup for 1x1 binning (some hard numbers)
;
; EXAMPLES:
;   esi_echtrcstd, esi, 1.0, /CHK
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtrcstd, esi, slit, NOFND=nofnd, NOSKY=nosky, CHK=chk, MXERR=mxerr,$
                              OFF=off, CUAR=cuar

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtrcstd, esi, slit, /NOFND, /NOSKY, /CHK, OFF=, MXERR= '
      prtin, '     /CUAR [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then ncoll = 10L
  if not keyword_set( MXERR ) then mxerr = 0.3

;;;;;;
;  Find standard star
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.slit EQ slit AND esi.type EQ 'STD', nindx)
  case nindx of 
      0 : begin
          print, 'esi_echtrcstrd: No standard star images! Returning..'
          return
      end 
      1 : print, 'esi_echtrcstd: Tracing standard star image --- ', $
        esi[indx].img_root
      else : begin
          print, 'esi_echtrcstd: Warning -- Multiple standard star images'
          indx = indx[0]
          print, 'esi_echtrcstd: Taking first one ', esi[indx].img_root
      end
  endcase
          
;;;;;;;;;
; Process
  esi_echproc, esi, indx, SUBSCAT=subscat

;;;;;;;;;
; FIND Obj + Create Obj Structure
  if not keyword_set( NOFND ) then $
    esi_echfndobj, esi, indx, /STD, SCICLM=off

;;;;;;;;;
; Sky subtract
  if not keyword_set( NOSKY ) then begin
      if keyword_set( CUAR ) then ordr = [0L, 8L] else ordr = [0L,9L]
      esi_echskysub, esi, indx, /STD, ORDR=ordr
  endif
  
;;;;;;;;;
; Open Stuff

  ;; Open Slit file
  c_s = esi_slitnm(esi[indx[0]].slit)
  sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
  if x_chkfil(sedg_fil+'*') EQ 0 then begin
      print, 'esi_echtrcstd: Slit edge file not found ', sedg_fil
      return
  endif
  print, 'esi_echtrcstd: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)

  ;; Open Image, Variance
  if x_chkfil(esi[indx].img_final+'*') EQ 0 then begin
      print, 'esi_echtrcstd: Slit edge file not found ', esi[indx].img_final
      return
  endif
  print, 'esi_echtrcstd: Opening images...'
  img = xmrdfits(esi[indx].img_final, 2, /silent)
  sz_img = size(img, /dimensions)
  var = xmrdfits(esi[indx].img_final, 1, /silent)

  ;; OBJ
  objfil = esi[indx].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'esi_echtrcstd: Object file not found ', objfil
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  subimg = fltarr(21,ncoll)
  subvar = fltarr(21,ncoll)
  nfit = sz_img[1]/ncoll
  xfit = dblarr(nfit)
  xsig = dblarr(nfit)
  yfit = dblarr(nfit)
  
  for qq=ordr[0],ordr[1] do begin
      print, 'esi_echtrcstd: Tracing order ', string(15L-qq, FORMAT='(i3)')
      ;; Zero out subimg
      subimg[*] = 0.
      svoff = 0.
      ;; Trace
      rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
      ;; Order 14 offset!
      if qq EQ 1 AND not keyword_set( NOFIX14 ) then $
        rnd_trc[0:500] = rnd_trc[0:500] - 7L
      ;; CHK
;      if keyword_set( CHK ) then begin
;          trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
;          tmp = img
;          tmp[trc_msk] = -1
;          xatv, tmp, /block
;          stop
;      endif
      case qq of
          0: begin
              jstrt =1500L
              jend = 3790L
              step = ncoll
          end
          9: begin
              jstrt = ncoll + 3L  
              jend = 2170
              step = ncoll
          end
          else: begin
              jstrt = sz_img[1]
              jend = ncoll + 3L  ;; First 3 rows are crummy
              step = -ncoll
          end
      endcase
      gdfit = 0L
      ;;;;;;;;;;;;;;;;;;;;
      ;; LOOP
      for j=jstrt,jend,step do begin
          ;; Create Subimg
          jp = 0L
          mn = min(rnd_trc[j-ncoll + lindgen(ncoll)]-10L)
          if mn LT 0L then begin
              print, 'esi_echtrcstd: Standard not well centered in slit'
              print, 'esi_echtrcstd: Consider using alternative standard'
              stop
          endif
          for jsub = (j-ncoll),j-1 do begin
              subimg[*,jp] = img[rnd_trc[jsub]-10L:rnd_trc[jsub]+10L,jsub]
              subvar[*,jp] = var[rnd_trc[jsub]-10L:rnd_trc[jsub]+10L,jsub]
              jp = jp+1
          endfor
          ;; Median
          mdn = djs_median(subimg,2)
          mdn_var = djs_median(subvar,2)
          ivar = 1./(mdn_var > 0.)
          xcen = 10.d + svoff
          ;; Find Centroid
          for k=0,19 do $
            xcen = trace_fweight(mdn, xcen, 0L, radius=3., xerr=xerr, invvar=ivar)
          ;; Keep the good points
          if xerr LT mxerr AND xerr GT 0.00000001 then begin
              xsig[gdfit] = xerr
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - 10.
              svoff = xcen-10.
              gdfit = gdfit + 1L
          endif
;          plot, mdn
;          oplot, [xcen,xcen], [-1000,1000]
;          print, xcen, xerr
;          stop
      endfor

      ;; CHK
      if gdfit LT 5 then stop

      ;; FIT
      trc_fit = x_setfitstrct(NITER=2L, NORD=9L, FLGREJ=1L, HSIG=5., LSIG=5., $
                              FUNC='POLY')
      if qq EQ 9 then trc_fit.nord = 5L
      new_trc = x_fitrej(yfit[0:gdfit-1], xfit[0:gdfit-1], SIG=xsig[0:gdfit-1], $
                         FITSTR=trc_fit)
      
      print, 'esi_echtrcstd: RMS = ', trc_fit.rms
      ;; objstr
      objstr[qq].trace[0:sz_img[1]-1] = x_calcfit(findgen(sz_img[1]), $
                                                  FITSTR=trc_fit)

  endfor

  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      for qq=ordr[0],ordr[1] do begin
          trc_msk = round(objstr[qq].trace[0:sz_img[1]-1]) $
            + lindgen(sz_img[1])*sz_img[0]
          tmp[trc_msk] = -1000
      endfor
      xatv, tmp, /block
  endif

; OUTPUT

  ;; STD
  mwrfits, objstr, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
  ;; TRACE
  std_trc = 'Extract/STD_ECH'+c_s+'_TRC.fits'
  mwrfits, objstr[0:9].trace, std_trc, /create, /silent
  

  return
end
              
      
      
