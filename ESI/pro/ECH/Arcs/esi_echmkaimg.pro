;+ 
; NAME:
; esi_echmkaimg   
;     Version 1.1
;
; PURPOSE:
;    Creates a wavelength image given the arc trace
;
; CALLING SEQUENCE:
;   
;  esi_echmkaimg, esi, slit, TRCSTR=, NX=, NY=, /CHK
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One wavelength map for the given slit size 
;     (e.g. Arcs/Arc/ECH_##IMG.fits)
;
; OPTIONAL KEYWORDS:
;   /CHK     - Manually check steps along the way
;   TRCSTR=  - Trace structure for fitting
;   NX, NY=  - Order of x,y 2D poly fit (default: 9,3)
;   /CUAR    - CuAr lamps only!
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Only setup for 1x1 binning (some hard numbers)
;  Zeros out image for rows [0-1500] in order 15; [2170-end] in order6
;
; EXAMPLES:
;   esi_echmkaimg, esi, 1.0, /CHK
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmkaimg, esi, slit, TRCSTR=trcstr, NX=nx, NY=ny, CHK=chk,$
                   CUAR=cuar

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echmkaimg, esi, slit, /CHK, NX=, NY=, TRCSTR=, /CUAR [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SZ_ARC ) then sz_arc = [2048L,4096L]
  if not keyword_set( MAP_FIL ) then map_fil = 'Maps/ECH_map.fits'

  c_s = esi_slitnm(slit)

; Open Map
  print, 'esi_echmkaimg: Reading map file: ', map_fil
  if x_chkfil(map_fil+'*') EQ 0 then begin
      print, 'esi_echmkaimg: Map file doesnt exist', map_fil, ' Returning...'
      return
  endif
  map = xmrdfits(map_fil, /silent)
  ;; Transpose the map
  map = transpose(map)

; Final Arc
  fin_arc = dblarr(sz_arc[0],sz_arc[1]) + 1.d

; Loop
  if not keyword_set( CUAR ) then qend = 9L else qend = 8L
  for qq=0L,qend do begin
      ;; Grab Arc Trace
      ordr = 15L - qq
      if ordr LT 10 then cordr = '0'+string(ordr, FORMAT='(i1)') $
      else cordr = string(ordr, FORMAT='(i2)')
      trcfil = 'Arcs/TRC/ArcECH_'+c_s+'trc'+cordr+'.fits'
      a = findfile(trcfil, count=na)
      if na EQ 0 then begin
          print, 'esi_echmkaimg: Trace ',trcfil,' doesnt exist. Run esi_echtrcarc!'
          return
      endif
      print, '------------------------------------------------'
      print, 'esi_echmkaimg: Order ',cordr
      print, 'esi_echmkaimg: Reading trace structure: ', trcfil
      trcstr = xmrdfits(trcfil, 1, /silent)
      if x_chkfil(trcfil+'*') EQ 0 then begin
          print, 'esi_echmkaimg: Trace file doesnt exist', trcfil, ' Returning...'
          return
      endif
      ntrc = n_elements(trcstr.wav)
      sz_trc = size(trcstr.xcen, /dimensions)
  
      ;; Setup Fit Structure
      fitstr = { fit2dstrct }
      fitstr.func = 'POLY'
      if keyword_set( NX ) then fitstr.nx = nx else begin
          if keyword_set( CUAR ) then begin
              case qq of
                  8L: fitstr.nx = 4
                  else: fitstr.nx = 5
              endcase
          endif else fitstr.nx = 9  ;; Default
      endelse
      if keyword_set( NY ) then fitstr.ny = ny else begin
          if keyword_set( CUAR ) then fitstr.ny = 2 $
          else fitstr.ny = 3  ;; Default
      endelse
      fitstr.niter = 1
      fitstr.lsig = 3.
      fitstr.hsig = 3.
      fitstr.flg_rej = 1

      ;; Setup Arrays
      ycen = lindgen(sz_trc[0]) # replicate(1., ntrc) 
      gd_cen = where(trcstr.xerr LE 0.2, ngd)  ;; Changed from 0.1
      xydat = dblarr(ngd,2)
      xydat[*,0] = trcstr.xcen[gd_cen]
      xydat[*,1] = ycen[gd_cen]
      ;; Set good y region
      mnxy = min(xydat[*,1],max=mxxy)
      mnxy = round(mnxy)
      mxxy = round(mxxy)

      wvdat = double(replicate(1, sz_trc[0]) # trcstr.wav)
  
      ;; Fit
      if not keyword_set( SILENT ) then $
        print, 'esi_echmkaimg: Fitting in 2D: ', systime()
      fit = x_fit2dsurf(xydat, wvdat[gd_cen], FITSTR=fitstr)
      res = fit - wvdat[gd_cen]
      print, 'esi_echmkaimg: RMS = ', fitstr.rms

      sz_subimg = [sz_arc[1], sz_trc[0]]
      ;; Map over good part of image
      numy = mxxy-mnxy+1
      adum = dindgen(sz_subimg[0]) # replicate(1., numy)
      adum2 = replicate(1., sz_subimg[0]) # (dindgen(numy)+mnxy)
  
      ;; Funny format
      xydat = dblarr(n_elements(adum), 2)
      xydat[*,0] = adum
      xydat[*,1] = adum2
      delvarx, adum, adum2
      
      ;; Mapping
      if not keyword_set( SILENT ) then $
        print, 'esi_echkmaimg: Mapping over the image  '
      tot_fit = x_calc2dfit(xydat, FITSTR=fitstr)
      tmp_img = reform(tot_fit,sz_subimg[0], numy)
;      xatv, aimg, /block
      delvarx, tot_fit, xydat

      ;; Place in fin_arc
      arc_img = dblarr(sz_subimg[0],sz_subimg[1]) + 1.d
      arc_img[*,mnxy:mxxy] = temporary(tmp_img)
      ;; KLUDGE for Order 6 and 15
      case qq of 
          0: arc_img[0:1500,mnxy:mxxy] = 0.d
          9: arc_img[2171:4095,mnxy:mxxy] = 0.d
          else:
      endcase
      arc_img = transpose(arc_img)

      ;; Into 'Big Image'
      fin_arc[trcstr.xoff:trcstr.xoff+sz_subimg[1]-1,*] = $
        fin_arc[trcstr.xoff:trcstr.xoff+sz_subimg[1]-1,*] * temporary(arc_img)

  endfor

  if not keyword_set( SILENT ) then $
    print, 'esi_echkmaimg: Mapping back to original frame  '
  fin_arc = transpose(fin_arc)
  fin_arc = transpose(x_invertarc(fin_arc, map, /DBL))
  a = where(fin_arc LT 3000.)
  fin_arc[a] = 0.

  ;; CHK
  if keyword_set( CHK ) then xatv, fin_arc, min=3800., max=11000., /block

; Vacuum wavelengths
  if not keyword_set( SILENT ) then $
    print, 'esi_echmkaimg: Converting to vacuum wavelengths'
  a = where(fin_arc GT 0.)
  tmpaimg = fin_arc[a]
  airtovac, tmpaimg
  fin_arc[a] = temporary(tmpaimg)
  
; Output
  outfil = 'Arcs/ArcECH_'+c_s+'IMG.fits'
  print, 'esi_echmkaimg: Writing file: ', outfil
  mwrfits, fin_arc, outfil, /create, /silent

  ;; Update arc_fil in structure for OBJ+STD
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.slit EQ slit AND $
               (strtrim(esi.type,2) EQ 'OBJ' or $
                strtrim(esi.type,2) EQ 'STD'), nindx)
  esi[indx].arc_fil = outfil

  ;; 
  print, 'esi_echmkaimg: All done!'
      
  return
end
