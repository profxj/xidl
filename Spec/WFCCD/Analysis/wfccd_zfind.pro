;+ 
; NAME:
; wfccd_zfind
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   wfccd_zfind, fspec_fil
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_zfind, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_zfind, fspec_fil, obj_nm, ZMIN=zmin, ZMAX=zmax, PLOT=plot, $
                 NPOLY=npoly, WVMNX=wvmnx

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_zfind, fspec_fil, [obj_nm], ZMIN=, ZMAX=, /PLOT, NPOLY= [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( npoly ) then npoly = 5L
  if not keyword_set( zmin ) then zmin = -0.03
  if not keyword_set( zmax ) then zmax = 0.5
  if not keyword_set( WVMNX ) then wvmnx = [3800., 9000.]
  if not keyword_set( SCL_VAR ) then scl_var = 1.0
  if not keyword_set( SKYLIN ) then begin
      nskylin = 6L
      skylin = fltarr(2,10)
      skylin[*,0] = [5569.5, 5593.0]
      skylin[*,1] = [5885.3, 5905.0]
      skylin[*,2] = [6290.3, 6317.0]
      skylin[*,3] = [6351.4, 6380.5]
      skylin[*,4] = [6855.5, 6906.]  ; A band
      skylin[*,5] = [7586.5, 7693.]  ; B band
  endif else begin
      sz_skylin = size(skylin, /dimensions)
      nskylin = sz_skylin[1]
  endelse

;  Read in

  wfccd_wrfspec, wffspec, fspec_fil, /read

;  Spec

  if not keyword_set( OBJ_NM ) then gdsp = where(wffspec.flg_anly NE 0, ngd) $
  else begin
      gdsp = lonarr(1)
      gdsp[0] = x_getobjnm(wffspec,obj_nm)
      ngd = 1L
  endelse
      

;  Loop

  nbad = 0L
  obj_bad = strarr(1000)

  for q=0L,ngd-1 do begin
      print, 'wfccd_zfind: Examining '+strtrim(wffspec[gdsp[q]].slit_id,2)+$
        strtrim(wffspec[gdsp[q]].obj_id,2)+';  Number '+strtrim(q,2)+$
        ' of '+strtrim(ngd-1,2)

      npix = wffspec[gdsp[q]].npix  

      spec_wv = wffspec[gdsp[q]].wave[0:npix-1] 
      spec_fx = wffspec[gdsp[q]].fx[0:npix-1]

      ; Trim
      mnwv = min( abs(spec_wv-wvmnx[0]), imnwv)
      mxwv = min( abs(spec_wv-wvmnx[1]), imxwv)
      gdvar = where(wffspec[gdsp[q]].var[0:npix-1] GT 0.)
      mnvar = min(gdvar, MAX=mxvar)
      
      mnpix = imnwv > mnvar
      mxpix = imxwv < mxvar

      ; Trim off last 3 pix (troubles with rebinning)
      mxpix = mxpix - 3L

      
      ; Final arrays
      npix = mxpix-mnpix+1
      fin_wv = spec_wv[mnpix:mxpix]
      fin_fx = double(spec_fx[mnpix:mxpix])
      fin_var = wffspec[gdsp[q]].var[mnpix:mxpix]
      
      fin_ivar = dblarr(npix)

      posvar = where(fin_var GT 0.)
      fin_ivar[posvar] = 1./(fin_var[posvar]*scl_var)

      ; Kill Sky Lines
      for i=0L,nskylin-1 do begin
          skypix = where(fin_wv GE skylin[0,i] AND $
                         fin_wv LE skylin[1,i], nskypix)
          if nskypix NE 0 then fin_ivar[skypix] = 0.d
      endfor

      
;      if keyword_set( PLOT ) then begin
;          clr = getcolor(/load)
;          !p.color = clr.black
;          !p.background = clr.white
;          posivar = where(fin_ivar GT 0.)
;          fin_sig = fltarr(npix)
;          fin_sig[posivar] = 1./sqrt(fin_ivar[posivar])
;          x_splot, fin_wv, fin_fx, YTWO=fin_sig, /block
;      endif

      ; Create the header
      mkhdr, head, 1, [5000L]
      sxaddpar, head, 'COEFF0', double(alog10(fin_wv[0]))
      sxaddpar, head, 'COEFF1', 100.d/2.9979e5/alog(10.d)

      ; Launch zfind
      zans = x_zfind(fin_fx, fin_ivar, hdr=head, $
                   eigenfile='wfEigenGal-52223.fits', $
                   zmin=zmin, zmax=zmax, doplot=plot, npoly=npoly, /silent)

      ; PLOT
      if keyword_set( PLOT ) then begin
          synth = synthspec(zans, loglam=alog10(fin_wv))
          x_splot, fin_wv, fin_fx, YTWO=synth, /block
      endif

      ; Check on zans
      if keyword_set( PLOT ) then wset, 0
      if zmax-zans.z LT 0.02 OR zans.z_err LT 0. then begin
          print, 'wfccd_zfind: zans ~ zmax!!  Retrying..'
          zans = x_zfind(fin_fx, fin_ivar, hdr=head, $
                       eigenfile='wfEigenGal-52223.fits', $
                       zmin=zmin, zmax=zmax+0.1, DOPLOT=plot, $
                       npoly=(npoly+2L)<6L, /silent)
      endif

      print, 'wfccd_zfind: ', zans.z, zans.z_err, zans.rchi2

      ; Check again
      if zmax-zans.z LT 0.02 OR zans.z_err LT 0. then begin
          print, 'wfccd_zfind: Still trouble in zfind!'
          nbad = nbad + 1
          obj_bad(nbad) = strtrim(wffspec[gdsp[q]].slit_id,2)+$
            strtrim(wffspec[gdsp[q]].obj_id,2)
           ; PLOT
          if keyword_set( PLOT ) then begin
              synth = synthspec(zans, loglam=alog10(fin_wv))
              x_splot, fin_wv, fin_fx, YTWO=synth, /block
          endif
      endif
      

      ; Save
      wffspec[gdsp[q]].zans = zans
  endfor

  ; bad fits
  if nbad GT 0L then begin
      print, 'wfccd_zfind: Bad fits!  Check output!'
      print, obj_bad[1:nbad]
  endif

  ; Write back
  wfccd_wrfspec, wffspec, fspec_fil

  print, 'wfccd_zfind: All done!'
  return
end
