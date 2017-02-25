;+ 
; NAME:
; x_origslit   
;    Version 1.0
;
; PURPOSE:
;    Given the slitstr and the map, find slit positions in the
;    original image.  The values are filled into the slit structure.
;
; CALLING SEQUENCE:
;   x_origslit, slitstr, map, /INVERSE
;
; INPUTS:
;   slitstr     - Slit structure
;   map         - y-distortion map (fits is ok)
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;  /INVERSE  -- The map is the inverse!
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_origslit, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_origslit, slitstr, map, INVERSE=inverse


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_origslit, slitstr, map, /INVERSE [v1.1]'
    return
  endif 


;  Optional Keywords

; Allow map to be fits file

  dat = x_readimg(map, /fscale)
  sz_map = size(dat, /dimensions)
  npix = sz_map[0]

; Slits
  nslit = n_elements(slitstr)

; Inverse Map
  if keyword_set( INVERSE ) then begin
      for ii=0L,nslit-1 do begin
          ; Top
          iy = long(slitstr[ii].yedg_flt[1])
          frac = slitstr[ii].yedg_flt[1] - iy
          ; Interpolate
          slitstr[ii].yedg_orig[0:npix-1,1] = $
            frac*map[*,iy+1] + (1-frac)*map[*,iy] + slitstr[ii].yedg_flt[1]
          ; Bottom
          iy = long(slitstr[ii].yedg_flt[0])
          frac = slitstr[ii].yedg_flt[0] - iy
          ; Interpolate
          slitstr[ii].yedg_orig[0:npix-1,0] = $
            frac*map[*,iy+1] + (1-frac)*map[*,iy] + slitstr[ii].yedg_flt[0]
      endfor
  endif else begin
      yy = findgen(sz_map[1])
      ; Loop on pixels
      for ii=0L,npix-1 do begin
          gy = yy + map[ii,*]   ; Map
          splin = spl_init(yy, gy, /double)
          ; Loop on slits
          for jj=0L,nslit-1 do begin
              for kk=0,1 do begin
                  mn = min(abs(slitstr[jj].yedg_flt[kk]-gy),imn)
                  slitstr[jj].yedg_orig[ii,kk] = $
                    x_fndspln(yy, gy, slitstr[jj].yedg_flt[kk], splin, $
                              IPX=((imn-5)>0), TOLER=1e-6)
              endfor
          endfor
      endfor
  endelse

  ; All done
  delvarx, dat

end

