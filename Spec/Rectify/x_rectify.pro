;+ 
; NAME:
; x_rectify   
;    Version 1.1
;
; PURPOSE:
;    Takes the y-distortion map and rectifies an image
;
; CALLING SEQUENCE:
;   
;   rect = x_rectify(img, map, /TRANSP, /SILENT, /DBL, /FIDL)
;
; INPUTS:
;   img       - Input image
;   map       - y-distortion map
;
; RETURNS:
;   rect      - Rectified image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /FIDL  -- Do this using IDL as opposed to a fast C program
;  /TRANSP -- Transpose the image first
;  /DBL   -- Return the rectified image in double precision
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   rect = x_rectify(flat, map)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   29-Jan-2002 Written by JXP
;   13-Aug-2002 Added transpose option (JXP)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_rectify, in_img, in_map, TRANSP=transp, SILENT=silent, DBL=dbl, FIDL=fidl

;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'rect = x_rectify(img, map, /SILENT, /DBL, /FIDL) [V1.1]'
    return,-1
  endif 

  sz_img = size(in_img)
  sz_map = size(in_map)
  if sz_img[0] NE 2 OR sz_map[0] NE 2 then begin
      print, 'x_rectify: need a 2-D image'
      return, -1
  endif
  if sz_img[0] NE sz_map[0] OR sz_img[1] NE sz_map[1] then begin
      print, 'x_rectify: Image and map must have the same size'
      return, -1
  endif

;  Optional Keywords


; Allow for transpose

  if keyword_set( TRANSP ) then begin
      img = transpose(in_img)
      map = transpose(in_map)
  endif else begin
      img = in_img
      map = in_map
  endelse
  

;  Rectify!
  if keyword_set( TRANSP ) then rect = dblarr(sz_img[2],sz_img[1]) $
  else rect = dblarr(sz_img[1],sz_img[2]) 

  ; Ugly For loop
  if not keyword_set(SILENT) then print, 'x_rectify: Rectifying'

  if keyword_set( FIDL ) then begin
      ; ip : Parital pixel (top) in rectified frame
      rows = replicate(1., sz_map[1]) # dindgen(sz_map[2]) 
      diff = rows+map
      ip = long(diff)
      f_ip = 1.d - (diff-ip)
      ; ipp : Partial pixel (bot)
      ipp = ip + 1
      f_ipp = 1.d - f_ip

      for i=0L, sz_img[1]-1 do begin
      ; IDL
          gdip = where(ip[i,*] GE 0 AND ip[i,*] LE (sz_img[2]-1),ngdip)
          if ngdip GT 0 then begin
              tmp = img[i,gdip] * f_ip[i,gdip]
              for j=0L,ngdip-1 do $
                rect[i,ip[i,gdip[j]]] = rect[i,ip[i,gdip[j]]]+tmp[j]
          endif
          gdipp = where(ipp[i,*] GE 0 AND ipp[i,*] LE (sz_img[2]-2),ngdipp)
          if ngdipp GT 0 then begin
              tmp = img[i,gdipp] * f_ipp[i,gdipp]
              for j=0L,ngdipp-1 do $
                rect[i,ipp[i,gdipp[j]]] = rect[i,ipp[i,gdipp[j]]]+tmp[j]
          endif
      endfor
  endif else begin ; C PROGRAM
      ndim = 2L
      soname = filepath('libxmath.' + idlutils_so_ext(), $
;      soname = filepath('libxmath.dylib', $
                        root_dir=getenv('XIDL_DIR'), $
                        subdirectory='/lib')
      retval = call_external(soname, 'rectify', $
                             ndim, size(rect,/dimensions), rect, double(img), $
                             double(map))
  endelse
  if keyword_set( TRANSP ) then rect = transpose(rect)
      

  if keyword_set( DBL ) OR size(img,/type) EQ 5 then return, rect $
  else return, float(rect)

end

