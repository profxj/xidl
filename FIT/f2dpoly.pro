;+ 
; NAME:
; f2dpoly
;   Version 1.1
;
; PURPOSE:
;    Creates basis functions of a 2dpoly for feeding into SVDFIT for
;    2D surface fitting. Requires an initial call with m=-1 to set 
;    up the common block.  This code even makes my head spin!
;
; CALLING SEQUENCE:
;   
;   fpoly = f2dpoly(s, m, XVAL=, YVAL=, FLG=)
;
; INPUTS:
;   s          - scalar or vector identifying the index number
;   m          - Total order (nx*ny) of the polynomial (-1 to
;                initialize; -2 to deconstruct)
;
; RETURNS:
;   fpoly      - Basis functions
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XVAL=  -- Dummy array used to initialize the common block
;  YVAL=  -- Dummy array used to initialize the common block
;  FLG=  -- Number of coefficients in the X direction.   
;           If FLG=0, then it is assumed that nx=ny=sqrt(m)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fpoly = f2dpoly(s, m)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function f2dpoly, s, m, XVAL=xval, YVAL=yval, FLG=flg

common f2dpoly_common, f2dpoly_x, f2dpoly_y, f2dpoly_flg

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'basis = f2dpoly(s, m) [v1.1]'
    return, -1
  endif 

;  Initialize

  if m EQ -1 then begin
      f2dpoly_x = xval
      f2dpoly_y = yval
      f2dpoly_flg = flg
      return, 0
  endif

;  Deconstruct

  if m EQ -2 then begin
      delvarx, f2dpoly_x, f2dpoly_y, f2dpoly_flg
      return, -1
  endif

;  Begin

  ns = n_elements(s)
  ls = long(s)

  case f2dpoly_flg of
      0 : begin
          sm = sqrt(m)
          if abs(sm - long(sm)) GT 0.001 then begin
              print, 'f2dpoly: This flag requires m to be a perfect square'
              return, -1
          endif
          nx = long(sm)
          ny = nx
      end
      else : begin
          nx = f2dpoly_flg
          ny = m/f2dpoly_flg
      end
  endcase

; Loop

  if ns GT 1 then begin
      f2dpoly = make_array(ns, m, type=size(f2dpoly_x,/type)) + 1.
      nyrep = replicate(1., ny)
      nxrep = replicate(1., nx)
      ; x basis first
      for i = 1,nx-1 do $
        f2dpoly[0:ns-1,i*ny:i*ny+ny-1] = f2dpoly[*,i*ny-1]*f2dpoly_x[ls] # nyrep
      ; y basis next
      dumi = lindgen(nx)*ny     ; x entries
      for j = 1,ny-1 do $
        f2dpoly[0:ns-1,(dumi+j)] = f2dpoly[*,(dumi+j-1)] * (f2dpoly_y[ls] # nxrep)
  endif else begin
      ls = ls[0]
      f2dpoly = make_array(m, type=size(s,/type)) + 1.
                                ; x basis first
      for i = 1,nx-1 do f2dpoly[i*ny:i*ny+ny-1] = f2dpoly[i*ny-1]*f2dpoly_x[ls]
                                ; y basis next
      dumi = lindgen(nx)*ny     ; x entries
      for j = 1,ny-1 do f2dpoly[(dumi+j)] = f2dpoly[(dumi+j-1)] * f2dpoly_y[ls]
  endelse

  return, f2dpoly
end

