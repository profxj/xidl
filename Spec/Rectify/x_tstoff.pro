;+ 
; NAME:
; x_tstoff
;    Version 1.1
;
; PURPOSE:
;    Does a quick check on the offset between a flat 
;       and object frame.  Uses FFT correlation
;
; CALLING SEQUENCE:
;   
;   x_tstoff, flat_in, obj_in, offset, CLINE=
;
; INPUTS:
;   flat_in       - Flat image (array or string)
;   obj_in        - Flat image (array or string)
;
; RETURNS:
;   offset - Vertical pixel offset
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ELINE - End line of data (routine takes 1/3, 1/2, 2/3)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_tstoff, flat, obj, offset
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_tstoff_spline, x

common x_tstoff_cmm, splin, spl_x, spl_y

  val = spl_interp(spl_x, spl_y, splin, x, /double)
  return, val[0]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_tstoff, flat_in, obj_in, offset, ELINE=eline

common x_tstoff_cmm

;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_tstoff, flat_in, obj_in, offset, ELINE= [V1.1]'
    return
  endif 

;  Read in the Images

  flat = x_readimg(flat_in, /fscale)
  obj  = x_readimg(obj_in, /fscale)

  if size(obj, /n_elements) NE size(flat, /n_elements) then begin
      print, 'x_tstoff: Images must be the same size!'
      return
  endif

  sz = size(flat, /dimensions)
 

;  Optional Keywords

  if not keyword_set( ELINE ) then eline = sz[0]

  cline = [eline/3, eline/2, 2*eline/3]

  all_off = dblarr(3)

  for qq=0,2 do begin

      ;  Collapse around cline
      flt_coll = djs_median( flat[cline[qq]-15:cline[qq]+15,*], 1)
      obj_coll = djs_median( obj[cline[qq]-15:cline[qq]+15,*], 1)


      ;  FFT's

      fft_flt = fft(flt_coll)
      fft_obj = fft(obj_coll)

      ;  Correlate

      corr = fft( fft_flt * conj(fft_obj), /inverse)

      ; Shift to make things easy

      spl_x = dindgen(20)
      spl_y = -extrac(shift(double(corr),10),0,20)

      ; Spline

      qmn = min(spl_y, imn)
      imn = double(imn)
      spl_y = spl_y + abs(qmn)
      splin = spl_init(spl_x, spl_y, /double)

      ; Find the max

      mn = x_golden('x_tstoff_spline', imn-1., imn, imn+1)

      all_off[qq] = mn-10.d
  endfor


  delvarx, obj, flat

  offset = median(all_off)

  return
end
