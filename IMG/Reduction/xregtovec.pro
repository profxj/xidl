;+ 
; NAME:
; xregtovec
;   Version 1.1
;
; PURPOSE:
;    Converts a string region into a vector for an image
;
; CALLING SEQUENCE:
;   
;   vect = xregtovec(sreg, size)
;
; INPUTS:
;   string     - String defining the region (e.g. '[1:10,*]')
;   size       - Size of the image
;
; RETURNS:
;   vect       - Integer vector defining the region in sreg
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   area = xregtovec('[1:10,*]',size(img))
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   21-June-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function xregtovec, sreg, size

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'xregtovec, sreg, [size] [v1.1]'
    return, -1
  endif 

  if strpos(sreg, '*') NE -1 and not keyword_set(SIZE) then return, -1

;  Optional Keywords

;  if keyword_set( OVSEC ) then    flgovsec    = 1 else flgovsec = 0

  area = lindgen(4)

;  Parsing x 

  if strmid(sreg,0,1) NE '[' then $
    stop, 'Region must begin with [ in xregtovec', sreg

  if strmid(sreg,1,1) EQ '*' then begin
      if strmid(sreg,2,1) NE ',' then $
        stop, 'Badly formed region in xregtovec', sreg
      cchr = 3
      area[0] = 0
      area[1] = size[1]-1
  endif else begin
      nchr = strpos(sreg,':')
      area[0] = fix(strmid(sreg,1,nchr-1))
      cchr = strpos(sreg,',')
      area[1] = fix(strmid(sreg,nchr+1,cchr-nchr-1))
  endelse
      
;  Parsing y

  yreg = strmid(sreg, cchr+1)
      
  if strmid(yreg,0,1) EQ '*' then begin
      if strmid(yreg,1,1) NE ']' then $
        stop, 'Badly formed region in xregtovec', yreg
      area[2] = 0
      area[3] = size[2]-1
  endif else begin
      nchr = strpos(yreg,':')
      area[2] = fix(strmid(yreg,0,nchr))
      cchr = strpos(yreg,']')
      area[3] = fix(strmid(yreg,nchr+1,cchr-nchr-1))
  endelse
      
  return, area
end

