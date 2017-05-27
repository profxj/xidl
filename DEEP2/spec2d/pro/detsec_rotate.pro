;+
; NAME:
;   detsec_rotate
;
; PURPOSE:
;   Rotate/flip an image according to NOAO mosaic detector standard
;
; CALLING SEQUENCE:
;   detsec_rotate, detsec, image 
;
; INPUTS:
;   detsec     - description of [x1,x2,y1,y2] (NOAO mosaic format)
;   image       - image to rotate/flip
;
; OUTPUTS:
;   same as inputs, overwritten
;
; RESTRICTIONS:
;
; EXAMPLES:
;   detsec_rotate, '[1:1024,8192:4097]', image
;
; COMMENTS:
;   We overwrite both detsec and image
;
; REVISION HISTORY:
;
;       Fri Mar 1 10:01:46 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
pro detsec_rotate, detsec, image

; -------- Parse detsec into one-indexed X1,X2,Y1,Y2
  c = (stregex(detsec, '\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)\]', $
               /subexp, /extract))[1:4] 
  
; rotdir specifies the rotation (or flip) to the rotate function
  rotdir = 0
  d = c
  if c[1] LT c[0] then begin 
     rotdir = rotdir+5          ; flip x
     d[1] = c[0]
     d[0] = c[1]
  endif 
  if c[3] LT c[2] then begin
     rotdir = rotdir+7          ; flip y
     d[3] = c[2]
     d[2] = c[3]
  endif 
  rotdir = rotdir MOD 10

;  print, 'Rotdir: ', rotdir

  if rotdir NE 0 then begin 
     image = rotate(image, rotdir)
     detsec = string('[', d[0], ':', d[1], ',', d[2], ':', d[3], ']')
  endif 

  return
end
