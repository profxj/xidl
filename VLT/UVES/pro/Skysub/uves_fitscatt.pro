;+ 
; NAME:
; uves_fitscatt
;     Version 1.2
;
; PURPOSE:
;    Fit scattered light image to the light in the order gaps
;
; CALLING SEQUENCE:
;   
;  uves_fitscatt, img, ivar, ordr_str, side
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   A fit to the scattered light 
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_fitgap
;  x_ordermask
;
; REVISION HISTORY:
;   02-Sep-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_fitscatt, img, ivar, ordr_str, side, SCATTTRIM=scatttrim, $
  NOMEDSCATT=nomedscatt, CHK=chk, NXBKPT=nxbkpt

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'uves_fitscatt, img, ivar, ordr_str, side [v1.1]'
      return, -1
  endif 

  if not keyword_set(SCATTTRIM) then scatttrim = 1.
  
  ;; Inter-order light fit (small for HIRES)
  case side of 
      1: begin
          if not keyword_set(NXBKPT) then nxbkpt = 4
          nybkpt = 4
          if not keyword_set(NOMEDSCATT) then MEDSCATT=1L
      end
      2: begin
          if not keyword_set(NXBKPT) then nxbkpt = 4
          nybkpt = 4
          MEDSCATT=0L
      end
  endcase

  sz_img = size(img, /dimensions)

  ;; Mask image (identifies the interorder gaps)
  maskimage = x_ordermask(sz_img[0], sz_img[1], ordr_str, trim=scatttrim)

  print, 'uves_fitscatt: Fitting inter-order light with x_fitgap', $
    nxbkpt, ' x ', nybkpt, ' bkpts..', format='(a,i3,a,i2,a)'
  scatt_img = x_fitgap(img, ivar, maskimage, $
                             nxbkpt=nxbkpt, nybkpt=nybkpt)

  ;; Median?
  if keyword_set(MEDSCATT) then begin
      meds = median(scatt_img)
      print, 'uves_fitscatt:  Using median of scattered light image', meds
      scatt_img[*] = meds
  endif

  ;; Error chk
  if min(scatt_img) LT -50. then begin
      print, 'uves_fitscatt: I suspect you have a problem here!, '
      print, ' If you continue, scattered light will have 0 set as minimum'
      xatv, scatt_img, /block
      stop  ;; Problem likely
      scatt_img = scatt_img > 0
  endif

  if keyword_set(CHK) then xatv, scatt_img, /bloc


  return, scatt_img
end

