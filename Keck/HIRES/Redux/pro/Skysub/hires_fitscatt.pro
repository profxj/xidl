;+ 
; NAME:
; hires_fitscatt
;     Version 1.2
;
; PURPOSE:
;    Fit scattered light image to the light in the order gaps.  The
;    main program is  x_fitgap
;
; CALLING SEQUENCE:
;   
;  hires_fitscatt, img, ivar, ordr_str, chip
;
; INPUTS:
;  img -- Image for subtracting
;  ivar -- Inverse variance of the image
;  ordr_str -- Order structure defining the echelle footprint
;  chip    -  Blue (1), Green (2) OR Red (3) chip
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SCATTRIM= -- Number of pixels of the edge of the slit to trim
;               before fitting scattered light
;  /NOMEDSCATT -- Override any attempt to characterize the scattered
;                 light by only the median value
;  NXBKPT= -- Number of coefficients to use in x direction for
;                scattered light subtraction.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_fitgap
;
; REVISION HISTORY:
;   02-Sep-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_fitscatt, img, ivar, ordr_str, chip, SCATTTRIM=scatttrim, $
                         NOMEDSCATT=nomedscatt, NXBKPT=nxbkpt

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_fitscatt, img, ivar, ordr_str, chip, NXBPKT= [v1.1]'
      return, -1
  endif 

  if not keyword_set(SCATTTRIM) then scatttrim = 1.
  
  ;; Inter-order light fit (small for HIRES)
  case chip of 
      1: begin
          if not keyword_set(NXBKPT) then nxbkpt = 3
          nybkpt = 3
          if not keyword_set(NOMEDSCATT) then MEDSCATT=1L
      end
      2: begin
          if not keyword_set(NXBKPT) then nxbkpt = 5
          nybkpt = 5
          MEDSCATT=0L
      end
      3: begin
          if not keyword_set(NXBKPT) then nxbkpt = 5
          nybkpt = 5
          MEDSCATT=0L
      end
  endcase

  sz_img = size(img, /dimensions)

  ;; Mask image (identifies the interorder gaps)
  maskimage = x_ordermask(sz_img[0], sz_img[1], ordr_str, trim=scatttrim)

  print, 'hires_fitscatt: Fitting inter-order light with x_fitgap', $
    nxbkpt, ' x ', nybkpt, ' bkpts..', format='(a,i3,a,i2,a)'
  scatt_img = x_fitgap(img, ivar, maskimage, $
                             nxbkpt=nxbkpt, nybkpt=nybkpt)

  ;; Median?
  if keyword_set(MEDSCATT) then begin
      meds = median(scatt_img)
      print, 'hires_fitscatt:  Using median of scattered light image', meds
      scatt_img[*] = meds
  endif

  ;; Error chk
  if min(scatt_img) LT -10. then begin
      print, 'hires_fitscatt: I suspect you have a problem here!, '
      print, ' If you continue, scattered light will have 0 set as minimum'
      xatv, scatt_img, /block
      stop  ;; Problem likely
      scatt_img = scatt_img > 0
  endif


  return, scatt_img
end

