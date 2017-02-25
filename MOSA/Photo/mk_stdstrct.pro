;+ 
; NAME:
; mk_stdstrct   
;   Version 1.1
;
; PURPOSE:
;    Create a direct image structure for a list of standard stars
;
; CALLING SEQUENCE:
;   mk_stdstrct, strct, outfil
;
; INPUTS:
;
; RETURNS:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Jan-2004 Written by JXP 
;-
;------------------------------------------------------------------------------
pro mk_stdstrct, strct, outfil

  if not keyword_set( OUTFIL ) then outfil = 'std_strct.fits'

  ;; Create structure
  xdimg_strct, strct, 'MOSA', 'KPNO-4m', LIST='Lists/std.list', $
    /NOFILE, /NOSTATS, IMG=img
  nimg = n_elements(img)

  ;; Update filter
  for q=0L,nimg-1 do begin
      case strtrim(strct[q].filter,2) of
          '1': strct[q].filter = 'U'
          '2': strct[q].filter = 'B'
          '3': strct[q].filter = 'V'
          '4': strct[q].filter = 'R'
          '5': strct[q].filter = 'I'
          else: stop
      endcase
  endfor

  ;; Type
  strct.type = 'STD'

  ;; Final image
  strct.img_final = img
  strct.flg_final = 1

  ;; Output
  mwrfits, strct, outfil, /create

  return
end
