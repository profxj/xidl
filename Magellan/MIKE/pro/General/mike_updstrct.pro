;+ 
; NAME:
; mike_updstrct   
;     Version 1.1
;
; PURPOSE:
;    Updates the MIKE structure with new images.  
;
; CALLING SEQUENCE:
;  mike_updstrct, struct, LIST=list, EDIT=edit
;
; INPUTS:
;   struct -- Original structure
;
; RETURNS:
;
; OUTPUTS:
;   struct  -  Updated structure list
;
; OPTIONAL KEYWORDS:
;   LIST=  - Image list to add to the mike file  (default: Read all
;            new ones from Raw/)
;  /EDIT   - Launch mike_editstrct
;
; OPTIONAL OUTPUTS:
;   FITS   - Name of FITS file to output MIKE structure into
;
; COMMENTS:
;
; EXAMPLES:
;   mike_updstrct, mike
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_updstrct, mike, LIST=list, FITS=fits, EDIT=edit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_updstrct, mike, LIST=, FITS=, /EDIT (v1.1)'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( FRMCRD ) then frmcrd = 'FRAMENO'
  
;  Grab the files
  if not keyword_set( LIST ) then img = findfile('Raw/*.fits*') $
  else readcol, list, img, FORMAT='A'

  nimg = n_elements(img)
  if nimg EQ 0 then stop

;  Loop on Indv images
  flg_fst = 0
  for q=0,nimg-1 do begin
      head = xheadfits(img[q], /silent)
      name = sxpar(head,'FILENAME')
      frameno = long(strmid(name,1,4))
      a = where(mike.frame EQ frameno, na)
      if na EQ 0 then begin
          if flg_fst EQ 0 then begin
              newlst = [img[q]]
              flg_fst = 1
          endif else newlst = [newlst, img[q]]
      endif
  endfor
          
; Run mike_strct

  if flg_fst EQ 1 then mike_strct, newmike, FILE_LIST=newlst, /NOFILE $
  else begin
      print, 'mike_updstrct: No new files!'
      return
  endelse

; Update

  mike = [mike, newmike]
  srt = sort(mike.frame)
  mike = temporary(mike[srt])

; Edit
  if keyword_set( EDIT ) then mike_editstrct, mike

; Output

  mike_wrstrct, mike, FITS=fits

  return
end
