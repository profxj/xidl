;+ 
; NAME:
; hires_updstrct   
;     Version 1.1
;
; PURPOSE:
;    Updates the MIKE structure with new images.  
;
; CALLING SEQUENCE:
;  hires_updstrct, struct, LIST=list, EDIT=edit
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
;   LIST=  - Image list to add to the hires file  (default: Read all
;            new ones from Raw/)
;  /EDIT   - Launch hires_editstrct
;
; OPTIONAL OUTPUTS:
;   FITS   - Name of FITS file to output MIKE structure into
;
; COMMENTS:
;
; EXAMPLES:
;   hires_updstrct, hires
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_updstrct, hires, LIST=list, FITS=fits, EDIT=edit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_updstrct, hires, LIST=, FITS=, /EDIT (v1.1)'
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
;      name = sxpar(head,'FILENAME')
      frameno = sxpar(head,FRMCRD)
      a = where(hires.frame EQ frameno, na)
      if na EQ 0 then begin
          if flg_fst EQ 0 then begin
              newlst = [img[q]]
              flg_fst = 1
          endif else newlst = [newlst, img[q]]
      endif
  endfor
          
; Run hires_strct

  if flg_fst EQ 1 then hires_strct, newhires, FILE_LIST=newlst, /NOFILE $
  else begin
      print, 'hires_updstrct: No new files!'
      return
  endelse

; Update

  hires = [hires, newhires]
  srt = sort(hires.frame)
  hires = temporary(hires[srt])

; Edit
  if keyword_set( EDIT ) then hires_editstrct, hires

; Output

  hires_wrstrct, hires, FITS=fits

  return
end
