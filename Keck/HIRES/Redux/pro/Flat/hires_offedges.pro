;+ 
; hires_offedges
;     Version 1.1
;
; PURPOSE:
;   Offset edges of the slit for specific orders
;
; CALLING SEQUENCE:
;   
;  hires_offedges, hires, setup, chip, off_fil
;
; INPUTS:
;   hires    -  HIRES structure
;   setup    -  Setup identifier 
;   chip     -  Chip
;   off_fil  -  File with offsets.  Format is Order#, left offset,
;              right offset
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
;   hires_offedges, hires, 1, 1, 'Flats/offset_01_B.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_getfil
;
; REVISION HISTORY:
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_offedges, hires, setup, chip, off_fil
;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_chktrcflat, hires, setup, chip, off_fil [v1.0]'
      return
  endif 

;  Optional Keywords

  ;; Read offset file
  readcol, off_fil, order, left_off, right_off, format='L,F,F'
  noff = n_elements(order)
  
  ;; Order structure
  ordr_str = hires_getfil('ordr_str', setup, CHIP=chip, FIL_NM=filnm)

  ;; Offset
  for qq=0L,noff-1 do begin
     idx = where(ordr_str.order EQ order[qq], nid)
     if nid EQ 0 then begin
        print, 'hires_offedges: No order ', order[qq]
        continue
     endif

     ordr_str[idx].lhedg = ordr_str[idx].lhedg + left_off[qq]
     ordr_str[idx].lhc = ordr_str[idx].lhc + left_off[qq]

     ;if order[qq] EQ 97 then stop
     ordr_str[idx].rhedg = ordr_str[idx].rhedg + right_off[qq]
     ordr_str[idx].rhc = ordr_str[idx].rhc + right_off[qq]
  endfor

  ;; Write
  print, 'hires_offedges: You are about to offset the orders.  This will be hard to undo. '
  print, 'hires_offedges: Continue with caution...'
  stop
  mwrfits, ordr_str, filnm, /create

  return
end
