;+ 
; NAME:
; hires_editord
;     Version 1.0
;
; PURPOSE:
;    Modifies order structure files (it is wise to back them up before
;     using this utility. E.g. OStr_B_01.fits)
;
; CALLING SEQUENCE:
;  hires_editord, hires, setup, chip
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup
;   chip    -  Chip to process
;
; RETURNS:
;
; OUTPUTS:
;  Modifies the selected (via inputs) ordr_str file.
;
; OPTIONAL KEYWORDS:
;  /LISTORDS -- List order IDs and physical IDs
;  DELORDS -- Delete orders listed in this array (by index by default)
;  /BYINDEX -- DELORDS lists order indexes
;  /BYPHYSID -- DELORDS lists physical order ids (may not be available 
;               early in the reduction)
;  PHYSIDS -- Array of physical order IDs to replace those in the resulting
;             order structure
; 
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_editord, hires, 1, 1, /LISTORDS
;   hires_editord, hires, 1, 1, /LISTORDS, DELORDS=[0L,1L], /BYINDEX, PHYSIDS=[85L,86L,87L,88L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Nov-2011 Written by ALM
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_editord, hires, setup, chip, LISTORDS=listords, $
  DELORDS=delords, BYINDEX=byindex, BYPHYSID=byphysid, PHYSIDS=physids

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_editord, hires, setup, chip, /LISTORDS, DELORDS='
      print, '        /BYINDEX, /BYPHYSID, PHYSIDS= [v1.0]'
      return
  endif
  
  if not keyword_set(BYINDEX) or not keyword_set(BYPHYSID) then byindex = 1
  
  ordr_fil = hires_getfil('ordr_str', setup, CHIP=chip, CHKFIL=chkf, /name)
  if CHKF EQ 0 then begin
      print, 'hires_editord: Order structure does not exist!'
      return
  endif
  
  ord_str = xmrdfits(ordr_fil, 1, /silent)
  iids = lindgen(n_elements(ord_str))
  pids = ord_str.order
  
  if keyword_set(LISTORDS) then begin
    print, 'RED <----- index -----> BLUE'
    print, iids
    print, 'RED <-- physical id --> BLUE'
    print, pids
    print, 'hires_editord: Nothing deleted.'
    print, 'hires_editord: Done!!!'
    return
  endif
  
  if keyword_set(DELORDS) then begin
    keep = iids
    if keyword_set(BYINDEX) then match2, iids, delords, delids, sub_delords
    if keyword_set(BYPHYSID) then match2, pids, delords, delids, sub_delords
    keep = where(delids eq -1)
    print, 'hires_editord: keeping orders: '
    print, 'RED <----- index -----> BLUE'
    print, keep
    print, 'RED <-- physical id --> BLUE'
    print, pids[keep]
    
    ord_str = ord_str[keep]
  endif
  
  if keyword_set(PHYSIDS) then begin
    if n_elements(PHYSIDS) eq n_elements(ord_str) then begin
      ord_str.order = physids
    endif else begin
      print, 'hires_editord: Number of physical order ids does not match' + $
        ' that in order structure!'
      return
    endelse
  endif
  
  if keyword_set(PHYSIDS) or keyword_set(DELORDS) then begin
    print, 'hires_editord: overwriting order structure file ', ordr_fil
    mwrfits, ord_str, ordr_fil, /create
  endif
  
  print, 'hires_editord: Done!!!'
  return
end
