;+ 
; NAME:
; hires_pipe   
;     Version 1.1
;
; PURPOSE:
;  Routine to run the full HIRES redux, primarily in concert with 
;   the UC/Keck Data Reduction website
;
; CALLING SEQUENCE:
;  hires_pipe, arc_fil, flat_fil, obj_fil, CHIP=, SETUP=, /NOSKY,
;           /NOOBJ
;
; INPUTS:
;   obj_fil -- Name of object file 
;   arc_fil -- Name of arc file
;   flat_fil -- Name of trace flat file
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CHIP= -- Restrict analyis to a single/pair of chips [Default: All]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_qckrdx, 'hires0023.fits',
;   ['hires0030.fits','hires0031.fits'], 'hires0082.fits', CHIP=1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-May-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_pipe, obj_fil, arcs, flats, CHIP=chip, SETUP=setup, _EXTRA=extra, $
                  NOFND=nofnd, NOSKY=nosky, OVERS=overs, $
                  HIRES=hires, NOOBJ=noobj, EXCHK=exchk

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_qckrdx, afil, ffil, ofil, CHIP=, SETUP=, /exchk [v1.1]'
      return
  endif 

  if not keyword_set(CHIP) then chip = [1L,2L,3L]

  if size(flats, /type) NE 7 then begin
     print, 'hires_pipe: You need to submit some flats!'
     return
  endif
  if size(arcs, /type) NE 7 then begin
     print, 'hires_pipe: You need to submit some arc frames!'
     return
  endif

  ;; List of images
  list = [arcs, flats, obj_fil]
  
  ;; Create structure
  hires_strct, hires, FILE_LIST=list
  nhires = n_elements(hires)

  ;; Over-ride structure (just in case)
  hires[0:n_elements(arcs)*3-1].type = 'ARC'
  hires[n_elements(arcs)*3:n_elements(flats)*3+n_elements(arcs)*3-1].type $
    = 'TFLT'
  hires[nhires-3:nhires-1].type = 'OBJ'

  ;; Setup
  hires_setup, hires, XTOLER=0.0015

  bad = where(hires.setup NE 1, nbad)
  if nbad NE 0 then begin
      if keyword_set(OVERS) then hires.setup = 1L $
      else begin
          print, 'hires_pipe:  More than one setup from these files.'
          print, 'hires_pipe:  Reduce input files to all have the same setup.'
          print, 'hires_pipe:  Check binning, ECHANGL, XDANGL, DECKER'
          return
      endelse
  endif
  if keyword_set(SETUP) then begin
      hires.setup = setup 
      hires_setarcfil, hires, setup
  endif else setup = (hires.setup)[0]

  ;; Flats
  hires_allflat, hires, setup, CHIP

  ;; Arcs
  hires_setarcfil, hires, setup
  hires_allarc, hires, setup, CHIP, _EXTRA=extra

  ;; Slitflat  (should avoid doing this)
  hires_slitflat, hires, setup 
;  if not keyword_set(DOPROF) then begin
;      proffil = $
;        getenv('XIDL_DIR')+'/Keck/HIRES/Redux/pro/Flat/hires_profile.fits'
;      prof0 = xmrdfits(proffil,0,/silent)
;      prof1 = xmrdfits(proffil,1,/silent)
;      ;; Order structure
;      for qq=0L,n_elements(chip)-1 do begin
;          ordr_str = hires_getfil('ordr_str', setup, chip=chip[qq],fil_nm=ofil)
;          ordr_str.profile0 = prof0
;          ordr_str.profile1 = prof1
;          mwrfits, ordr_str, ofil, /create
;      endfor
;  endif

  ;; Object
  if keyword_set(NOOBJ) then return

  objs = where(hires.type EQ 'OBJ' and hires.setup EQ setup AND $
               hires.flg_anly EQ 1, nobjs)
  all_objid = hires[objs[uniq(hires[objs].obj_id, $
                              sort(hires[objs].obj_id))]].obj_id
  nobjid = n_elements(all_objid)

  for kk=0L,nobjid-1 do begin
      obj_id = all_objid[kk]
      for qq=0L,n_elements(CHIP)-1 do begin
          ;; Process
          hires_proc, hires, SETUP=setup, OBJ=obj_id, CHIP=chip[qq]
          ;; CR
          hires_objcr, hires, setup, obj_id, chip[qq]
          
          ;; Trace
          if not keyword_set(NOFND) then $
            hires_fntobj, hires, setup, obj_id, chip[qq], CHK=chk
          
          ;; Skysub
          if not keyword_set(NOSKY) then $
            hires_skysub, hires, setup, obj_id, chip[qq], CHK=chk
          
          ;; Extract
          hires_extract, hires, setup, obj_id, chip[qq], $
                         CHK=(keyword_set(chk) or $
                              keyword_set(EXCHK)), _EXTRA=EXTRA, $
                         FIN_TRC=0L
      endfor
  endfor

  return
end

