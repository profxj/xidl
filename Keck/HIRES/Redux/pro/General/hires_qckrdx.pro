;+ 
; NAME:
; hires_qckrdx   
;     Version 1.1
;
; PURPOSE:
;  Simple routine to reduce an Arc, Flat and Object frame
;
; CALLING SEQUENCE:
;  hires_qckrdx, arc_fil, flat_fil, obj_fil, CHIP=, SETUP=, /NOSKY,
;           /NOOBJ
;
; INPUTS:
;   arc_fil -- Name of arc file
;   flat_fil -- Name of trace flat file
;   obj_fil -- Name of object file 
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
pro hires_qckrdx, afil, ffil, ofil, CHIP=chip, SETUP=setup, _EXTRA=extra, $
                  NOFND=nofnd, NOSKY=nosky, OVERS=overs, HIRES=hires, NOOBJ=noobj, $
                  EXCHK=exchk

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_qckrdx, afil, ffil, ofil, CHIP=, SETUP=, /exchk [v1.1]'
      return
  endif 

  if not keyword_set(CHIP) then chip = [1L,2L,3L]

  ;; List of images
  list = [afil, ffil, ofil]
  
  ;; Create structure
  hires_strct, hires, FILE_LIST=list
  nhires = n_elements(hires)

  ;; Override structure (just in case)
  hires[0:n_elements(afil)*3-1].type = 'ARC'
  hires[n_elements(afil)*3:n_elements(ffil)*3+n_elements(afil)*3-1].type $
    = 'TFLT'
  hires[nhires-3:nhires-1].type = 'OBJ'

  ;; Setup
  hires_setup, hires, XTOLER=0.0015

  bad = where(hires.setup NE 1, nbad)
  if nbad NE 0 then begin
      if keyword_set(OVERS) then hires.setup = 1L else stop
  endif
  if keyword_set(SETUP) then begin
      hires.setup = setup 
      hires_setarcfil, hires, setup
  endif else setup = (hires.setup)[0]

  ;; Flats
  hires_allflat, hires, setup, CHIP, IEXTRAP=[0L,0L], /nogain  
  ;; Take the input gain value

  ;; Arcs
  hires_allarc, hires, setup, CHIP, _EXTRA=extra

  ;; Slitflat  (should avoid doing this)
  if not keyword_set(DOPROF) then begin
      proffil = $
        getenv('XIDL_DIR')+'/Keck/HIRES/Redux/pro/Flat/hires_profile.fits'
      prof0 = xmrdfits(proffil,0,/silent)
      prof1 = xmrdfits(proffil,1,/silent)
      ;; Order structure
      for qq=0L,n_elements(chip)-1 do begin
          ordr_str = hires_getfil('ordr_str', setup, chip=chip[qq],fil_nm=ofil)
          ordr_str.profile0 = prof0
          ordr_str.profile1 = prof1
          mwrfits, ordr_str, ofil, /create
      endfor
  endif

  ;; Object
  if keyword_set(NOOBJ) then return
  obj_id = 1L
  for qq=0L,n_elements(CHIP)-1 do begin
      ;; Process
      hires_proc, hires, SETUP=setup, OBJ=obj_id, CHIP=chip[qq]
      ;; Trace
      if not keyword_set(NOFND) then $
        hires_fntobj, hires, setup, obj_id, chip[qq], CHK=chk
      ;; Skysub
      if not keyword_set(NOSKY) then $
        hires_skysub, hires, setup, obj_id, chip[qq], CHK=chk
      ;; Extract
      hires_extract, hires, setup, obj_id, chip[qq], $
                     CHK=(keyword_set(chk) or keyword_set(EXCHK)), _EXTRA=EXTRA, $
        FIN_TRC=0L
        
  endfor

  return
end

