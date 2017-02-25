;+ 
; NAME:
; dla_updkin
;  V1.0
;
; PURPOSE:
;    Given a list of DLA .dat files, measure the kinematic
;    characteristics, modify the .dat files and write everything out. 
;
; CALLING SEQUENCE:
;   dla_updkin, dla, /CHK
;
; INPUTS:
;  list -- List of DLA .dat files.
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CORR_VAL= -- Value to offset ESI measurements by
;  /CHK  -- Plot the transitions and kinematic characteristics
;  /CIV  -- Run CIV kinematics too
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_updabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP 
;   14-Feb-2012 Added kbin keyword (MN)
;-
;------------------------------------------------------------------------------
; dla_updkin, 'Lists/tst_dla.lst', SPEC_ROOT=getenv('DLA')+'/Kin/', ROOT=getenv('DLA')
pro dla_updkin, list, ROOT=root, CHK=chk, CORR_VAL=corr_val, $
                SPEC_ROOT=spec_root, kbin=kbin, CIV=civ

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_updkin, list, /CHK, ROOT=, CORR_VAL=  [v1.0]'
    return
  endif 

  ;; Fix Zn check
  if keyword_set(FIXZN) then begin
      print, 'dla_updabd: Warning you are about to fiddle with Zn!'
      stop
   endif



  close, /all
  ;; Parse
  parse_dlalst, dla, list, /noelm, ROOT=root
  ndla = n_elements(dla)

  ;; Loop
  for nn=0L,ndla-1 do begin
      ;; Low-ion kinematics
      print, dla[nn].qso, ' ', dla[nn].zabs
      ism_lowkin, dla, nn, CHK=chk, SPEC_ROOT=spec_root, KBIN=kbin
      if keyword_set(CIV) then $
         ism_civkin, dla, nn, CHK=chk, SPEC_ROOT=spec_root, KBIN=kbin

      ;; Correction?
      if keyword_set(CORR_VAL) and dla[nn].flglw GT 0 then begin
          ;; Check for ESI
          pos = strpos(dla[nn].lwfil, 'ESI')+strpos(dla[nn].lwfil, 'esi')+1
          if pos GT 0 then begin
              print, 'Correcting ', dla[nn].qso, ' for ESI by ', corr_val, 'km/s'
              dla[nn].lwfvel = dla[nn].lwfvel - CORR_VAL
          endif
      endif

      ;; Correction?
      if keyword_set(CIV) then begin
         if keyword_set(CORR_VAL) and dla[nn].flgciv GT 0 then begin
            ;; Check for ESI
            pos = strpos(dla[nn].civfil, 'ESI')+strpos(dla[nn].civfil, 'esi')+1
            if pos GT 0 then begin
               print, 'Correcting ', dla[nn].qso, ' for ESI by ', corr_val, 'km/s'
               dla[nn].civfvel = dla[nn].civfvel - CORR_VAL
            endif
         endif
      endif
  endfor


  ;; Write
  dla_writestr, dla
  print, 'dla_updkin: All done'

  return
end
