;+ 
; NAME:
; hires_getarcfil
;     Version 1.1
;
; PURPOSE:
;   Sets the name of the wavelength files for all of the 
;
; CALLING SEQUENCE:
;   hires_getarcfil, hires, io, aindx, ARC_IMG=, RAW_FIL=
;
; INPUTS:
;   hires  -  HIRES structure
;   i0     -  Setup numbers
;   aindx  -  
;
; RETURNS:
;
; OUTPUTS:
;  A series of files related to Arc calibration
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
;   15-Aug-2003 Written by JXP
;
;  Usage:
;
;-
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_getarcfil, hires, i0, aindx, ARC_IMG=arc_img, RAW_FIL=raw_fil

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'hires_getarcfil, hires, indx [v1.1]'
      return, -1
  endif 

  ;; Deal with multiple files
  if n_elements(i0) GT 1 then begin
      for ii=0L,n_elements(i0)-1 do begin
          if ii EQ 0 then begin
              arcfil = hires_getarcfil(hires, i0[ii], aindx, $
                                       ARC_IMG=arc_img, RAW_fil=raw_fil)
          endif else begin
              afil = hires_getarcfil(hires, i0[ii], aindx, $
                                       ARC_IMG=aimg, RAW_fil=rfil)
              arcfil = [arcfil, afil]
              arc_img = [arc_img, aimg]
              raw_fil = [raw_fil, rfil]
          endelse
      endfor 
  endif

  ;; Set the index
  if not keyword_set(AINDX) then $
    aindx = where(hires.setup EQ hires[i0].setup AND $
                  hires.flg_anly NE 0 AND $
                  hires.chip EQ hires[i0].chip AND $
                  hires.type EQ 'ARC')

  ;; Find the closest in UT
  diff = abs( (hires[i0].date + $
               hires[i0].exp/2./3600./24.) - $
              (hires[aindx].date + hires[aindx].exp/2./3600./24.) )
  mn = min(diff, gdarc)

  ;; Grab the name
  tchip = hires[i0].chip
  arc_fil = hires_getfil('arc_fil', $
                         FRAME=hires[aindx[gdarc]].frame, $
                         CHIP=tchip, /name, CHKFIL=chkf)
  
  arc_img = hires_getfil('arc_mkaimg', $
                         FRAME=hires[aindx[gdarc]].frame, $
                         CHIP=tchip, /name, CHKFIL=chkf)
  raw_fil = strtrim(hires[aindx[gdarc]].rootpth,2)+ $
    strtrim(hires[aindx[gdarc]].img_root,2)

  return, arc_fil
end
