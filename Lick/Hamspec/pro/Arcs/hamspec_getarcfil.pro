;+ 
; NAME:
; hamspec_getarcfil
;     Version 1.1
;
; PURPOSE:
;   Sets the name of the wavelength files for all of the 
;
; CALLING SEQUENCE:
;   hamspec_getarcfil, hamspec, io, aindx, ARC_IMG=, RAW_FIL=
;
; INPUTS:
;   hamspec  -  HIRES structure
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

function hamspec_getarcfil, hamspec, i0, aindx, ARC_IMG=arc_img, RAW_FIL=raw_fil

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'hamspec_getarcfil, hamspec, indx [v1.1]'
      return, -1
  endif 

  ;; Deal with multiple files
  if n_elements(i0) GT 1 then begin
      for ii=0L,n_elements(i0)-1 do begin
          if ii EQ 0 then begin
              arcfil = hamspec_getarcfil(hamspec, i0[ii], aindx, $
                                       ARC_IMG=arc_img, RAW_fil=raw_fil)
          endif else begin
              afil = hamspec_getarcfil(hamspec, i0[ii], aindx, $
                                       ARC_IMG=aimg, RAW_fil=rfil)
              arcfil = [arcfil, afil]
              arc_img = [arc_img, aimg]
              raw_fil = [raw_fil, rfil]
          endelse
      endfor 
  endif

  ;; Set the index
  if not keyword_set(AINDX) then $
    aindx = where(hamspec.setup EQ hamspec[i0].setup AND $
                  hamspec.flg_anly NE 0 AND $
                  hamspec.type EQ 'ARC')

  ;; Find the closest in UT
  diff = abs( (hamspec[i0].date + $
               hamspec[i0].exp/2./3600./24.) - $
              (hamspec[aindx].date + hamspec[aindx].exp/2./3600./24.) )
  mn = min(diff, gdarc)

  ;; Grab the name
  arc_fil = hamspec_getfil('arc_fil', $
                         hamspec[aindx[gdarc]].setup, $
                         /name, CHKFIL=chkf)
  
  arc_img = hamspec_getfil('arc_mkaimg', $
                         hamspec[aindx[gdarc]].setup, $
                         /name, CHKFIL=chkf)
  raw_fil = strtrim(hamspec[aindx[gdarc]].rootpth,2)+ $
    strtrim(hamspec[aindx[gdarc]].img_root,2)

  return, arc_fil
end
