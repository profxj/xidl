;+ 
; NAME:
; uves_getarcfil
;     Version 1.1
;
; PURPOSE:
;   Sets the name of the wavelength files for all of the 
;
; CALLING SEQUENCE:
;   uves_getarcfil, uves, io, aindx, ARC_IMG=, RAW_FIL=
;
; INPUTS:
;   uves  -  HIRES structure
;   i0     -  Obj index
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

function uves_getarcfil, uves, i0, aindx, ARC_IMG=arc_img, RAW_FIL=raw_fil

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'uves_getarcfil, uves, indx [v1.1]'
      return, -1
  endif 

  ;; Deal with multiple files
  if n_elements(i0) GT 1 then begin
      for ii=0L,n_elements(i0)-1 do begin
          if ii EQ 0 then begin
              arcfil = uves_getarcfil(uves, i0[ii], aindx, $
                                       ARC_IMG=arc_img, RAW_fil=raw_fil)
          endif else begin
              afil = uves_getarcfil(uves, i0[ii], aindx, $
                                       ARC_IMG=aimg, RAW_fil=rfil)
              arcfil = [arcfil, afil]
              arc_img = [arc_img, aimg]
              raw_fil = [raw_fil, rfil]
          endelse
      endfor 
  endif

  ;; Set the index
  if not keyword_set(AINDX) then $
    aindx = where(uves.setup EQ uves[i0].setup AND $
                  uves.flg_anly NE 0 AND $
                  uves.side EQ uves[i0].side AND $
                  uves.type EQ 'ARC')


  ;; WCEN
  wcen = round(uves[i0].xdangl)

  ;; Find the closest in UT
  diff = abs( (uves[i0].date + $
               uves[i0].exp/2./3600./24.) - $
              (uves[aindx].date + uves[aindx].exp/2./3600./24.) )
  mn = min(diff, gdarc)

  ;; Grab the name
  arc_fil = uves_getfil('arc_fil', uves[i0].setup, $
                         FRAME=uves[aindx[gdarc]].frame, $
                         WCEN=wcen, /name, CHKFIL=chkf)
  
  arc_img = uves_getfil('arc_mkaimg', uves[i0].setup, $
                         FRAME=uves[aindx[gdarc]].frame, $
                         WCEN=wcen, /name, CHKFIL=chkf)
  raw_fil = strtrim(uves[aindx[gdarc]].rootpth,2)+ $
    strtrim(uves[aindx[gdarc]].img_root,2)

  return, arc_fil
end
