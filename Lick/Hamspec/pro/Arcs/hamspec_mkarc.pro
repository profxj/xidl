;+ 
; NAME:
; hamspec_echmkarc   
;     Version 1.1
;
; PURPOSE:
;    Process and combine arc files  
;
; CALLING SEQUENCE:
;   
;  hamspec_echmkarc, hamspec, slit, /CLOBBER
;
; INPUTS:
;   hamspec     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_ECH##.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_echmkarc, hamspec, 0.5
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;   29-Apr-2007 Removed option to flat field arc (JFH)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_mkarc, hamspec, setup, CLOBBER=clobber, BIASFIL=biasfil $
                  , SEDG_FIL = SEDG_FIL, FLATFIL=flatfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_mkarc, hamspec, slit, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
;; Binning

  arcs = where(hamspec.setup EQ setup AND hamspec.flg_anly NE 0 AND $
              strtrim(hamspec.type,2) EQ 'ARC', narc)
  if narc EQ 0 then begin
      print, 'hamspec_echmkarc: No Arcs found! Returning' 
      return
  endif

  ;; Check for existing file
  outfil = hamspec_getfil('arc_fil', setup, /name)
  a = findfile(outfil+'*', count=na)
  if na NE 0 AND not keyword_set( CLOBBER ) then begin
     print, 'hamspec_mkarc: Arc ', outfil, ' exists.  Returning'
     return
  endif

  ;; Sum em up
  
  for qq=0L,narc-1 do begin
     ;; Process
     rslt = hamspec_subbias_sngl( hamspec[arcs[qq]].img_root, /CLOBBER, $
                                  USEBIAS=usebias, BADROW=badrow, OVIMG=ovimg, $
                                  CBIN=cbin, RBIN=rbin, FRAME=frame, NAMP=hamspec[arcs[qq]].amp)
     sz = size(ovimg, /dimen)
     if qq EQ 0 and narc GT 1 then begin 
        msk = intarr(sz[0], sz[1], narc)
        msk[*] = 1
        sv_arc = fltarr(sz[0], sz[1], narc)
     endif 
     if narc GT 1 then begin 
        ;; Bad pix
        satur = where( ovimg GT 40000., nsat)
        tmpmsk = intarr(sz[0],sz[1])
        tmpmsk[*] = 1
        if nsat GT 0 then tmpmsk[satur] = 0
        msk[*,*,qq] = tmpmsk
        ;; Frame
        sv_arc[*,*,qq] = ovimg
     endif 
  endfor

  ;; Sum
  if narc EQ 1 then fin_arc = ovimg else begin
     ;; 
     fin_arc = total( sv_arc * msk, 3)
     fin_msk = total( msk, 3)
     ;; All rejected?
     allrej = where( total(msk, 3) LE 0, nallr, complement=gdm)
     if nallr GT 0 then begin
        tmp = sv_arc[*, *, narc-1]
        fin_arc[allrej] = tmp[allrej]
     endif
     ;; Normalize a bit
     fin_arc[gdm] = fin_arc[gdm] / fin_msk[gdm]
  endelse
      

  ;; Output
  print, 'hamspec_echmkarc:  Creating ', outfil
  mwrfits, fin_arc, outfil, head, /create, /silent
  ;mwrfits, fin_var, outfil, /silent

  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

  return
end

