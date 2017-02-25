;+ 
; NAME:
; uves_nrm2dspec
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'uves_*.fits' and passes back the uves structure.
;
; CALLING SEQUENCE:
;   
;  hires = hires_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./hires_*fits*)
;
; RETURNS:
;    hires -  MIKE structure
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
;   hires = hires_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   hires_rslvall
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_nrm2dspec, uves, setup, obj_id, BLAZE=blaze

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_nrm2dspec, uves, setup, obj, /BLAZE [v1.1]'
      return
  endif 

  if not keyword_set(OBJ_NM) then obj_nm = 'a'

  ;; Dat fil
  allexp = where((uves.type EQ 'OBJ' OR uves.type EQ 'STD') AND $
                 uves.flg_anly NE 0 AND $
                 uves.setup EQ setup AND uves.obj_id EQ obj_id )
  subfil = strcompress(strtrim(uves[allexp[0]].Obj,2), $
                       /remove_all)+obj_nm
  clrc = strtrim(round(uves[allexp[0]].xdangl),2)
  datfil = 'FSpec/'+subfil+'_'+clrc+'.fits'

  spec = xmrdfits(datfil,1,/silent)

  ;; Conti fil
  contifil = 'FSpec/'+subfil+'_'+clrc+'c.fits'
  if x_chkfil(contifil+'*') NE 1 then begin
      print, 'uves_nrm2dspec: No file named ', contifil
      stop
  endif
  conti = xmrdfits(contifil, /silent)
      
  ;; Normalize
  spec.fx = spec.fx / conti
  spec.var = spec.var / conti / conti
  
  ;; Write
  newout = 'FSpec/'+subfil+'_'+clrc+'_2DF.fits'
  mwrfits, spec, newout, /create
  spawn, 'gzip -f '+newout
  print, 'uves_nrm2dspec: Wrote ', newout

  return 

end
      
