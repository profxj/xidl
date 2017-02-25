 ;+ 
; NAME:
; hires_nrm2dspec
;     Version 1.1
;
; PURPOSE:
;   After one has produced a combined spectrum with hires_combspec and
;   a continuum with hires_conti2d, this code will normalize the data.
;
; CALLING SEQUENCE:
;  hires_nrm2dspec, hires, setup, obj_id, chip, BLAZE=
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
;
; RETURNS:
;
; OUTPUTS:
;  Normalized data
;
; OPTIONAL KEYWORDS:
;  /BLAZE -- Use a blaze normalized file (bz.fits extension)
;  ZERO_FIL=  -- Array of file names that store the zero level
;                information.  This zero level is subtracted prior to
;                normalization
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_nrm2dspec, hires, setup, obj_id, chip, BLAZE=blaze $
                     , OBJ_NM = OBJ_NM, FIL_NM=fil_nm, OUT_FIL=out_fil, $
                     CONTI_FIL=conti_fil, ZERO_FIL=zero_fil

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_nrm2dspec, hires, setup, obj, [chip], /BLAZE, ZERO_FIL= [v1.1]'
      return
  endif 

  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(CHIP) then chip = [1,2,3]
  nchip = n_elements(chip)

  if keyword_set(FIL_NM) then begin
      if not keyword_set(CONTI_FIL) then begin
          print, 'hires_nrm2dspec: Must set conti_fil!'
          return
      endif
      if not keyword_set(OUT_FIL) then begin
          print, 'hires_nrm2dspec: Must set out_fil!'
          return
      endif
      nchip = 1
  endif

  for qq=0L,nchip-1 do begin
      ;; Dat fil
      allexp = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') AND $
                     hires.flg_anly NE 0 AND $
                     hires.setup EQ setup AND hires.obj_id EQ obj_id $
                     AND hires.chip EQ chip[qq])
      subfil = strcompress(strtrim(hires[allexp[0]].Obj, 2)+obj_nm, /rem)
      case chip[qq] of
          -1: clrc = '_S' 
          1: clrc = '_B' 
          2: clrc = '_G' 
          3: clrc = '_R' 
          else: stop
      endcase

      ;; Data
      if keyword_set(FIL_NM) then datfil = fil_NM else begin
          if not keyword_set(BLAZE) then datfil = 'FSpec/'+subfil+clrc+'.fits' $
          else datfil = 'FSpec/'+subfil+clrc+'bz.fits'
      endelse
      spec = xmrdfits(datfil,1,/silent)

      ;; Zero Level
      if keyword_set(ZERO_FIL) then begin
         if strlen(zero_fil[qq]) GT 0 then begin
            zero_lvl = xmrdfits(zero_fil[qq]) 
            ;; Subtract
            print, 'hires_nrm2dspec: Subtracting zero level: ', zero_fil[qq]
            spec.fx = spec.fx - zero_lvl
         endif
      endif

      ;; Conti fil
      if keyword_set(CONTI_FIL) then contifil = conti_FIL else $
        contifil = 'FSpec/'+subfil+clrc+'c.fits'
      if x_chkfil(contifil+'*') NE 1 then begin
          print, 'hires_nrm2dspec: No file named ', contifil
          stop
      endif
      conti = xmrdfits(contifil, /silent)
      
      ;; Normalize
      spec.fx = spec.fx / conti
      spec.var = spec.var / conti / conti
      IF TAG_EXIST(spec, 'NOVAR') THEN spec.novar = spec.novar/conti/conti
      
      ;; Write
      if keyword_set(OUT_FIL) then newout = OUT_FIL else $
        newout = 'FSpec/'+subfil+clrc+'F.fits'
      mwrfits, spec, newout, /create
      print, 'hires_nrm2dspec: Wrote ', newout
  endfor


  return 

end
      
