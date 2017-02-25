;+ 
; NAME:
; hires_calibstd   
;   Version 1.1
;
; PURPOSE:
;    Create a sensitivity function given a standard star and its
;    appropriate calibration file.  The program mainly does some
;    bookkeeping and then calls x_echfitstd
;
; CALLING SEQUENCE:
; hires_calibstd, hires, indx, HSTFIL=, CHKFIT=, ESOFIL=
;   
; INPUTS:
;   hires   -  MIKE structure
;   indx   -  Index of the Standard star in the MIKE structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Sensitivity file in the 'Extract/' directory
;
; OPTIONAL KEYWORDS:
;  OUTFIL - Name of the output file (default:
;           'Extract/Sens_{hiresroot}.fits' )
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_calibstd, hires, 1, ESOFIL='fhr4469.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_echfitstd
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_calibstd, hires, indx, HSTFIL=hstfil, CHKFIT=chkfit, $
                    ESOFIL=esofil,  OBJFIL=objfil, $
                    SV_ALLFIT=sv_allfit

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_calibstd, hires, indx, HSTFIL=, FXFIL=, NCOEFF= [v1.1]'
    return
  endif 

;  Optional Keywords
  if n_elements(indx) NE 1 then stop
  if not keyword_set( OBJFIL ) then objfil = hires[indx].obj_fil
  if not keyword_set( OUTFIL ) then begin
      slen = strlen(hires[indx].img_root)
      outfil = hires_getfil('sens_fil', chip=hires[indx].chip, $
                            frame=hires[indx].frame, /name)
;      outfil = 'Extract/Sens_'+hires[indx].img_root
  endif

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
     if strlen(HSTFIL) LT 30 then $
         hstfil = getenv('XIDL_DIR')+'/Spec/Flux/HSTFIL/'+hstfil
       readcol, hstfil, swv, sfx
  endif

  if keyword_set( ESOFIL ) then begin
      if strlen(ESOFIL) LT 30 then $
        esofil = getenv('XIDL_DIR')+'/Spec/Flux/ESOFIL/'+esofil
      readcol, esofil, swv, sfx
;
;     let's keep it at 1d-16
;      sfx = sfx*1d-16
  endif

  objstr = xmrdfits(objfil, 1, STRUCTYP='hiresobjstrct', /silent)
; Create summed array
  
;  velpix = (hires[indx].side EQ 1 ? 1.50d : 2.10d) * hires[indx].rowbin
;  loglam = alog10(1.0d + velpix / 299792.458d)
;  wave0  = alog10(3000.0d)
;  npix = 100000L
;  wvtot = 10^(wave0 + dindgen(npix)*loglam)
;  fxtot = dblarr(npix)

  ;; Added by JFH
  if not keyword_set(BOX) then begin
      objstr.box_wv  = objstr.wave
      objstr.box_fx  = objstr.fx
      objstr.box_var = objstr.var
  endif
  
  ;; Normalize exposure
  objstr.box_fx = objstr.fx / hires[indx].exp
  objstr.box_wv = objstr.wave
  objstr.box_var = objstr.var / (hires[indx].exp^2)
  objstr.sig = sqrt(objstr.var / (hires[indx].exp^2))
;  objstr.npix = objstr.nrow

  ;; Call the fit program
  x_echfitstd, objstr, swv, sfx, outfil
  

  print, 'hires_calibstd:  Writing ', outfil


  print, 'hires_calibstd: All Done!'
  return
end
