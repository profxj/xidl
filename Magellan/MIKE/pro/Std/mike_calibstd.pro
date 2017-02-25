;+ 
; NAME:
; mike_calibstd   
;   Version 1.1
;
; PURPOSE:
;    Create a sensitivity function given a standard star and its
;    appropriate calibration file.  
;
; CALLING SEQUENCE:
; mike_calibstd, mike, indx, HSTFIL=, CHKFIT=, ESOFIL=
;   
; INPUTS:
;   mike   -  MIKE structure
;   indx   -  Index of the Standard star in the MIKE structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Sensitivity file in the 'Extract/' directory
;
; OPTIONAL KEYWORDS:
;  OUTFIL - Name of the output file (default:
;           'Extract/Sens_{mikeroot}.fits' )
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_calibstd, mike, 1, ESOFIL='fhr4469.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_echfitstd
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_calibstd, mike, indx, HSTFIL=hstfil, CHKFIT=chkfit, $
                BSPLIN=bsplin, SWV=swv, SFX=sfx, ESOFIL=esofil, $
                BSET=bset, YFIT=yfit, SENS=sens, EVERYN=everyn, $
                OBJFIL=objfil, indx2=indx2, INTER=INTER, NCOEFF=ncoeff, $
                   SV_ALLFIT=sv_allfit

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_calibstd, mike, indx, HSTFIL=, FXFIL=, NCOEFF= [v1.1]'
    return
  endif 

;  Optional Keywords
  if n_elements(indx) NE 1 then stop
  if not keyword_set( EVERYN ) then everyn = 36/mike[indx].rowbin
  if not keyword_set( OBJFIL ) then objfil = mike[indx].obj_fil
  if not keyword_set( OUTFIL ) then begin
      slen = strlen(mike[indx].img_root)
      outfil = mike_getfil('sens_fil', SUBFIL=mike[indx].img_root, /name)
;      outfil = 'Extract/Sens_'+mike[indx].img_root
  endif
  if not keyword_set( NCOEFF ) then ncoeff= 7

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
     if strlen(HSTFIL) LT 30 then $
         hstfil = getenv('MIKE_DIR')+'/pro/Std/HSTFIL/'+hstfil
       readcol, hstfil, swv, sfx
  endif

  if keyword_set( ESOFIL ) then begin
      if strlen(ESOFIL) LT 30 then $
        esofil = getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+esofil
      readcol, esofil, swv, sfx
;
;     let's keep it at 1d-16
;      sfx = sfx*1d-16
  endif

  objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)

; Create summed array
  
;  velpix = (mike[indx].side EQ 1 ? 1.50d : 2.10d) * mike[indx].rowbin
;  loglam = alog10(1.0d + velpix / 299792.458d)
;  wave0  = alog10(3000.0d)
;  npix = 100000L
;  wvtot = 10^(wave0 + dindgen(npix)*loglam)
;  fxtot = dblarr(npix)

  ;; Normalize exposure
  objstr.box_fx = objstr.box_fx / mike[indx].exp

  ;; Call the fit program
  x_echfitstd, objstr, swv, sfx, outfil
  

  print, 'mike_calibstd:  Writing ', outfil


  print, 'mike_calibstd: All Done!'
  return
end
