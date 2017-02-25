;+ 
; NAME:
; kast_calibstd   
;   Version 1.1.5
;
; PURPOSE:
;    Creates a sensitivity function from a standard star; sensfunc can
;    be used by kast_flux to flux-calibrate science-object spectra
;
; CALLING SEQUENCE:
;  kast_calibstd, kast, indx, outfil, STFIL=, STTYPE=
;
; INPUTS:
;   kast  -- Kast IDL structure
;   indx  -- Index of standard star (not obj_id)
;
; REQUIRED KEYWORDS:
;   STFIL= --  Spectrophotometric file for the standard star to be
;             processed
;          --  Units are F_lambda vs wavelength
;   STTYPE= --  Type of standard file - tells the program whether to
;              read it using xmrdfits (STTYPE=fits, vectors are
;              st.wavelength and st.flux) or readcol (STTYPE=ascii,
;              columns are wavelength, f_lambda)
;
; RETURNS:
;
; OUTPUTS:
;   outfil -- Name of sensitivity function (written with mwrfits,
;             readable with xmrdfits)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_calibstd, kast, 83, 'Calibs/sensfunc_test',
;   stfil='kastsens_g191g1d46', sttype=fits
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Mar-2003 Written by JXP
;   25-Jan-2005 editing by SLM
;   12-Nov-2005 revert code to use indx correctly, KLC
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_calibstd, kast, indx, outfil, STFIL=stfil, STTYPE=sttype

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'kast_calibstd, kast, indx, outfil, STFIL=, STTYPE= [v1.1.5]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( EVERYN ) then everyn = 5

; Open objfil
;filno=indx+1

;objfil='Extract/Obj_b'+strtrim(filno,1)+'.ccd.gz'
  objfil = kast[indx[0]].obj_fil 
; if x_chkfil(objfil+'*') EQ 0 then begin
;      print, 'kast_calibstd: No Obj file ', objfil
;      return
;  endif

;stop

  objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
  npix = objstr[0].npix

; Grab standard file
  ;; Choose data type
  if(STTYPE eq 'fits') then begin

      stfil=strtrim(stfil,2)+'.fits'

      ;;Read file
      st = xmrdfits(stfil, 1, /silent)
      wv = st.wavelength
      fx = st.flux
  endif else begin

      if(STTYPE eq 'ascii') then begin

          stfil=strtrim(stfil,2)+'.dat'
          readcol,stfil,F='X,D,D',wv,fx
          ;readcol,stfil,F='F,D',wv,fx
      endif
  endelse


;  if keyword_set(STFIL) then begin
;      ;; Read file
;      st = xmrdfits(stfil, 1, /silent)
;      ;; 
;      wv = st.wavelength
;      fx = st.flux
;  endif
  
  if not keyword_set(fx) then stop

  ;; BSpline
  bset = bspline_iterfit(wv,fx, yfit=yfit, everyn=everyn)
  
  ;; Calculate at wavelength
  sens = bspline_valu(objstr[0].wave[0:npix-1], bset) / $
    (objstr[0].fx[0:npix-1] / kast[indx[0]].exp)
  
  ;; BSpline
  bset = bspline_iterfit(objstr[0].wave[0:npix-1], sens, $
                         yfit=yfit, everyn=everyn)
  
  x_splot, sens, ytwo=yfit, /block
  
  ;; Output
  mwrfits, bset, outfil, /create
      
  print, 'kast_calibstd: All Done!'
  return
end
