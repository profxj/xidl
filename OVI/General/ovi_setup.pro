;+ 
; NAME:
; ovi_setup   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of WFCCD
;    spectroscopic images
;
; CALLING SEQUENCE:
;  ovi_setup, struct, phot_fil, FSPEC_FIL=, IMG_FIL=,$
;              OUTFIL=, FLG_PHOT=, SURV_FIL=
;
; INPUTS:
;  phot_fil -- File containing the photometry of the galaxies
;
; RETURNS:
;  struct -- A galaxy survey structure containing all of the relevant
;            info for the set of observations
;
; OUTPUTS:
;  OUTFIL=  -- Name of galaxy structure output file 
;              [default: 'ovi_strct.fits']
;
; OPTIONAL KEYWORDS:
;  SURV_FIL=  -- List of ID, RA and DEC for the galaxies
;  flg_phot=  -- Type of survey (1= LCO)
;  fspec_fil= -- Galaxy spectroscopic FITS file
;  IMG_FIL=   -- List of image files 
;  /WFCOEFF -- Setup the galaxy coefficients from the SDSS fits
;  ID_KLUDGE=  -- Offset ID values by +1 for all IDs larger than this
;                 input value (Kludge for 3C273)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ovi_setup, nght1_strct, 'WFTek5', 'LCO-100', /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ovi_setup, struct, phot_fil, FSPEC_FIL=fspec_fil, IMG_FIL=img_fil,$
               OUTFIL=outfil, FLG_PHOT=flg_phot, SURV_FIL=surv_fil, $
               WFCOEFF=wfcoeff, ID_KLUDGE=id_kludge

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'ovi_setup, struct, phot_fil, FSPEC_FIL=, IMG_FIL= OUTFIL=, FLG_PHOT=, SURV_FIL=, ID_KLUDGE= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OUTFIL ) then outfil = 'ovi_strct.fits'
  if not keyword_set( FLG_PHOT ) then flg_phot = 1L
  
; Read Phot file and create struct

  tmpstrct = { galsurveystrct }
  tmpstrct.field = ' '
  tmpstrct.obj_id = ' '
  tmpstrct.filter[0] = ' '
  tmpstrct.img_fil[0] = ' '
  tmpstrct.gal_type = ' '
  tmpstrct.fspec_fil[0] = ' '
  case flg_phot of 
      1: begin  ; LCO low z
          readcol, phot_fil, id, xpix, ypix, B, Bs, R, Rs, area, karea, star, $
            FORMAT='L,F,F,F,F,F,F,F,F,F'
          nobj = n_elements(ID)
          struct = replicate(tmpstrct, nobj)
          struct.id = id
          if keyword_set(ID_KLUDGE) then begin
              ;; 
              print, 'ovi_setup: Kludging the ID values!!'
              kludge = where(id GE ID_KLUDGE)
              struct[kludge].id = struct[kludge].id - 1
          endif
          struct.obj_id = 'a'
          struct.flg_anly = 1
          struct.mag[0] = B
          struct.mag[1] = R
          struct.magerr[0] = Bs
          struct.magerr[1] = Rs
          struct.filter[0] = 'B'
          struct.filter[1] = 'R'
          ;; xy
          struct.xypix[0] = xpix
          struct.xypix[1] = ypix
          ;; star, gal
          struct.stargal = star
          ;; Area
          struct.area = karea
          ;; Redshift
          struct.z = -1.
          nfilt = 2L
          all_nm = strtrim(struct.id,2)+'a'
      end
      else:
  endcase

; Image Fil

  if keyword_set( IMG_FIL ) then begin
      readcol, img_fil, filt, imgnm, FORMAT='A,A', /silent
      if n_elements(imgnm) NE nfilt then stop
      for q=0L,n_elements(filt)-1 do begin
          indx = where(struct[0].filter EQ strtrim(filt[q],2), nindx)
          if nindx EQ 1 then struct.img_fil[indx] = imgnm[q]
      endfor
  endif

; Fspec Fil

  if keyword_set( FSPEC_FIL ) then begin
      for i=0L,n_elements(fspec_fil)-1 do begin
          wfccd_wrfspec, fspec, fspec_fil[i], /read
          nspec = n_elements(fspec)
          for j=0L,nspec-1 do begin
              ;; Rej flg_anly = 0
              if fspec[j].flg_anly EQ 0 then continue
              ;; Name
              obj_nm = strtrim(fspec[j].slit_id,2)+strtrim(fspec[j].obj_id,2)
              indx = where(all_nm EQ obj_nm)
              if indx[0] EQ -1 then begin
                  all_nm = [all_nm, obj_nm]
                  struct = [struct, tmpstrct]
                  indx = n_elements(all_nm) - 1
                  struct[indx].id = fspec[j].slit_id
                  struct[indx].obj_id = fspec[j].obj_id
                  struct[indx].z = -1.
              endif else indx = indx[0]
              ;; Increment Flag
              if struct[indx].flg_anly MOD 4 LT 2 then $
                struct[indx].flg_anly = struct[indx].flg_anly + 2
              ;; Redshift
              struct[indx].z = fspec[j].zans.z
              ;; Fit coeff (WFCCD/SDSS)
              if keyword_set(WFCOEFF) then begin
                  nnrm = [0.985064d,  0.164177,  0.0492532,  0.0164177]
                  coeff = double(fspec[j].zans.theta[0:3]) * nnrm
                  coeff = coeff / sqrt(total(coeff^2))
                  struct[indx].gal_coeff[0:3] = coeff
                  ;; Early, late
                  struct[indx].gal_coeff[4] = $
                    total(coeff * [0.,1,-1,-1] ) / 1.5
              endif
              ;; Fspec
              a = where( strlen(struct[indx].fspec_fil) EQ 0)
              struct[indx].fspec_fil[a[0]] = fspec_fil[i]
          endfor
      endfor
  endif
  
; Increment redshift flag

  a = where(struct.z GT -0.0999 AND ((struct.flg_anly MOD 8) LT 4), na)
  if na NE 0 then struct[a].flg_anly = struct[a].flg_anly + 4

; Survey File
  if keyword_set(SURV_FIL) then begin
      readcol, surv_fil, id_surv, ra, dec, FORMAT='L,F,F', /silent
      ;; Loop
      for j=0L,n_elements(id_surv)-1 do begin
          a = where(struct.id EQ id_surv[j],na)
          if na NE 0 then begin
              struct[a].flg_survey = 1
              struct[a].ra = ra[j]
              struct[a].dec = dec[j]
          endif
      endfor
  endif
      
; Sort
  srt = sort(struct.id)

; Write the structure to FITS

  mwrfits, struct[srt], outfil, /create

; All done
  print, 'ovi_setup: All done!'

end
