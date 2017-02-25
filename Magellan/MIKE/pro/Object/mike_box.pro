;+ 
; NAME:
; mike_box
;     Version 1.1
;
; PURPOSE:
;   Extract 1D spectra from the 2D images.  For each order, a boxcar
;   and an optimal extraction is performed.  For the latter, an object
;   profile is derived and both the object flux and sky are fit throughout
;   the order (i.e. not row by row).  The main driver is
;   mike_box_sngl.  For the optimal extraction, the data is extracted
;   to a specific set of vacuum wavelengths, chosen to be the same for
;   every spectrum to facilitate coadding without rebinning.  Here are
;   the steps in detail:
;
;    1.  Perform a boxcar extraction using extract_box
;    2.  Estimate the SNR per order from the boxcar extraction
;      -- LOOP ON ORDERS IN DECREASING SNR --
;    3.  Fit the boxcar extraction with a bspline
;    4.  Calculate the object profile
;          a. bspline_iterfit the flux vs position on slit
;          b. Force the profile to be positive everywhere and have a
;          sensible FWHM
;    5.  Fit the order using the profile and sky (bspline_extract)
;
; CALLING SEQUENCE:
;   
;  mike_box, mike, setup, obj_id, side, [exp], /RESCHK, /CHK
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) or Red (2) side
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   Fills up the 1D spectral tags in the Object structure.
;
; OPTIONAL KEYWORDS:
;   /OCHK   - Plot the extracted flux (optimal) for each order
;   /CHK    - Live plots of the profile and trace
;   /RESCHK - Check the residuals of the 2D, fully extracted image
;   /STD    - Extraction should be set for a standard star
;   /DEBUG  - Stop within extraction routine to check stuff
;   /SKIPSKYSUB - Perform extraction on the processed but NOT sky
;                 subtracted image.
;   HIGHSNR - Value of SNR^2 of the data for a given order which when
;             exceeded mike_box uses an additional parameter for the
;             profile shape.  (Default:  500 corresponding to SNR=22)
;             Lowering this parameter may improve extraction.
;   ORDRS   - Input array of physical order numbers to extract
; aasdf
; Optional OUTPUTS:
;
; COMMENTS:
;  The program extracts the orders in order of decreasing SNR.  If the
;  SNR is lower than lowsnr (default: 2.49) then the optimal
;  extraction is performed using the profile parameters from the
;  previous order(s).
;
; EXAMPLES:
;   mike_box, mike, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_box_sngl
;  extract_boxcar
;  smooth_mask
;  bspline_extract
;
; REVISION HISTORY:
;   26-Aug-2003 Written by SMB
;-
;------------------------------------------------------------------------------


function smooth_mask, outmask, indx, ncol
   
     if N_PARAMS() LT 3 then message, 'smooth_mask(outmask, indx, ncol)' 
     x = indx mod ncol
     y = indx /  ncol

     nx = max(x) - min(x) + 5L
     ny = max(y) - min(y) + 5L
     new_indx = (y-min(y)+2) *nx + x - min(x) + 2 
     
     temp_arr  = fltarr(nx,ny)
     temp_arr[new_indx] = 1.0 * (outmask EQ 0)
    
     kernel = [[0,0,1,0,0],[0,1,1,1,0],[1,1,1,1,1],[0,1,1,1,0],[0,0,1,0,0]] 
     smooth_arr = convol(temp_arr,kernel)
     
     return_mask = smooth_arr[new_indx] EQ 0
   
return, return_mask
end 
;---------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_box, mike, setup, obj_id, side, exp, STD=std, $
              ORDRS=ordrs, HIGHSNR=highsnr, skipskysub=skipskysub, NOVAC=novac, $
              readsky=readsky, _EXTRA=EXTRA, MODELOUT=modelout, NOWAVE=nowave, $
              NOHELIO=nohelio

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_box, mike, setup, obj_id, side, [exp], ' + $
        '/DEBUG, /CHK, /OCHK, /RESCHK, /READSKY, ORDRS=,'
      print, '          /STD, HIGHSNR= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( REJSIG ) then rejsig = 7.
  if not keyword_set( uniformsky ) then uniformsky = 0
  if not keyword_set( FIN_TRC ) then fin_trc = 2

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   mike.side EQ side AND mike.obj_id EQ obj_id AND $
                  (strtrim(mike.type,2) EQ 'OBJ' OR $ 
                   strtrim(mike.type,2) EQ 'STD'), nindx)
  
      if nindx EQ 0 then begin
          print, 'mike_boxextrct: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STD star
;      stop
      indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   mike.side EQ side AND $
                   mike.obj_id EQ obj_id AND strtrim(mike.type,2) EQ 'STD', nindx)
      if nindx EQ 0 then begin
          print, 'mike_box: No standards to find for obj_id!', obj_id
          return
      endif
      radius = 40L
  endelse

;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, side=side)
  nordr = n_elements(ordr_str)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop

  for q=0L,n_elements(exp)-1 do begin

      ;; Offset edges
      xyoff = mike[indx[exp[q]]].arc_xyoff

      ;; Shift ORDR_STR
      ordr_shift = ordr_str
      shft = mike_shifti(xyoff, OSTR=ordr_shift)
      rnd_edg = [[round(ordr_shift.lhedg)], [round(ordr_shift.rhedg)]]
      sz = size(rnd_edg, /dimensions)
      rnd_edg = reform(rnd_edg,sz[0],nordr,2)

      print, 'mike_boxextrct: Reading files...'

      ;;;;;;;;;;;;;;
      ;; Open Obj file
      objfil = mike[indx[exp[q]]].obj_fil
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'mike_box: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
    
;      nobj = n_elements(objstr)

      ;;;;;;;;;;;;;;
      ;; SKY SUB Fil 
      imgfil = mike_getfil('fin_fil', setup, SUBFIL=mike[indx[exp[q]]].img_root, $
                           /name)
;      imgfil = mike[indx[exp[q]]].img_final
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'mike_box: No Image file!  Returning...'
          return
      endif
      print, 'mike_box: Image file -- ', imgfil
      head = xheadfits(imgfil)
      img = xmrdfits(imgfil, 0, /silent)   ; non-SKY subtracted
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(ivar, /dimensions)

      scatlight_file = mike_getfil('scatt_fil',setup, $
                     subfil=mike[indx[exp[q]]].img_root, /name, chkfil=chkfil)

      scattered_light = 0

      if chkfil EQ 1 then scattered_light = mike_getfil('scatt_fil',setup, $
                                subfil=mike[indx[exp[q]]].img_root) $
      else begin
          if not keyword_set( SCATT_SKIP ) then begin
              print, 'mike_box: Fitting inter-order light ' + $
                'with x_fitgap  8x10 bkpts'
              maskimage = x_ordermask(sz_img[0], sz_img[1], $
                                         ordr_shift, trim=0.)
              
              ;; Inter-order light fit
              scattered_light = x_fitgap(img, ivar, maskimage,$
                                            nxbkpt=8,nybkpt=10)
          endif
      endelse

      img_sub = img - scattered_light

      if keyword_set(skipskysub) then skysub = img_sub $
      else skysub = xmrdfits(imgfil, 2, /silent) ; SKY subtracted

      ;; Need var for most of my routines
;      var = 1./(ivar + (ivar EQ 0)) * (ivar GT 0)
      
      if (size(img))[0] EQ 0 then begin
        print, 'mike_box: No skysubtracted image found'
        print, 'mike_box: Either run mike_skysub or mike_box ' + $
          'with the /skipsky option'
        continue
      endif

      ;;  Read Arc
      img_arc = 0
      hel_corr = 1.0d
      helio = 0.
      if NOT keyword_set(NOWAVE) then begin
        arc_img = strtrim(mike[indx[exp[q]]].arc_img,2)
        print, 'mike_box: Arc -- ', arc_img
        if x_chkfil(arc_img+'*') EQ 0 then begin
            print, 'mike_boxextrct: No Arc file!  Returning...', arc_img
            return
        endif
        img_arc = xmrdfits(arc_img, /silent) 

        ;; Vacuum
        if not keyword_set( NOVAC ) then begin
            print, 'mike_box: Converting to vacuum wavelengths'
            a = where(img_arc GT 0.)
            tmpaimg = 10^img_arc[a]
            airtovac, tmpaimg
            img_arc[a] = temporary(alog10(tmpaimg))
        endif

        ;;  HELIO correction
        if not keyword_set(NOHELIO) then begin
            radeg = sxpar(head,'RA-D')
            decdeg = sxpar(head,'DEC-D')
            helio = (-1.) * x_keckhelio(radeg, decdeg, $
                                        mike[indx[exp[q]]].equinox, $
                                        jd=mike[indx[exp[q]]].date, $
                                        altitude=2282., $
                                        longitude=360.-70.70, latitude=-29.08333)
            print, 'mike_box: heliocentric correction :', helio,' km/s' , $
                   format='(a,f8.3,a)'
            hel_corr = sqrt( (1.d + helio/299792.458d) / (1.d - helio/299792.458d) )
            img_arc = img_arc + alog10(hel_corr)
        endif else hel_corr = 1.d
      endif

      ;; Check for Original sky fit
      if keyword_set( readsky ) then begin
;         skyfil = 'Sky/sky_'+mike[indx[exp[q]]].img_root
          skyfil = mike_getfil('sky_fil', setup, $
                               SUBFIL=mike[indx[exp[q]]].img_root, /name)
;         skyfil = findfile(skyfil)
      endif         

      ;; 
      print, systime()
      mike_slit = 5.  ; arcsec

      ;; Velocity pix
      if keyword_set(img_arc) then begin
          sz = size(img_arc, /dimens)
          nrow = sz[1]
          velpix = (side EQ 1 ? 1.50d : 2.10d) * round(4096./nrow)
      endif else begin
          stop
      endelse

      x_extechopt, img_sub, skysub, ivar, ordr_shift, objstr, velpix, $
        chk=chk, img_arc=img_arc, HIGHSNR=highsnr,$
        ORDRS=ordrs, skyfil=skyfil, $
        helio=alog10(hel_corr), OBJ_NAME=mike[indx[exp[q]]].Obj, $
        _EXTRA=EXTRA, ORDERMASK = ordermask, SLIT_LEN=mike_slit, $
        MODEL_OBJ=model_obj, MODEL_SKY=model_sky, MODEL_PROF=model_prof, $
        FIN_TRC=fin_trc

;      mike_box_sngl, img_sub, skysub, ivar, ordr_shift, objstr, chk=chk, $
;        img_arc=img_arc, SIDE=side, HIGHSNR=highsnr,$
;        SKIPSKYSUB=skipskysub, ORDRS=ordrs, skyfil=skyfil, $
;        helio=alog10(hel_corr), OBJ_NAME=mike[indx[exp[q]]].Obj, $
;        _EXTRA=EXTRA, ORDERMASK = ordermask, $
;        MODEL_OBJ=model_obj, MODEL_SKY=model_sky, MODEL_PROF=model_prof

      sxaddpar, head, 'HELIO', helio
;      xmodfits, imgfil, 0, head
      print, 'mike_box: Helio correction applied -- ', helio, hel_corr
      
      print, systime()
      ;; Ouptut Spectra
      print, 'mike_box: Output spectrum in -- ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

      ;; Output the image
      mike_taghead, head
      mwrfits, img, imgfil, head, /create
      mwrfits, ivar, imgfil
      mwrfits, skysub, imgfil
      mwrfits, img_sub - model_sky, imgfil
      spawn, 'gzip -f '+imgfil

      ;; Optional Image 
      if keyword_set(modelout) then begin
        modelfil = strmid(objfil,0,rstrpos(objfil,'.fits')) + '-Model.fits'
        mwrfits,  img_arc, modelfil, head, /create, /silent
        mwrfits,  ordermask, modelfil, /silent
        mwrfits,  model_prof, modelfil, /silent
        mwrfits,  model_obj, modelfil, /silent
        mwrfits,  model_sky, modelfil, /silent
        spawn, 'gzip -f '+modelfil
      endif

        
  endfor
  
;  DONE
  print, 'mike_box: All done! '
  return
end

