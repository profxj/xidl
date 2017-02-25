;+ 
; NAME:
; mike_chksky   
;     Version 1.1
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_mike0020.fits).  The program also outputs
;    information related to the sky fit in the Sky directory.  The
;    user has control over which order to subtract (ORDR=), the
;    wavelength image to use (AIMG=) and can input a break point file
;    for detailed sky line subtraction (SKLFIL=).
;
; CALLING SEQUENCE:
;   
;  mike_chksky, mike, setup, side, obj_id, [exp], /CHK, /STD, ORDR=, /NOVAC,
;  SKLFIL=, BORDR=, AIMG=
;
; INPUTS:
;   mike   -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Sky subtract standard star
;  /CHK     - Show steps along the way
;  AIMG=    - Use alternate Wavelength map for sky subtraction
;             (string filename)
;  ORDR=    - Sky subtract a single order only (e.g.  [5L,5L])
;  /NOVAC   - Do not perform vacuum wavelength correction
;  BORDR=   - Order to begin bspline fitting (default: 5L)
;  SKLFIL=  - ASCII file setting breakpoints around sky lines (string)
;  /CLOBBER - Overwrite any previos sky image
;  /USEOLD  - Overwrite only the new orders into the old sky sub image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_chksky, mike, 1L, [0L], /CHK, ORDR=7L   
;           {Sky sub exposure 0 and order 7 only}
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_chksky, mike, setup, obj_id, side, exp, STD=std, ORDR=ordr, $
                 NOVAC=novac, CCHI=CCHI, AIMG=aimg, $
                 USEOLD=useold, BCHK=bchk

;
  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'mike_chksky, mike, setup, obj_id, side, exp, /CCHI, /STD, '
      print, '          BORDR=, AIMG=, /USEOLD, /BCHK [v1.0]'
      return
  endif 
  
;  Optional Keywords

;  Setup and setup
  case side of 
      1: ostep = 1
      2: ostep = -1
      else: stop
  endcase

      
;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(mike.flg_anly NE 0 AND mike.side EQ side AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   strtrim(mike.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'mike_chksky: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse
      
; TRIM
;  if not keyword_set(TRIM) then trim = 4L / mike[indx[0]].colbin

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)


;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, side=side)
  if (where(tag_names(ordr_str) EQ 'PROFILE0'))[0] NE -1 then begin
     if total(ordr_str.profile0) GT 0 then profile_info = 1
  endif

  nordr = n_elements(ordr_str)

; Orders
  if not keyword_set( ORDR ) then begin  ;; Orders to sky subtract
      ordr = [ordr_str[0].order,ordr_str[nordr-1].order]
      flg_ordr = 0
  endif else begin
      if keyword_set( USEOLD ) then flg_ordr = 1 else flg_ordr = 0
      if n_elements(ordr) NE 2 then begin
          print, 'mike_chksky: ORDR must be a 2 element array'
          return
      endif
  endelse

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      ;;  Read Arc
      if not keyword_set(AIMG) then arc_img = mike[indx[exp[q]]].arc_img $
      else arc_img = aimg
      ;; READ
      if x_chkfil(arc_img+'*') EQ 0 then begin
          print, 'mike_chksky: Arc file doesnt exist!', arc_img
          return
      endif
      print, 'mike_chksky: Reading Arc file ', arc_img
      img_arc = xmrdfits(arc_img, /silent) 
      sz_arc = size(img_arc, /dimensions)

      ;; xoff
      xoff = mike[indx[exp[q]]].arc_xyoff[0]

;;;;;;;   Big Change HERE, let's just add xoff to the ordr_str params...  ;;;;;;;;;;;;;;;;;;;;;;;;;
      ordr_str.lhedg = ordr_str.lhedg + xoff
      ordr_str.rhedg = ordr_str.rhedg + xoff
      ordr_str.lhc[0] = ordr_str.lhc[0] + xoff
      ordr_str.rhc[0] = ordr_str.rhc[0] + xoff


      ;; Open Obj file
      objfil = mike[indx[exp[q]]].obj_fil
      a = findfile(objfil+'*', count=na)
      if na EQ 0 then begin
          print, 'mike_chksky: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)

      ;; IMG+VAR Fil 
      imgfil = strtrim(mike[indx[exp[q]]].img_final,2)
      skyfil = 'Sky/sky_'+mike[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'mike_chksky: Image file doesnt exist!', imgfil
          return
      endif
      print, 'mike_chksky: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Mask
      msk = lonarr(sz_img[0],sz_img[1])
      maskimage = mike_ordermask(sz_img[0], sz_img[1], ordr_str, trim=0.)

      ;; Inter-order light fit
      if keyword_set( profile_info ) then begin
          print, 'Fitting inter-order light with mike_fitgap  5x5 bkpts'
          scattered_light = mike_fitgap(img, ivar, maskimage, nxbkpt=5, nybkpt=5)
          img_sub = img - scattered_light
      endif else img_sub = img

;
;
      obj_temp = ordr_str
      obj_temp.lhedg = objstr.trace[0:sz_img[1]-1] - $
                     replicate(1,sz_img[1]) # objstr.aper[0]
      obj_temp.rhedg = objstr.trace[0:sz_img[1]-1] + $
                     replicate(1,sz_img[1]) # objstr.aper[1]

      mask_temp_obj = mike_ordermask(sz_img[0], sz_img[1], obj_temp, trim=0.)
;
;      Object pixels are order number + 10000
;
;
      maskimage = maskimage + (mask_temp_obj GT 0)*10000L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop on Orders
      for qq=ordr[0],ordr[1],ostep do begin
          ;;
          mm = where(ordr_str.order EQ qq, nmm)
          if nmm EQ 0 then stop else mm = mm[0]
          skystr = xmrdfits(skyfil, mm+1, /silent)

          ;; 
          print, 'mike_chksky: Checking order ', string(qq, FORMAT='(i3)')

          if keyword_set(profile_info) then begin

                  
              inorder = where(maskimage mod 10000L EQ qq AND img_arc GT 0)
              wave_sort = inorder[sort(img_arc[inorder])]
              xstart = wave_sort mod sz_img[0]
              ystart = wave_sort /   sz_img[0]
              slit_length = ordr_str[mm].rhedg[sz_img[1]/2] - $
                ordr_str[mm].lhedg[sz_img[1]/2]
              ordrcen  =  (ordr_str[mm].lhedg[ystart] + $
                           ordr_str[mm].rhedg[ystart])/2.0 ; Offset
              slit_frac  = 2.0* (xstart- ordrcen) / slit_length
              
              
              profile = mike_profile_return(xstart, ystart, ordr_str[mm])
              
              img_norm = (img_sub)[wave_sort] / $
                (profile + (profile EQ 0)) * (profile GT 0)
              
              img_norm_ivar = ivar[wave_sort] * profile^2 * (profile GT 0)
              wave_norm =  img_arc[wave_sort]
              
              sky_pix = (maskimage[wave_sort] LT 200 AND profile GT 0.3)
              sky = where(sky_pix)
              obj_pix = (maskimage[wave_sort] GT 10000L)
              
              nsky_pix = total(sky_pix)
              nobj_pix = total(obj_pix)
              
              if nsky_pix LT 0.5*nobj_pix then $
                print, 'mike_chksky: sky pixel to object pixel ratio is lower than 0.5'
                
              if nobj_pix LT 2*sz_img[1] then $
                message, 'maybe a short order '+string(nobj_pix)
              
              bkpt_s = 0    
              
              ;; calculate chi^2 in sky pixels
              
              img_fit = bspline_valu(wave_norm, skystr)
              img_model = img_fit*profile+scattered_light[wave_sort]
              chi = (img[wave_sort]-img_model)*  sqrt(ivar[wave_sort])
;              badsky = where(pixmsk EQ 0)
                
;
;   Calculate first moment of chi^2 in slit_frac
;
              numer = slit_frac[sky] * ((chi[sky] > (-5.0)) < 5)
              guess_pix_shift = total(numer)/total(abs(chi) < 5) $
                * (slit_length/2.0)
              
              use_shift = (guess_pix_shift > (-2.0)) < 2.0
              print, 'The order residuals suggest a shift relative to ' + $
                'the order trace of ', guess_pix_shift, $
                ', Using a shift of ', use_shift,' pixels to center search.'
                

             if keyword_set( BCHK ) then begin
                 x_splot, 10^wave_norm[sky], img[sky], $
                   ytwo=profile[sky]*bspline_valu(wave_norm[sky], skystr)+$
                   scattered_light[sky], /block
             endif
             
             if keyword_set( CCHI ) then begin
                 if skystr.badsky[0] NE -1 then $
                   x_splot, 10^wave_norm[sky], chi[sky], psym1=1, /block, $
                   xtwo=10^wave_norm[sky[skystr.badsky]], $
                   ytwo=chi[sky[skystr.badsky]], $
                   psym2=2 $
                 else x_splot, 10^wave_norm[sky], chi[sky], psym1=3, /block
             endif
             
;             img_new[wave_sort] = img[wave_sort] - img_model
         endif else begin
             stop
         endelse
             
     endfor
 endfor
  
;  DONE
  print, 'mike_skysub: All done! '
  return
end
