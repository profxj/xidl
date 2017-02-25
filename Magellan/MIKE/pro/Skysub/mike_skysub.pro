;+ 
; NAME:
; mike_skysub   
;     Version 1.2
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_mike0020.fits).  The program also outputs
;    information related to the sky fit in the Sky directory.  The
;    user has control over which order to subtract (ORDR=), the
;    wavelength image to use (AIMG=).  The program preforms the
;    following steps:
;
;    1.  Offset the Order structure by the appropriate shift
;    2.  Subtract scattered light (mike_fitgap)
;     -- LOOP ON ORDERS --
;    3.  After masking the object, the code finds the sky pixels and
;    normalizes them by the slit profile.
;    4.  It bspline fit to the sky
;    5.  It assesses the fit and adds in breakpoints in a somewhat
;    smart way.
;    6.  It performs a second bspline fit
;    7.  The coefficients of the fit are saved and the image is
;    subtracted 
;    8.  The code compares the extracted sky spectrum against the UVES
;    sky spectrum to guess at a shift.  This is saved to the Object
;    structure. 
;     -- END OF LOOP
;    9.  The final sky subtracted image is written to the fits file
;    (extension #2)
;
; CALLING SEQUENCE:
;   
;  mike_skysub, mike, setup, obj_id, side, [exp], /CHK, /STD, ORDR=,
;  AIMG=, /USEOLD, /FCHK
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
;   A sky subtracted image added as extension #2 to the fits file.
;   A file in the directory Sky/ describing the fit
;
;
; OPTIONAL KEYWORDS:
;  /STD     - Sky subtract standard star
;  /CHK     - Show steps along the way
;  AIMG=    - Use alternate Wavelength map for sky subtraction
;             (string filename)
;  ORDR=    - Sky subtract a single order only (e.g.  [5L,5L])
;  /CLOBBER - Overwrite any previos sky image
;  /USEOLD  - Overwrite only the new orders into the old sky sub image
;  /FCHK    - Show the final sky subtracted image
;  /DOALL   - ??
;
; OPTIONAL OUTPUTS:
;  SCATTERED_LIGHT= -- Scattered light image.  It may be useful to
;                      save and use this if you are subtracting individual
;                      orders
;
; COMMENTS:
;   Note added by J. Moustakas: If /STD is used, then the code expects
;     the **index number** of the standard star (within the MIKE
;     structure), rather than the **object ID number**, which is what
;     MIKE_PROC and MIKE_FNTOBJ expect.  This isn't catastrophic, but
;     a bit inconsistent.
;
; EXAMPLES:
;   mike_skysub, mike, 1L, 1L, 1L, /FCHK
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_fitgap
;  mike_ordermask
;  bspline_iterfit
;  mike_getfil
;  mike_skyshift
;  mike_skysub_write
;
; REVISION HISTORY:
;   16-May-2003 Written by JXP
;   ??-2003     Modified extensively by SMB
;   26-Jun-2004 Update arc shifting
;-
;------------------------------------------------------------------------------

pro mike_skysub_write, mm, bset, skyfil, flg_ordr
  ;; Write Sky info (bset or spectrum)
  if flg_ordr EQ 0 then begin  ;; ENTIRE SPECTRUM
      if mm EQ 0 then mwrfits, bset, skyfil, /create, /silent $
      else mwrfits, bset, skyfil, /silent 
  endif else begin  ;; SELECT ORDERS
      a = findfile(skyfil+'*',count=na)
      if na EQ 0 then stop
      spawn, 'cp '+skyfil+' tmp_sky.fits'
      case mm of 
          0: begin
              mwrfits, bset, skyfil, /create,/silent 
              for i=1L,nordr-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
          end
          else: begin
              tmp_sky = xmrdfits('tmp_sky.fits', 1, /silent)
              mwrfits, tmp_sky, skyfil, /silent, /create
              ;; MIDDLE
              for i=1L,mm-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
              ;; NEW
              mwrfits, bset, skyfil, /silent 
              ;; REST
              for i=mm+1L,nordr-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
          endelse
      endcase
      spawn, '\rm -f tmp_sky.fits'
  endelse

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_skysub, mike, setup, obj_id, side, exp, CHK=chk, STD=std, $
                 ORDR=ordr, AIMG=aimg, USEOLD=useold, FCHK=fchk, $
                 DEBUG=debug, DOALL=doall, SCATTERED_LIGHT=iscattered_light

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_skysub, mike, setup, obj_id, side, [exp], /CHK, /STD, '
      print, '      ORDR=, AIMG=, /USEOLD, /FCHK, /DOALL [v1.2]'
      return
  endif 
  
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

;  Optional Keywords

;  Setup and setup
  case side of 
      1: ostep = 1
      2: ostep = -1
      else: stop
  endcase

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(mike.flg_anly NE 0 AND mike.side EQ side[0] AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                    strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_skysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse
      
;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

  skylinefile = getenv('MIKE_DIR')+'/pro/Skysub/UVES_SKY/uves.linelist'
  if (findfile(skylinefile))[0] NE '' then begin
     print, 'MIKE_SKYSUB: Reading ', +strcompress(skylinefile)
     readcol, skylinefile[0], uvesindx, uveswave, uvesint, uvesfwhm, uvesflux,$
           format='(i,d,f,f,f)'
     uveslog = alog10(uveswave)
     nuves= n_elements(uveslog)
  endif

;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, side=side)
  sv_ostr = ordr_str

  if (where(tag_names(ordr_str) EQ 'PROFILE0'))[0] NE -1 then begin
     if total(ordr_str.profile0) EQ 0 then $
       message, 'Need to run mike_slitflat to obtain slit profile'
  endif
  nordr = n_elements(ordr_str)

; Orders
  if not keyword_set( ORDR ) then begin  ;; Orders to sky subtract
      ordr = [ordr_str[0].order,ordr_str[nordr-1].order]
      flg_ordr = 0
  endif else begin
      if keyword_set( USEOLD ) then flg_ordr = 1 else flg_ordr = 0
      if n_elements(ordr) NE 2 then begin
          print, 'mike_skysub: ORDR must be a 2 element array'
          return
      endif
  endelse

  clr = getcolor(/load)


;  Loop

  for q=0L,n_elements(exp)-1 do begin
    
      thisone = indx[exp[q]]
      ;;  Read Arc
      if not keyword_set(AIMG) then arc_img = mike[thisone].arc_img $
      else arc_img = aimg
      ;; READ
      if x_chkfil(arc_img+'*') EQ 0 then begin
          print, 'mike_skysub: Arc file doesnt exist!', arc_img
          stop
          return
      endif
      print, 'mike_skysub: Reading Wave image file ', arc_img
      img_arc = xmrdfits(arc_img, /silent) 
      sz_arc = size(img_arc, /dimensions)

      ;; xoff
      xyoff = mike[thisone].arc_xyoff

      ;; Offset ordr_str
      ordr_str = sv_ostr
      shft = mike_shifti(xyoff, OSTR=ordr_str)

      ;; Open Obj file
;      objfil = mike[thisone].obj_fil
      objfil = mike_getfil('obj_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
      if keyword_set( DOALL ) then begin
          a = where((objstr.flg_anly MOD 2) NE 1, na)
          if na NE 0 then objstr[a].flg_anly = objstr[a].flg_anly + 1
      endif

      ;; IMG+VAR Fil 
;      imgfil = strtrim(mike[thisone].img_final,2)
      imgfil = mike_getfil('fin_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      skyfil = mike_getfil('sky_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      scattfil = mike_getfil('scatt_fil', setup, $
                             SUBFIL=mike[thisone].img_root, /name)
      print, 'mike_skysub: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Mask
      msk = lonarr(sz_img[0],sz_img[1])
      model = fltarr(sz_img[0],sz_img[1])
      maskimage = mike_ordermask(sz_img[0], sz_img[1], ordr_str, trim=0.)

      ;; Inter-order light fit
      nxbkpt = side EQ 1 ? 4 : 12
      nybkpt = side EQ 1 ? 4 : 8

      qafil = mike_getfil('qa_skysub', setup, SUBFIL=mike[thisone].img_root) ; added by jm08apr13nyu
      x_echskysub, img, ivar, img_arc, objstr, ordr_str, img_new, $
        CHK=chk, $
        ORDR=ordr, USEOLD=useold, FCHK=fchk, $
        DEBUG=debug, DOALL=doall, SCATTERED_LIGHT=scattered_light, $
        NXBKPT=nxbkpt, NYBKPT=nybkpt, OSTEP=ostep, SKYFIL=skyfil, $
        COLBIN=mike[thisone].colbin, QAFIL=qafil, MEDSCATT=medscatt, $
        SCATTTRIM=scatttrim, OBJTRIM=objtrim

      ;; Ouptut Obj Struct (with sky arrays)
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

      ;; CHK
      print, 'mike_skysub: Writing output to: ', imgfil
      if keyword_set( CHK ) or keyword_set( FCHK ) then $
        xatv, img_new, WVIMG=10^img_arc, /block, min=-50., max=50.

      ;; Ouptut New Image
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, ivar, imgfil, /silent
      mwrfits, img_new, imgfil, /silent

      ;; Scattered light image
      print, 'mike_skysub: Writing scattered image: ', scattfil
      mwrfits, scattered_light, scattfil, /create, /silent

      ;; COMPRESS
      print, 'mike_skysub: Compressing...'
      spawn, 'gzip -f '+imgfil
      spawn, 'gzip -f '+scattfil
  endfor
  
;  DONE
  print, 'mike_skysub: All done! '

  return
end

