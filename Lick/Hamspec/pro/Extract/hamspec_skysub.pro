;+ 
; NAME:
; hamspec_skysub   
;     Version 1.2
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_hamspec0020.fits).  The program also outputs
;    information related to the sky fit in the Sky directory.  The
;    user has control over which order to subtract (ORDR=), the
;    wavelength image to use (AIMG=).  
;
;  The main call is to x_echskysub. This program preforms the following steps:
;
;    1.  Offset the Order structure by the appropriate shift
;    2.  Subtract scattered light (hamspec_fitgap)
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
;  hamspec_skysub, hamspec, setup, obj_id, [exp], /CHK, /STD, ORDR=,
;  AIMG=, /USEOLD, /FCHK
;
; INPUTS:
;   hamspec   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
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
;  IOBJTRIM - Input trimming of the slit edge.  Default is 1.5 pixels
;             for G,R and 0 for B.  Note, you need to input a negative
;             value to trim and a positive to add pixels.
;
; OPTIONAL OUTPUTS:
;  SCATTERED_LIGHT= -- Scattered light image.  It may be useful to
;                      save and use this if you are subtracting individual
;                      orders
;  X_ECHSKYSUB KEYW ==  /SCTCHK
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_skysub, hamspec, 1L, 1L, 1L, /FCHK
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_fitgap
;  hamspec_ordermask
;  bspline_iterfit
;  hamspec_getfil
;  hamspec_skyshift
;  hamspec_skysub_write
;
; REVISION HISTORY:
;   16-May-2003 Written by JXP
;   ??-2003     Modified extensively by SMB
;   26-Jun-2004 Update arc shifting
;   04-Feb-2013   Modified from HIRES by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_skysub, hamspec, setup, obj_id, exp, CHK=chk, STD=std, $
                  ORDR=iordr, AIMG=aimg, FCHK=fchk, SCTCHK=sctchk, $
                  DEBUG=debug, DOALL=doall, IOBJTRIM=iobjtrim, $
                  SCATTERED_LIGHT=scattered_light,$
                  NOMEDSCATT=nomedscatt, _EXTRA=extra

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hamspec_skysub, hamspec, setup, obj_id, [exp], /CHK, /STD, '
      print, '      ORDR=, AIMG=, /FCHK, [v1.2]'
      return
  endif 
  if not keyword_set(SCATTTRIM) then scatttrim = 1.
  if keyword_set(CHK) then window, 0, title='hamspec_skysub'
  
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  ;; Scattered light parameters
  nxbkpt = 5
  nybkpt = 5
  MEDSCATT=0L


  ;;  Find all relevant obj
  if not keyword_set( STD ) then begin
     indx = where(hamspec.flg_anly NE 0 AND $
                  hamspec.obj_id EQ obj_id AND hamspec.setup EQ setup AND $
                  (strtrim(hamspec.type,2) EQ 'OBJ' OR $
                   strtrim(hamspec.type,2) EQ 'STD'), nindx)
     if nindx EQ 0 then begin
        print, 'hamspec_skysub: No images to sky subtract!', obj_id
        return
     endif
  endif else begin
     indx = obj_id[0]
     nindx = 1L
  endelse
      
  ;; Exposures
  if size(exp,/type) EQ 0 then exp = lindgen(nindx)
      
  ;; Read in order structure
  ordr_str = hamspec_getfil('ordr_str', setup)
  sv_ostr = ordr_str
      
  if (where(tag_names(ordr_str) EQ 'PROFILE0'))[0] NE -1 then begin
     if total(ordr_str.profile0) EQ 0 then $
        message, 'Need to run hamspec_slitflat to obtain slit profile'
  endif
  nordr = n_elements(ordr_str)
      
  ;;  Loop
  for q=0L,n_elements(exp)-1 do begin
     
     thisone = indx[exp[q]]
     ;;  Read Arc
     if not keyword_set(AIMG) then begin
        arc_img = hamspec_getfil('arc_mkaimg', setup, /name, CHKFIL=chkf) 
     endif else arc_img = aimg
     ;; READ
     if x_chkfil(arc_img+'*') EQ 0 OR strlen(arc_img) EQ 0 then begin
        print, 'hamspec_skysub: Arc file doesnt exist!', arc_img
        stop
        return
     endif
     print, 'hamspec_skysub: Reading Wave image file ', arc_img
     img_arc = xmrdfits(arc_img, /silent) 
     sz_arc = size(img_arc, /dimensions)
          
     ;; OBJTRIM
     if not keyword_set(iobjtrim) then objtrim = 0.  $
     else objtrim = iobjtrim
          
     ;; Offset ordr_str
     ordr_str = sv_ostr
     ordr_str.lhedg = ordr_str.lhedg +  hamspec[thisone].arc_xyoff[0]
     ordr_str.rhedg = ordr_str.rhedg +  hamspec[thisone].arc_xyoff[0]
     ordr_str.lhc[0] = ordr_str.lhc[0] +  hamspec[thisone].arc_xyoff[0]
     ordr_str.rhc[0] = ordr_str.rhc[0] +  hamspec[thisone].arc_xyoff[0]
;      shft = hamspec_shifti(xyoff, OSTR=ordr_str)
          
     ;; Shift arc?
     rnshft = round(hamspec[thisone].arc_xyoff[0])
     if rnshft NE 0 then begin
        ;; Shift
        img_arc = shift(img_arc, rnshft, 0)
        ;; Zero out edge
        if rnshft GT 0 then img_arc[0:rnshft-1,*] = 0. $
        else img_arc[sz_arc[0]-rnshft:*,*] = 0.
     endif

     ;; Open Obj file
;      objfil = hamspec[thisone].obj_fil
     objfil = hamspec_getfil('obj_fil', /name, $
                             FRAME=hamspec[indx[exp[q]]].frame)
     objstr = xmrdfits(objfil, 1, STRUCTYP='hamspecobjstrct', /silent)
     if keyword_set( DOALL ) then begin
        a = where((objstr.flg_anly MOD 2) NE 1, na)
        if na NE 0 then objstr[a].flg_anly = objstr[a].flg_anly + 1
     endif
     
     ;; IMG+VAR Fil 
;      imgfil = strtrim(hamspec[thisone].img_final,2)
     imgfil =hamspec_getfil('fin_fil', /name, CHKFIL=chkf, $
                            FRAME=hamspec[thisone].frame)
     skyfil =hamspec_getfil('sky_fil', /name, CHKFIL=chkf, $
                            FRAME=hamspec[thisone].frame)
     scattfil =hamspec_getfil('scatt_fil', /name, CHKFIL=chkf, $
                                 FRAME=hamspec[thisone].frame)
     print, 'hamspec_skysub: Reading Image files... ', imgfil
     img = xmrdfits(imgfil, 0, head, /silent)
     ivar = xmrdfits(imgfil, 1, /silent)
          
     ;; QAFIL
     qafil = hamspec_getfil('qa_skysub', setup, FRAME=hamspec[thisone].frame)
     
     if keyword_set(IORDR) then ordr=iordr
     x_echskysub, img, ivar, img_arc, objstr, ordr_str, img_new, $
                  CHK=chk, ORDR=ordr, USEOLD=useold, FCHK=fchk, $
                  DEBUG=debug, DOALL=doall, SCATTERED_LIGHT=scattered_light, $
                  NXBKPT=nxbkpt, NYBKPT=nybkpt, SKYFIL=skyfil, $
                  COLBIN=hamspec[thisone].colbin, QAFIL=qafil, $
                  MEDSCATT=medscatt, SCTCHK=sctchk, $
                  SCATTTRIM=scatttrim, OBJTRIM=objtrim, _EXTRA=extra
     delvarx, ordr
          
     ;; Ouptut Obj Struct (with sky arrays)
     mwrfits, objstr, objfil, /create
     spawn, 'gzip -f '+objfil
          
     ;; Ouptut New Image
     print, 'hamspec_skysub: Writing output to: ', imgfil
     if keyword_set( CHK ) or keyword_set( FCHK ) then $
        xatv, img_new, WVIMG=10^img_arc, /block, min=-50., max=50.
     
     mwrfits, img, imgfil, head, /create, /silent
     mwrfits, ivar, imgfil, /silent
     mwrfits, img_new, imgfil, /silent
     
     ;; Scattered light image
     print, 'hamspec_skysub: Writing scattered image: ', scattfil
     mwrfits, scattered_light, scattfil, /create, /silent
     
     ;; COMPRESS
     print, 'hamspec_skysub: Compressing...'
     spawn, 'gzip -f '+imgfil
     spawn, 'gzip -f '+scattfil
  endfor
  
;  DONE
  print, 'hamspec_skysub: All done! '

  return
end

