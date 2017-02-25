;+ 
; NAME:
; mike_boxextrct   
;     Version 1.1
;
; PURPOSE:
;    Extract flux from 2D image to create ten 1D spectra (1 per order)
;    Output is written to the object structure (e.g. Extract/Obj_mike0024.fits)
;    The code only does boxcar extraction for now.
;
; CALLING SEQUENCE:
;   
;  mike_boxextrct, mike, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
;  RADIUS=
;
; INPUTS:
;   mike   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK    - Show final trace
;   /STD    - Extraction should be set for a standard star
;   /DEBUG  - Stop within extraction routine to check stuff
;   APER=   - Set aperture by hand (e.g. [5., 7.] )
;   RADIUS= - Size of window for setting aperture size (default: 20L)
;   ORDRS=  - Orders to extract (default: [0L,9L])
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  1)  The program begins extracting in order 9L (physical 6) and will
;  automatically calculate an aperture for that order.  If there is
;  insufficient flux in the following orders, it will adopt the value
;  from the next higher order. 
;
; EXAMPLES:
;   mike_boxextrct, mike, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_boxextrct, mike, setup, obj_id, side, exp, $
                    DEBUG=debug, CHK=chk, SVAPER=svaper, $
                    STD=std, APER=aper, RADIUS=radius, ORDRS=ordrs, $
                    SEDG_FIL=sedg_fil

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_boxextrct, mike, setup, obj_id, side, [exp], ' + $
        'RADIUS=, APER=, /DEBUG, /CHK'
      print, '          /STD, ORDRS= [v1.0]'
      return
  endif 
  
;  Optional Keywords
;  if not keyword_set( EXTREG ) then extreg = [25., 25.]
  if not keyword_set( RADIUS ) then radius = 20L
  if not keyword_set(ORDRS) then ordrs=[0L,9L]
  if not keyword_set( REJSIG ) then rejsig = 7.

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   mike.side EQ side AND $
                   mike.obj_id EQ obj_id AND strtrim(mike.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'mike_boxextrct: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STD star
      indx = obj_id[0]
      nindx = 1L
;      extreg = [45., 45.]
      radius = 40L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Setup and setup
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  case side of 
      1: begin
          nm = 'B'
          redblue = 0
      end
      2: begin
          nm = 'R'
          redblue = 1
      end
      else: stop
  endcase

;  Read in order structure
  ordr_fil = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
  if x_chkfil(ordr_fil+'*',/silent) EQ 0 then begin
      print, 'mike_boxextrct: Order structre doesnt exist. ' + $
        'Run mike_fndobj first!'
      return
  endif
  ordr_str = xmrdfits(ordr_fil,1,/silent)
  nordr = n_elements(ordr_str)
  if side EQ 1 then ordrs = [ordr_str[nordr-1].order,ordr_str[0].order] $
  else ordrs = [ordr_str[0].order,ordr_str[nordr-1].order]


;  Wavelength solution (put in a file)
  print, 'mike_boxextrct: Setting up wavlengths'
  mike_setwav, side, ordr_str, WV_ENDS=wv_mnx, BIN=mike[indx[0]].rowbin, $
    ALL_CRVAL=all_crval, CDELT=cdelt, NPIX=npix

; Wavelength regions for collapsing (finding profile)
  ;; Use inner 60%
  all_coll = wv_mnx
  all_coll[*,0] = wv_mnx[*,0] + (wv_mnx[*,1] - wv_mnx[*,0])*0.2
  all_coll[*,1] = wv_mnx[*,1] - (wv_mnx[*,1] - wv_mnx[*,0])*0.2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop

  for q=0L,n_elements(exp)-1 do begin

      ;; Offset edges
      xoff = mike[indx[exp[q]]].arc_xyoff[0]
      rnd_edg = [[round(ordr_str.lhedg+xoff)], [round(ordr_str.rhedg+xoff)]]
      sz = size(rnd_edg, /dimensions)
      rnd_edg = reform(rnd_edg,sz[0],nordr,2)

      print, 'mike_boxextrct: Reading files...'

      ;;;;;;;;;;;;;;
      ;; Open Obj file
      objfil = mike[indx[exp[q]]].obj_fil
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'mike_boxextrct: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
;      nobj = n_elements(objstr)

      ;;;;;;;;;;;;;;
      ;; SKY SUB Fil 
      imgfil = objstr[0].spec2d_fil
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'mike_boxextrct: No Image file!  Returning...'
          return
      endif
      print, 'mike_boxextrct: Image file -- ', imgfil
      head = xheadfits(imgfil)
      img = xmrdfits(imgfil, 2, /silent) ; SKY subtracted
      ivar = xmrdfits(imgfil, 1, /silent)
      ;; Need var for most of my routines
      var = 1./ivar
      sz_img = size(img, /dimensions)

      ;;  Read Arc
      arc_img = strtrim(mike[indx[exp[q]]].arc_img,2)
      print, 'mike_boxextrct: Arc -- ', arc_img
      if x_chkfil(arc_img+'*') EQ 0 then begin
          print, 'mike_boxextrct: No Arc file!  Returning...', arc_img
          return
      endif
      img_arc = xmrdfits(arc_img, /silent) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  HELIO correction
      ;; THIS NEEDS TO BE UPDATED FOR MIKE
      helio = x_keckhelio(mike[indx[exp[q]]].RA, mike[indx[exp[q]]].dec, $
                          mike[indx[exp[q]]].equinox, jd=mike[indx[exp[q]]].date)
      hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
      img_arc = img_arc + alog10(hel_corr)
      sxaddpar, head, 'HELIO', helio
;      xmodfits, imgfil, 0, head
      print, 'mike_boxextrct: Helio correction applied -- ', helio, hel_corr
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     LOOP ON ORDERS
      print, systime()
      for qq=ordrs[1], ordrs[0] do begin
          ;; Index
          mm = where(ordr_str.order EQ qq, nmm)
          if nmm EQ 0 then stop else mm = mm[0]
          ;;
          print, '-----------------------------------------------'
          print, 'mike_boxextrct: Extracting order ', string(qq,2)
           
          ;; Mask Sky
          msk = lonarr(sz_img[0],sz_img[1])
          lhs = (rnd_edg[*,mm,0]+2L) > 0L  
          rhs = (rnd_edg[*,mm,1]-2L) < (sz_img[0]-1)
          mni = min(round(lhs))
          mxi = max(round(rhs))
          a = where(lhs LE rhs, na)
          for ii=0L,na-1 do begin
              j = a[ii]
              msk[lhs[j]:rhs[j],j] = 1
          endfor
          
          ;; Transpose
          timg = transpose(img[mni:mxi,*]*msk[mni:mxi,*])
          tarc = 10^transpose(img_arc[mni:mxi,*]*msk[mni:mxi,*])
          tvar = transpose(var[mni:mxi,*]*msk[mni:mxi,*])
          
          ;; Run x_extobjbox
          x_extobjbox, timg, tarc, $
            [objstr[mm].xcen, objstr[mm].ycen-mni], fin_spec, $
            VAR=tvar, WVMNX=wv_mnx[mm,*], $
            APER=aper, DEBUG=debug, COLLMNX=all_coll[mm,*], $
            CRVAL1=all_crval[mm], CDELT=cdelt, NPIX=npix, $
            /REJ_CR, TOT_TRC=objstr[mm].trace[0:sz_img[1]-1]-mni, $
            REJSIG=rejsig, /REBINC, RADIUS=radius, BKAPER=svaper, $
            REDBLUE=redblue

          ;; Save the aperture
          svaper=aper

          if keyword_set( DEBUG ) then stop
          ;; CHK
          if keyword_set( CHK ) AND fin_spec.npix NE 0 then begin
              x_splot, fin_spec.wv, fin_spec.fx, YTWO=sqrt(fin_spec.var > 0.), $
                /block
          endif

          ;; Write to structure
          if fin_spec.npix NE 0 then begin
              objstr[mm].npix = fin_spec.npix
              objstr[mm].wave[0:fin_spec.npix-1] = fin_spec.wv
              objstr[mm].fx[0:fin_spec.npix-1] = fin_spec.fx
              objstr[mm].var[0:fin_spec.npix-1] = fin_spec.var
              ;; Aper
              objstr[mm].aper = fin_spec.aper
              ;; Flag
              objstr[mm].flg_anly = 1
          endif else objstr[mm].flg_anly = 0
      endfor

      print, systime()
      ;; Ouptut Spectra
      print, 'mike_boxextrct: Output spectrum in -- ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil
  endfor
  
;  DONE
  print, 'mike_boxextrct: All done! '
  return
end
