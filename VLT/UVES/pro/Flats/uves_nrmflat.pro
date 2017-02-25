;+ 
; NAME:
; uves_nrmflat
;     Version 1.1
;
; PURPOSE:
;  To normalize the flat (remove the blaze) to create an image useful
;  for pixel-to-pixel corrections.  Ideally, I recommend you use one
;  of the archived flat frames as the flats created here have the
;  following problems:
;  1.  They fail at the slit edges
;  2.  They contain scattered light which cannot be removed
;  3.  Large chip defects cause more trouble
;
; CALLING SEQUENCE:
;   
;  uves_nrmflat, uves, setup, [side], /CHK, 
;
; INPUTS:
;   uves     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  A normalized flat image
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_nrmflat, uves, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro uves_nrmflat, uves, setup, side, CHK=chk, CLOBBER=clobber, DEBUG=debug, $
                  SCATT=scatt_img, XSEP=xsep, AFRAME=aframe

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_nrmflat, uves, setup, [side], /CLOBBER,'
      print, '           /CHK, /DEBUG, XSEP= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( AFRAME ) then stop
  if not keyword_set( SIDE ) then side = [1L]
  if not keyword_set( NRMMED ) then NRMMED = 5L
  if not keyword_set( XSEP ) then xsep = -1L

;  x_psclose
;  !p.multi=[0,1,1]

; Loop on chip

  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      idx = where(uves.setup EQ setup and uves.flg_anly NE 0 and $
                  uves.side EQ qq and uves.type EQ 'TFLT',nidx)
      if nidx EQ 0 then stop
      rbin = uves[idx[0]].rowbin
      wcen = uves[idx[0]].xdangl
      twcen = strtrim(round(wcen),2)

      ;; Vignetting due to CCD cover
;      upedg = round(3990./rbin) 
      ;; Chip
      case qq of 
          1: begin
              print, 'uves_nrmflat: Normalizing BLUE Flat' 
              nxbkpt=5L 
              nybkpt=5L
          end
          2: begin
              print, 'uves_nrmflat: Normalizing GREEN Flat' 
              nxbkpt=5L 
              nybkpt=5L
          end
          3: begin
              print, 'uves_nrmflat: Normalizing RED Flat' 
              nxbkpt=5L 
              nybkpt=5L
          end
          else: stop
      endcase

      ;; Check output file
      out_fil = uves_getfil('nqtz_fil', setup, WCEN=wcen, /name, CHKFIL=chkf)
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'uves_nrmflat: Normalized flat exists, moving on..'
          continue
      endif

      ;; Read in flat 
      flat_fil = uves_getfil('qtz_fil', setup, WCEN=wcen,/name)
      flat = xmrdfits(flat_fil, 0, /silent)
      flat_ivar = xmrdfits(flat_fil, 1, /silent)

      sz = size(flat, /dimensions)
      nrm_flat = fltarr(sz[0],sz[1])
      nrm_ivar = fltarr(sz[0],sz[1])

      ;; Read in Arc
      img_arc = uves_getfil('arc_img', setup, FRAME=AFRAME, WCEN=wcen)

      ;; Read in Order structure
      ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen)

      ;; Fit parameters
      slit_cen = round((ordr_str.lhedg + ordr_str.rhedg)/2.)
      med_img = fltarr(sz[1], NRMMED*2+1L) 

      ;; Subtract scattered light
;      x_modelslit, tflat, tflativar, ordr_str, scat_model=gapfit
;      model_slit, tflat, tflativar, ordr_str, scat_model=gapfit
;      tflat_sub = tflat  - transpose(x_medianrow(transpose(gapfit), 31))
      maskimage = x_ordermask(sz[0], sz[1], ordr_str, trim=2)
      scatt_img = x_fitgap(flat, flat_ivar, maskimage, $
                           nxbkpt=nxbkpt, nybkpt=nybkpt)
      flat = flat - scatt_img
      if keyword_set(CHK) then xatv, flat, /blo

      ;; Loop on Orders
      nordr= n_elements(ordr_str)
      for mm=0L,nordr-1 do begin
          inorder = where(maskimage mod 10000L EQ ordr_str[mm].order $
                          AND img_arc GT 3. AND flat_ivar GT 0.)
          wave_sort = inorder[sort(img_arc[inorder])]
          xstart = 1.0d*(wave_sort mod sz[0])
          ystart = wave_sort /   sz[0]
          slit_length = ordr_str[mm].rhedg[ystart] - ordr_str[mm].lhedg[ystart]

          ;;Ordrcen is just center of flat-field order
          ordrcen  =  (ordr_str[mm].lhedg[ystart] + $
                       ordr_str[mm].rhedg[ystart])/2.0 ; Offset
          
          ywave = x_qckwav(xstart-ordrcen, ystart, ordr_str[mm].arc_m, $
                           arc_slope=arc_slope, slit_dist=slit_dist)
          
          slit_frac  = frac_order(ordr_str[mm], xstart, ywave)
          profile = x_slitprofile_return(slit_frac, ystart, ordr_str[mm])
          
          nrm_flat[wave_sort] = (flat)[wave_sort] / $
                     (profile + (profile EQ 0)) * (profile GT 0)
          nrm_ivar[wave_sort] = (flat_ivar)[wave_sort] / $
                     (profile + (profile EQ 0)) * (profile GT 0)
          
          ;; Fit the central flux
          scen = (ordr_str[mm].lhedg + ordr_str[mm].rhedg) / 2.
          gd = where(scen GT NRMMED and (scen+NRMMED) LT (sz[0]-1), ngd)
          frstj = gd[0]
          lstj = gd[ngd-1] ;< upedg
          for j=frstj,lstj do med_img[j,*] = nrm_flat[scen[j]-NRMMED: $
                                                      scen[j]+NRMMED,j]
          med_nrm = djs_median(med_img, 2)

          bset = bspline_iterfit(float(gd),med_nrm[frstj:lstj], $
                                 yfit=fit, everyn=round(15L/float(rbin)))

          ;; Divide it out!
          nrm = bspline_valu(ystart, bset)
          nrm_flat[wave_sort] = nrm_flat[wave_sort] / nrm
          nrm_ivar[wave_sort] = nrm_ivar[wave_sort] / nrm

      endfor

      ;; CHK
      if keyword_set(CHK) then xatv, nrm_flat, min=0.9, max=1.1, /block

      ;; Write 
      print, 'uves_nrmflat: Writing ', out_fil
      mwrfits, nrm_flat, out_fil, /create
      mwrfits, nrm_ivar, out_fil
      mwrfits, scatt_img, out_fil
      spawn, 'gzip -f '+out_fil

          
  endfor

  ;; 
  print, 'uves_nrmflat: All done!'

  return
end

