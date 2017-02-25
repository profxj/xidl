pro mage_quicklook, fitsfile, psout, SNR=snr

    ARCHIVE = getenv('MAGE_DIR')+'/ARCHIVE/'
  
  slitfile = strtrim(ARCHIVE,2)+'mage_archive_orders.fits' 
  sensfile = strtrim(ARCHIVE,2)+'GD108_mar09_std.sav' 
  wavefile = strtrim(ARCHIVE,2)+'arcimg.fits.gz' 
  pixfile  = strtrim(ARCHIVE,2)+'Piximg.fits' 
  ostruct  = strtrim(ARCHIVE,2)+'OStr_mage.fits'
  filstd   = strtrim(ARCHIVE,2)+'ObjStr_GD108.fits'
  ;  Generate a mask of the order numbers
  tset_slits = mrdfits(slitfile, 1)
  slitmask = long_slits2mask(tset_slits)
  ordermask = 0*slitmask
  ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L

  ; Read in the files
  exptime=0
  filenames = [fitsfile]
  nfiles = n_elements(filenames)
  IF nfiles GT 0 THEN BEGIN
     sciimg = 0
     var_tot = 0
     FOR ifile = 0L, nfiles-1L DO BEGIN
        mage_proc, filenames[ifile], sciimg1, scivar1, hdr=header ;, pixflatfile=pixflatfile
        exptime+=fxpar(header,'EXPTIME')
        sciimg = sciimg + sciimg1
        var_tot = var_tot + 1.0D/(scivar1 + (scivar1 EQ 0))
        var_tot = var_tot*(scivar1 GT 0)
     ENDFOR
     scivar = (var_tot GT 0)/(var_tot + (var_tot EQ 0))
  ENDIF ELSE  begin
     mage_proc, filenames, sciimg1, scivar, hdr=header ; pixflatfile=pixflatfile
     exptime+=fxpar(header,'EXPTIME')
  endelse
  
  ximg = long_slits2x(tset_slits, edgmask = edgmask,TOL_EDG=3)
  waveimg = xmrdfits(wavefile, 0, /fscale)
  piximg = xmrdfits(pixfile, 0)

  ; Mask out the object for first pass sky subtraction
  bsp = 0.6
  FWHM = 3.0
  print, "Masking out the object for first pass sky subtraction..."
  objstruct1 = long_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
                            , skymask = skymask1, objmask = objmask1 $
                            , nperslit = 1L, peakthresh = reduxthresh $
                            , fwhm = FWHM, ISLIT = ISLIT, /silent)
  print, "Performing 2D sky subtraction"
  skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask1, edgmask $
                         , bsp = bsp, ISLIT = ISLIT)

  
  box_rad = 5L
  velpix  = 22.0d

;; The order struct is output by running mage_traceorders.pro
  ordr_strct = xmrdfits(ostruct, 1)
  obj_strct  = m_mkobjstr(15)
  for i=0,14 do obj_strct[i].exp=exptime
  
  sciivar = 0.0*sciivar
  sciivar[where(skyimage GT 0)] = 1.0d/(skyimage[where(skyimage GT 0)])
  

  ; Object Finding
  print, "Fancy object finding..."
  ;m_fntobj, obj_strct, ordr_strct, sciimg-skyimage, sciivar, qafil
  obj_strct=mage_findobj(sciimg-skyimage,sciivar,waveimg,tset_slits $
                         ,filstd=filstd,chk=chk)
  ; Sanity check
  ;xatv, (sciimg-skyimage)*(slitmask GT 0), min = -20.0, max = 200.0 $
   ;     , wv = waveimg, /block
;  for qq=0L,14 do xatvplot, obj_strct[qq].trace[0:2048-1], findgen(2048)

  ; Quick boxcar extract
  
  print, "Boxcar extracting object..."
  ;m_extechopt, sciimg, sciimg-skyimage, sciivar, ordr_strct, obj_strct, $
  ;             velpix, img_arc=waveimg, skyfil=skyimage, helio=0.0, $
  ;             obj_name="Quicklook", ordermask=ordermask, /boxonly
  outfil='ql_obj.fits'
 mage_echextobj,sciimg,sciivar,header,skyimage $
                       ,piximg,waveimg,tset_slits $
                       ,obj_strct,outfil=outfil $
                       ,/boxcar, box_rad=box_rad
 obj_strct.exp=float(sxpar(header,'EXPTIME'))
 obj_strct.wave = obj_strct.box_wv
 obj_strct.fx = obj_strct.box_fx
 obj_strct.var = obj_strct.box_var
 mage_flux, sensfile, obj_strct, rej=0.05

 mage_combspec, obj_strct, fspec

 set_plot, "ps"
 device, filename=psout, /landscape, /color
 colors=getcolor(/load)
 res = 299792.458/4100.*0.7
 print, "Combining to 1D spectrum..."
 mage_1dspec, fspec, "ql.fits", "ql_sig.fits", "ql_comb.fits", resvel= res ;, /rebin

  t1 = mrdfits('ql.fits',0,hdr)
  t2 = mrdfits('ql_sig.fits')
  wv=10^(sxpar(hdr,"CRVAL1")+dindgen(n_elements(t1))*sxpar(hdr, "CDELT1"))

  !p.multi=[0,1,2]  

  sf = sort(t1)
  ninety = t1[sf[0.97*n_elements(sf)]]

  plot, wv, t1, xrange=[3200,6500], yrange=[-0.1,1.7*ninety], /xsty, /ysty, xtitle="Wavelength (A)", ytitle="Relative Flux"
  oplot, wv, t2, color=colors.green
  oplot, wv, t1-t1, color=colors.red

  plot, wv, t1, xrange=[6500,10000], yrange=[-0.1,1.5*ninety], /xsty, /ysty, xtitle="Wavelength (A)", ytitle="Relative Flux"
  oplot, wv, t2, color=colors.green
  oplot, wv, t1-t1, color=colors.red

  !p.multi=[0,2,2]


  for i=0, 14 do begin
     if not keyword_set(SNR) then begin
        sf = sort(obj_strct[i].fx)
        ninety = obj_strct[i].fx[sf[0.95*n_elements(sf)]]
        medfx = median(obj_strct[i].fx)
        plot, obj_strct[i].wave, obj_strct[i].fx, yrange=[-0.1*medfx, 1.4*ninety], xtitle="Wavelength (A)", ytitle="Counts"
        oplot, obj_strct[i].wave, sqrt(obj_strct[i].var), color=colors.green
     endif else begin
        sf = sort(obj_strct[i].fx)
        ninety = obj_strct[i].fx[sf[0.95*n_elements(sf)]]
        medfx = median(obj_strct[i].fx)
        mednoise = median(sqrt(obj_strct[i].var))
        plot, obj_strct[i].wave, obj_strct[i].fx/sqrt(obj_strct[i].var), yrange=[-0.1*medfx/mednoise, 1.4*ninety/mednoise]
     endelse
  endfor
  device, /close
  set_plot, "x"
  !p.multi=0

  x_specplot, "ql.fits", "ql_sig.fits", /block

;  stop

end

