;+
; NAME:
; hires_1dspec
;    Version 1.1
;
; PURPOSE:
;   Combines orders of a echfspec structure into a 1d spectrum.  
;   Outputs flux and error arrays into separate files (_F.fits, _E.fits)
;
; CALLING SEQUENCE:
;   hires_1dspec, hires, setup, obj_id, chip
;
; INPUTS:
;   hires    - HIRES structure
;   setup
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   1d flux      -   (fits file; FSpec/name_ech_F.fits)
;   1d error     -   (fits file; FSpec/name_ech_E.fits)
;
; OPTIONAL KEYWORDS:
;  /KEEPALL -- Recommended and required currently.
;  LROWBIN -- Binning of CCD in spectral direction. Required for list
;             input.  The final spectrum has a pixels size of 1.3km/s
;             * ROWBIN.
;  /NOVAR --  Also write out NOVAR variance without object flux noise. 
;  /BZF   -- Blaze normalized and fluxed!  [rarely done]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hires_1dspec, hires, setup, obj_id, chip
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Jun-2004 Written by GEP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------


pro hires_1dspec, hires, setup, obj_id, chip, SILENT=silent, LIST=list, $
                 OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, CHK=chk, ORDNM=ordnm, $
                 KEEPALL=keepall, FLUX=flux, LROWBIN=lrowbin, $
                 LSIDE=lside, NOVAR = NOVAR, PIX_OFF=pix_off, $
                  BZF=bzf

  if N_params() LT 3 and not keyword_set(LIST) then begin 
    print, 'Syntax - ' + $
      'hires_1dspec, hires, setup, obj_id, [chip], /SILENT ' + $
      'OBJ_NM=, ORDNM=, /CHK, /BZF '
    print, '     LIST= [v1.0]'
    return
  endif

  if not keyword_set(PIX_OFF) then pix_off = 0L
  if not keyword_set(CHIP) then chip = [1,2,3]
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if keyword_set (LIST) then begin
      stop
      if not keyword_set (LROWBIN) then stop
  endif

  if keyword_set (LIST) then $ 
    print, 'hires_1dspec:  2d --> 1d-ing the 1st file in the following list -- ', list
  
; Set exp
  if not keyword_set( LIST) then begin
    ;if not keyword_set( STD ) then begin
      exp = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') $
                     AND hires.flg_anly NE 0 AND $
                     hires.setup EQ setup AND hires.obj_id EQ obj_id $
                     AND hires.chip EQ chip[0])

      nexp = n_elements(chip)

  endif else begin 
    obj_id = 99L
    readcol, list, files, FORMAT='A'
    nexp = n_elements(files)
    if nexp LT 3 then begin
        print, 'hires_1dspec LIST format requires at least 3 files to be listed in ' + $
                'list file:  first file is object file to be 1d-ed, second file ' +$
                'is flux output file, third is sigma output file.  Returning.'
        return
    endif
  endelse

;;;;;;;;;;;;;
;; Open file

  if not keyword_set( SILENT ) then $
    print, 'hires_1dspec: Loading up the files....'
  
  if not keyword_set( LIST ) then begin
      if keyword_set ( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(hires[exp[0]].Obj, 2)+obj_nm, /rem)
      for jj=0L,nexp-1 do begin
          case chip[jj] of 
              -1: clrc = '_S' 
              1: clrc = '_B' 
              2: clrc = '_G'
              3: clrc = '_R'
           endcase
          app = 'F'
          if keyword_set(BZF) then begin
             if bzf[jj] EQ 1 then app='bzf' 
          endif
          infil = 'FSpec/'+subfil+clrc+app+'.fits'
          print, 'hires_1dspec: Reading ', infil
          if jj EQ 0 then begin
             spec = xmrdfits(infil,1,/silent) 
             flg_novar = 0
          endif else begin
              tspec = xmrdfits(infil,1,/silent)
              gdo = where(tspec.phys_ordr NE 0,nuse)
              spec.phys_ordr[gdo] = gdo
              spec.wave[*,gdo] = tspec.wave[*,gdo] 
              spec.fx[*,gdo] = tspec.fx[*,gdo] 
              spec.var[*,gdo] = tspec.var[*,gdo] 
              if tag_exist(spec,'NOVAR') then begin
                  spec.novar[*, gdo] = tspec.novar[*, gdo]
                  flg_novar = 1
              endif else flg_novar = 0
              spec.npix[gdo] = tspec.npix[gdo] 
          endelse
      endfor
      imgfil = hires_getfil('fin_fil', setup, $
                            chip=chip[0], frame=hires[exp[0]].frame, $
                            /name)  
  endif else begin
      ;; List
      for jj=0L,nexp-3 do begin
          if jj EQ 0 then spec = xmrdfits(files[jj],1,/silent) $
          else begin
              tspec = xmrdfits(files[jj],1,/silent)
              gdo = where(tspec.phys_ordr NE 0,nuse)
              spec.phys_ordr[gdo] = gdo
              spec.wave[*,gdo] = tspec.wave[*,gdo] 
              spec.fx[*,gdo] = tspec.fx[*,gdo] 
              spec.var[*,gdo] = tspec.var[*,gdo]
              if tag_exist(spec,'NOVAR') then begin
                  spec.novar[*,gdo] = tspec.novar[*, gdo]
                  flg_novar = 1
              endif else flg_novar = 0
              spec.npix[gdo] = tspec.npix[gdo] 
          endelse
      endfor
      imgfil=files[0]
  endelse
  head = xheadfits(imgfil)
;  hires_taghead, head

;;;;;;;;;;;;
;;; begin 1d

  use = where(spec.phys_ordr NE 0,nuse)
  
  gd = where(spec.wave NE 0.)
  minwv = min(spec.wave[gd])
  maxwv = max(spec.wave[gd])
  
  if not keyword_set(ORDNM) then ordnm = 2
  
;  if not keyword_set( LIST ) then $
;    velpix = (side EQ 1 ? 1.50d : 2.10d) * double(hires[exp[0]].rowbin) $
;  else velpix = (side EQ 1 ? 1.50d : 2.10d) * double(lrowbin)
  if keyword_set(LIST) then rbin = lrowbin else rbin = hires[exp[0]].rowbin
  velpix = 1.30d  * double(rbin)
  
  cdelt = alog10(1.0d + velpix/299792.458d)
  npix = 320000L
  tot_wave = 10^(alog10(3000.0d) + (dindgen(npix)+PIX_OFF)*cdelt)
  weight = fltarr(npix)
  weight[*] = 0.
  tot_flux = fltarr(npix)
  sig = replicate(-1.,npix)
  nosig = replicate(-1., npix)

  sxaddpar, head, 'CDELT1', cdelt
  sxaddpar, head, 'CRPIX1', 1
  sxaddpar, head, 'CTYPE1', 'LINEAR'
  sxaddpar, head, 'DC-FLAG', 1
  sxaddpar, head, 'BITPIX', -32
  sxaddpar, head, 'NAXIS', 1
  sxaddpar, head, 'NAXIS1', n_elements(tot_flux)
  sxdelpar, head, 'NAXIS2'
  
    ;; Get S/N for each order overlap, zero out fx, sig in bad blaize areas

    edge_weight = 1.0*(spec.var[*,use] GT 0)

    if not keyword_set(KEEPALL) then begin
        stop  ;; Use /KEEPALL until futher notice
        for i=0L,nuse-1L do begin
            
            ordr = spec.phys_ordr[use[i]]
            hup_wv = (chip EQ 1 ? 356100. : 343400.) / (ordr - 0.5)
            hlo_wv = (chip EQ 1 ? 356100. : 343400.) / (ordr + 0.5)
            fsr = hup_wv - hlo_wv
            
            gd = where(spec.wave[*,use[i]] LE hlo_wv and $
                       spec.wave[*,use[i]] NE 0. ,ngd)
            
            if ngd GT 0 then begin
                fixfun = (1.0 - (hlo_wv  - spec.wave[gd,use[i]])/ fsr )  > 0
                edge_weight[gd,i] = fixfun
            endif
            
            gd = where(spec.wave[*,use[i]] GE hup_wv,ngd)
            if ngd GT 0 then begin
                fixfun = (1.0 - (spec.wave[gd,use[i]] - hup_wv)/ fsr )  > 0
                edge_weight[gd,i] = fixfun
            endif
            
        endfor
    endif
  

    if keyword_set( CHK ) then begin
        for ii=1,nuse-2 do begin
            title = strtrim(string(use[ii]),2)
            u1 = where(spec.wave[*,use[ii]] GT 0)
            u2 = where(spec.wave[*,use[ii-1]] GT 0)
            u3 = where(spec.wave[*,use[ii+1]] GT 0)
            x_splot, spec.wave[u1,use[ii]], spec.fx[u1,use[ii]], $
                     ytwo=spec.fx[u2,use[ii-1]], $
                     xtwo=spec.wave[u2,use[ii-1]], $
                     ythr=spec.fx[u3,use[ii+1]], $
                     xthr=spec.wave[u3,use[ii+1]], $
                     psym3=-3, /block, TITLE=title
        endfor
    endif
    

   sw = spec.wave[*,use]
   sv = spec.var[*,use]
   if flg_novar then snov = spec.novar[*,use]
   sf = spec.fx[*,use]
   
   gd = where(sw GT 0 AND sv GT 0 AND sf NE 0.0000, ngd)
   sg = gd[sort(sw[gd])]
   min_ii = (where(abs(tot_wave-min(sw[gd])) LT 0.006))[0]
   if min_ii LT 0 then begin
       print, 'hires_1dspec: It is likely your spectrum goes bluer than 3000Ang.'
       print, 'hires_1dspec: Try using a negative PIX_OFF value'
       print, 'hires_1dspec: Returning...'
       return
   endif
   max_ii = (where(abs(tot_wave-max(sw[gd])) LT 0.006))[0]


    ll = 0L
    for ii=min_ii,max_ii do begin
       uu = (ll + 6L) <  (ngd-1L)
       a = where(abs(tot_wave[ii]-sw[sg[ll:uu]]) LT 0.006, na)

       if na EQ 0 then continue

       if na EQ 1 then begin
         tot_flux[ii] = sf[sg[a[0]+ll]] 
         sig[ii] = sqrt(sv[sg[a[0]+ll]])
         if flg_novar then nosig[ii] = sqrt(snov[sg[a[0]+ll]])
       endif else begin
         h = sg[a + ll]
         wtmp = edge_weight[h] / sv[h]
         wtot = total(wtmp) 
         if wtot GT 0 then begin

;  Flux = Sum (W * fx) / Sum(W)
;  Var = Sum(W^2 * var) / Sum(W)^2

           tot_flux[ii] = total(wtmp * sf[h]) / wtot
           sig[ii] =  sqrt(total(wtmp^2 * sv[h]))/wtot
           if flg_novar then nosig[ii] = sqrt(total(wtmp^2 * snov[h]))/wtot
           outlier = where(sf[h] GT tot_flux[ii] + 5.0*sqrt(sv[h]), $
                    complement=keep)
           if outlier[0] NE -1 then begin
             if n_elements(keep) EQ 1 then begin
                tot_flux[ii] = sf[h[keep[0]]]
                sig[ii] = sqrt(sv[h[keep[0]]])
                if flg_novar then nosig[ii] = sqrt(snov[h[keep[0]]])
             endif else begin
                h = h[keep]
                wtmp = edge_weight[h] / sv[h]
                wtot = total(wtmp) 
                tot_flux[ii] = total(wtmp * sf[h]) / wtot
                sig[ii] =  sqrt(total(wtmp^2 * sv[h]))/wtot
                if flg_novar then nosig[ii] = sqrt(total(wtmp^2 * snov[h]))/wtot
             endelse
           endif
         endif
       endelse
       
       ll = max(a) + ll

      endfor

    minwv = tot_wave[min_ii]
    tot_pix = max_ii - min_ii + 1L
    a = lindgen(tot_pix) + min_ii
   
    if keyword_set( CHK ) then $
      x_splot, tot_wave[a], tot_flux[a], ytwo=sig[a], /block

    sxaddpar, head, 'CRVAL1', alog10(minwv)

    if not keyword_set( LIST ) then begin
     if keyword_set ( OUTNM ) then subfil = outnm+obj_nm $
     else subfil = strcompress(strtrim(hires[exp[0]].Obj, 2) + obj_nm, /rem)
     outfil = 'FSpec/'+subfil+'_f.fits'
     sigfil = 'FSpec/'+subfil+'_e.fits'
     if flg_novar then nosigfil = 'FSpec/'+subfil+'_n.fits'
    endif else begin
      outfil = files[nexp-2]
      sigfil = files[nexp-1]
    endelse

    if not keyword_set( SILENT ) then $
      print, 'hires_1dspec: Writing results to files: ', outfil, ',', sigfil

    mwrfits, tot_flux[a], outfil, head, /create, /silent
    mwrfits, sig[a], sigfil, head, /create, /silent
    IF (FLG_NOVAR) THEN mwrfits, nosig[a], nosigfil, head $
      , /create, /silent

    return
end
