;+
; NAME:
; mike_1dspec
;    Version 1.1
;
; PURPOSE:
;   Combines orders of a echfspec structure into a 1d spectrum.  
;   Outputs flux and error arrays into separate files (_F.fits, _E.fits)
;
; CALLING SEQUENCE:
;  
;   mike_1dspec, mike, setup, obj_id, side
;
; INPUTS:
;   mike      - MIKE structure
;   obj_id   -  Object ID  (e.g. 0L, 1L, etc)
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   1d flux      -   (fits file; FSpec/name_ech_F.fits)
;   1d error     -   (fits file; FSpec/name_ech_E.fits)
;
; OPTIONAL KEYWORDS:
;    /SILENT   - No text output
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file
;    ENDTRM=   - Trim order edges by this when calculating ratios
;    MINPIX1=  - Minimum 'good' pixels to calculate fitting profile
;    ORDNM=    - 0th or 1st order fitting (1 = 1st, 2 = 0th (default))
;    SNRMIN=   - Minimum S/N per pixel for consideration in fitting
;    MINPIX2=  - Minimum 'good' pixels to calculate simple 0th order fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   mike_1dspec, mike, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Jun-2004 Written by GEP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------


pro mike_1dspec, mike, setup, obj_id, side, exp_id, SILENT=silent, LIST=list, $
                 OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, CHK=chk, ORDNM=ordnm, $
                 KEEPALL=keepall, FLUX=flux, LROWBIN=lrowbin, NOFLUX=noflux,$
                 LSIDE=lside, MIDORDR=midordr

  if N_params() LT 4 and not keyword_set(LIST) then begin 
    print, 'Syntax - ' + $
      'mike_1dspec, mike, setup, obj_id, side, [exp_id], /SILENT ' + $
      'OBJ_NM=, ORDNM=, /CHK '
    print, '     LIST= [v1.0]'
    return
  endif

  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(MIDORDR) then midordr = 1L
  if keyword_set (LIST) then begin
      if keyword_set (LSIDE) then side = lside else side = 1L
      if not keyword_set (LROWBIN) then lrowbin = 2L
  endif

;  if side EQ 1 then begin
;      crosswave = [4983., 4914., 4847.2, 4782., $
;                   4719.,4657.,4597.,4538.,4482.,4426.,4371.,4318.,4268., $
;                   4217.,4168.,4120.,4072.,4026.,3981.,3938.,3894.,3852., $
;                   3811.,3771.,3731.,3693.,3654.,3618.,3581.,3545.,3510., $
;                   3484.,3450.,3415.,3386.,3353.,3320.,3292.]
;  endif else begin
;    crosswave = [9165., 8925, 8693.,8478.,8275.,8081.,7896.,7718., $
;                 7550.,7388.,7233.,7084., $
;                 6942.,6804.,6672.,6546.,6422.,6305.,6192.,6083.,5976.,5874., $
;                 5776.,5680.,5587.,5498.,5411.,5328.,5246.,5167.,5088.,5016., $
;                 4944.,4873.,4803.]
;  endelse

  if keyword_set (LIST) then $ 
    print, 'mike_1dspec:  2d --> 1d-ing the 1st file in the following list -- ', list
  
; Set exp
  if not keyword_set( LIST) then begin
    ;if not keyword_set( STD ) then begin
      allexp = where((mike.type EQ 'OBJ' OR mike.type EQ 'STD') $
                     AND mike.flg_anly NE 0 AND $
                     mike.setup EQ setup AND mike.obj_id EQ obj_id $
                     AND mike.side EQ side)
      if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
      nexp = n_elements(exp)
    ;endif else begin ;; STD
      ;exp = obj_id[0]
      ;nexp = 1
      ;obj_id = 99L
    ;endelse 
  endif else begin 
    obj_id = 99L
    readcol, list, files, FORMAT='A'
    nexp = n_elements(files)
    if nexp NE 3 then begin
        print, 'mike_1dspec LIST format requires 3 files to be listed in ' + $
                'list file:  first file is object file to be 1d-ed, second file ' +$
                'is flux output file, third is sigma output file.  Returning.'
        return
    endif
  endelse

;;;;;;;;;;;;;
;; Open file

  if not keyword_set( SILENT ) then $
    print, 'mike_1dspec: Loading up the files....'
  
  if not keyword_set( LIST ) then begin
      if keyword_set ( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(mike[exp[0]].Obj,2),/remove_all) + obj_nm
      if side EQ 1 then clrc = '_b' else clrc = '_r'
      infil = 'FSpec/'+subfil+clrc+'.fits'
      print, 'mike_1dspec: Reading ', infil
      spec = xmrdfits(infil,1)
;    head = xheadfits(infil)
      imgfil = mike_getfil('fin_fil', setup, SUBFIL=mike[allexp[0]].img_root, $
                           /name)  
  endif else begin
      spec = xmrdfits(files[0],1)
;    head = xheadfits(files[0])
      imgfil=files[0]
  endelse
  head = xheadfits(imgfil)
  mike_taghead, head

;;;;;;;;;;;;
;;; begin 1d

  use = where(spec.phys_ordr NE 0,nuse)
; gd = where(spec.wave NE 0.)
  gd = where(spec.wave[*,use] NE 0.,ngd) ; jm08apr14nyu - bug fix!
  if (ngd eq 0L) then message, 'Problem here!'

  minwv = min(spec.wave[gd])
  maxwv = max(spec.wave[gd])
  
  if not keyword_set(ORDNM) then ordnm = 2
  
  if not keyword_set( LIST ) then $
    velpix = (side EQ 1 ? 1.50d : 2.10d) * double(mike[exp[0]].rowbin) $
  else velpix = (side EQ 1 ? 1.50d : 2.10d) * double(lrowbin)
  
  cdelt = alog10(1.0d + velpix/299792.458d)
  npix = 150000L
  tot_wave = 10^(alog10(3000.0d) + dindgen(npix)*cdelt)
  weight = fltarr(npix)
  weight[*] = 0.
  tot_flux = fltarr(npix)
  tot_flux[*] = 0.
  sig = replicate(-1.,npix)

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
       for i=0L,nuse-1L do begin
          
          ordr = spec.phys_ordr[use[i]]
          hup_wv = (side EQ 1 ? 356100. : 343400.) / (ordr - 0.5)
          hlo_wv = (side EQ 1 ? 356100. : 343400.) / (ordr + 0.5)
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
            x_splot, spec.wave[u1,use[ii]], spec.fx[u1,use[ii]], ytwo=spec.fx[u2,use[ii-1]], xtwo=spec.wave[u2,use[ii-1]], ythr=spec.fx[u3,use[ii+1]], xthr=spec.wave[u3,use[ii+1]], psym3=-3, /block, TITLE=title
        endfor
    endif


   sw = spec.wave[*,use]
   sv = spec.var[*,use]
   sf = spec.fx[*,use]
   
   gd = where(sw GT 0 AND sv GT 0 AND sf NE 0.0000, ngd)
   sg = gd[sort(sw[gd])]
   min_ii = (where(abs(tot_wave-min(sw[gd])) LT 0.006))[0]
   max_ii = (where(abs(tot_wave-max(sw[gd])) LT 0.006))[0]
   if max_ii EQ -1 then mx = min(abs(tot_wave-max(sw[gd])),max_ii)

   ;; Combine
   ll = 0L
   for ii=min_ii,max_ii do begin
      uu = (ll + 6L) <  (ngd-1L)
;      a = where(abs(tot_wave[ii]-sw[sg[ll:uu]]) LT 0.006, na)
      a = where(abs(tot_wave[ii]-sw[sg[ll:uu]]) LT (0.006*double(mike[exp[0]].rowbin)), na)
      
;      if tot_wave[ii] GT 3466. then stop
      if na EQ 0 then begin
         print, 'No pixels at ', tot_wave[ii]
         continue
      endif
      
;      if na GT 2 then stop
      if na EQ 1 then begin
         tot_flux[ii] = sf[sg[a[0]+ll]] 
         sig[ii] = sqrt(sv[sg[a[0]+ll]])
      endif else begin
         h = sg[a + ll]
         wtmp = edge_weight[h] / sv[h]
         wtot = total(wtmp) 
         if wtot GT 0 then begin
            
;  Flux = Sum (W * fx) / Sum(W)
;  Var = Sum(W^2 * var) / Sum(W)^2
            
            tot_flux[ii] = total(wtmp * sf[h]) / wtot
            sig[ii] =  sqrt(total(wtmp^2 * sv[h]))/wtot
            outlier = where(sf[h] GT tot_flux[ii] + 5.0*sqrt(sv[h]), $
                            complement=keep)
            if outlier[0] NE -1 then begin
               if n_elements(keep) EQ 1 then begin
                  tot_flux[ii] = sf[h[keep[0]]]
                  sig[ii] = sqrt(sv[h[keep[0]]])
               endif else begin
                  h = h[keep]
                  wtmp = edge_weight[h] / sv[h]
                  wtot = total(wtmp) 
                  tot_flux[ii] = total(wtmp * sf[h]) / wtot
                  sig[ii] =  sqrt(total(wtmp^2 * sv[h]))/wtot
               endelse
            endif
         endif
;         print, tot_wave[ii], tot_flux[ii]
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
     else subfil = strcompress(strtrim(mike[exp[0]].Obj,2),/remove_all) + obj_nm
      if side EQ 1 then clrc = '_b' else clrc = '_r'
      outfil = 'FSpec/'+subfil+clrc+'_F.fits'
      sigfil = 'FSpec/'+subfil+clrc+'_E.fits'
    endif else begin
      outfil = files[1]
      sigfil = files[2]
    endelse

    if not keyword_set( SILENT ) then $
      print, 'mike_1dspec: Writing results to files: ', outfil, ',', sigfil

    mwrfits, tot_flux[a], outfil, head, /create, /silent
    mwrfits, sig[a], sigfil, head, /create, /silent

    return
end
