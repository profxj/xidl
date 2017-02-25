function tspec_wpix2image, tset2d, tset_slits, wave_struct = wave_struct $
                           , xshift = xshift, waveimg = waveimg

if N_PARAMS() LT 2 then begin
    print, 'Syntax: wpix = wpix2image(tset2d, tset_slits)'
    return, 0
endif

if n_elements(tset_slits) NE 2 then begin
    splog, 'tset_slits must be a 2 element structure of the slit edges'
    return, 0
endif


traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
if (keyword_set(xshift)) then begin
    left_edge = left_edge + xshift
    right_edge = right_edge + xshift
endif

ordermask = tspec_ordermask(tset_slits, order_vec=order_vec)

;   maskim  = lris_slits2mask(tset_slits, xshift = xshift)
slitpos = long_slits2x(tset_slits, xshift = xshift)

edge_sep = right_edge - left_edge

med_width = djs_median(edge_sep, 1)
slit_order = reverse(sort(med_width))
norders = n_elements(med_width)

if norders NE n_elements(tset2d) then begin
    splog, 'WARNING: Nslits in tset_slits does not equal the number is tset2d'
endif

dims = tset2d[0].dims
pix_image = dblarr(dims[0], dims[1])
IF KEYWORD_SET(WAVE_STRUCT) THEN BEGIN
    waveimg =  dblarr(dims[0], dims[1])
    ordr_str = replicate(create_struct('FLG_ANLY', 0L, 'ORDER', 0L) $
                         , norders)
    ordr_str.flg_anly = 1L
    ordr_str.ORDER = order_vec
ENDIF

for iorder = 0, norders-1 do begin
    coeff = tset2d[iorder].coeff2d
;    only take non-zero coefficients
    nxcoeff = max(where(total(abs(coeff), 1) GT 0)) + 1  
    nycoeff = max(where(total(abs(coeff), 2) GT 0)) + 1
    
    in = where(ordermask EQ order_vec[iorder], nin)
;    if nin EQ 0 then begin
;        splog, 'WARNING: This slit has no pixels mapped', iorder
;        continue
;    endif
    if nxcoeff EQ 0 OR nycoeff EQ 0 then begin
        splog, 'WARNING: coefficient numbers make no sense', $
               order_vec[iorder], nxcoeff, nycoeff 
        splog, '         cannot correct wavelength tilt'
        yimg = replicate(1.0, dims[0])  # findgen(dims[1])
        pix_image[in] = yimg[in]
    endif else begin
    
        slit_frac = slitpos[in]
        t = 2.0D*(double(in/dims[0]) - double(tset2d[iorder].xmin))/ $
          double(tset2d[iorder].xmax - tset2d[iorder].xmin) - 1.0D
        
        y = 2.0D*(double(slit_frac) - double(tset2d[iorder].ymin))/ $
          double(tset2d[iorder].ymax - tset2d[iorder].ymin) - 1.0D
        
        if tset2d[iorder].func EQ 'legendre' then begin
            tbasis = flegendre(t, nxcoeff)
            ybasis = flegendre(y, nycoeff)
        endif else begin
            splog, 'Not sure which basis function is being used'
            stop
            continue
        endelse
        
        for ix = 0, nxcoeff-1 do $
          for iy = 0, nycoeff-1 do $
          pix_image[in] = pix_image[in] + $
          tbasis[*, ix] * ybasis[*, iy] * coeff[ix, iy]        
    ENDELSE
    
    if keyword_set(wave_struct) then begin
;       create pixel input        
        npix = n_elements(pix_image[in])
        pix_nrm = 2. * (pix_image[in] - wave_struct.nrm[0])/wave_struct.nrm[1]
        worky = flegendre(pix_nrm[*], wave_struct.ny)
;       create order input        
        tsub = replicate(float(order_vec[iorder]), npix)
        t_nrm = 2. * (tsub - wave_struct.nrmt[0])/wave_struct.nrmt[1]
        ;; work2d and wv
        work2d = dblarr(npix, wave_struct.ny*wave_struct.no)
        workt = flegendre(t_nrm[*], wave_struct.no)

        for i = 0, wave_struct.no-1 do begin
            for j = 0, wave_struct.ny-1 do begin
                work2d[*, j*wave_struct.no+i] = worky[*, j] * workt[*, i]
            endfor
        endfor
        wv = dblarr(npix)
        wv[*] = work2d # wave_struct.res / ordr_str[iorder].order
        waveimg[in] = wv
     endif
     
 endfor


   return, pix_image
end
