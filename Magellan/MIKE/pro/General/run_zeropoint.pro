
; get overscan subtracted...
; mike_strct, mike
; mike.setup =1 
; mike.flg_anly=1 
; mike_setup, mike
; mike_subbias, mike, lindgen(n_elements(mike)), /arc

; newgain = 0.42
; redgain = 1.06

; read in blue
; bluefile = 'OV/ov_b0044.fits'
; arcfile  = 'OV/ov_b0047.fits'
; ostep = 1
  

pro run_zeropoint, bluefile, arcfile, gain, side, ostep, fwhm=fwhm, $
           nfind=nfind, reverse=reverse, ystart=ystart, bkspace=bkspace

  if N_PARAMS() LT 5 then begin
     print, 'run_zeropoint, bluefile, arcfile, gain, side, ostep,i' + $
               ' fwhm=fwhm, nfind = nfind, reverse=reverse'
     return
  endif

  if NOT keyword_set(fwhm) then fwhm = 3.3*2.3
  if NOT keyword_set(nfind) then nfind = side EQ 1 ? 35 : 35

  img = mrdfits(bluefile, 0, hdr)
  ncol   =(size(img))[1]
  nrow   =(size(img))[2]

  if nrow EQ 4096 then begin
     rebin = 1
     img = total(reform(img, ncol, 2, 2048), 2)
     nrow = 2048L
  endif

  if NOT keyword_set(ystart) then ystart = long(0.4*nrow)
  xstart = find_npeaks(djs_median(img[*,ystart-10:ystart+10],2), $
        nfind= nfind, minsep=2*fwhm, width=1.5*fwhm)
  print, min(xstart), max(xstart)
  xcen = trace_crude(8*img, yset=ycen, radius=fwhm, xerr=xerr, $
                     ystart=ystart,xstart=xstart[sort(xstart)], maxerr=1.0)
  norders=(size(xcen))[2]

  xy2traceset, ycen, xcen, tset, ncoeff=4, yfit=xfit1, invvar=(xerr LT 900), $
      maxdev=1.0, outmask=outmask
  xfit = mike_basis(tset, xcen, outmask, ncoeff=4, eigenvec=eigenvec, $
                         msktrc=msktrc, outstr=pca_out)

  boxcen = extract_boxcar(img, xfit, ycen, radius=fwhm)
  boxl   = extract_boxcar(img, xfit-1.5*fwhm, ycen, radius=fwhm/2.)
  boxr   = extract_boxcar(img, xfit+1.5*fwhm, ycen, radius=fwhm/2.)

stop

  arc = mrdfits(arcfile)
  if keyword_set(rebin) then arc = total(reform(arc, ncol, 2, 2048),2)

  guessarc = mike_getfil('guess_arc', SIDE=side, SZ=[1024,2048], /name)
  restore, guessarc
  boxarc = sv_aspec*0.0
  boxarc[*,0:norders-1] = extract_boxcar(arc, xfit, ycen, radius=fwhm)

  if keyword_set(reverse) then boxarc=reverse(boxarc)

  guess_fft = fft([sv_aspec[*], sv_aspec[*]*0])
  sv_fft    = fft([boxarc[*], boxarc[*]*0])
  cc        = float(fft( sv_fft * conj(guess_fft) ))
  ntot_fft = n_elements(cc)
  cc_peak = max(cc,total_shift)

  if total_shift GT ntot_fft/2 then total_shift = total_shift - ntot_fft
  djs_iterstat, cc, sigma=cc_sig 
  print, 'Max peak in full correlation ',cc_peak,' compared to RMS of ',cc_sig

  ordr_shift = round(double(total_shift) / nrow )
  row_shift = total_shift - ordr_shift* nrow 
 
  ;;;wavelengths....
  boxwave = boxcen*0.0
  for i=0,norders-1 do $
    if i+ordr_shift  GE 0 AND i+ordr_shift LT n_elements(all_arcfit) then $
      if keyword_set(all_arcfit[i+ordr_shift].ffit) then $
         boxwave[*,i] = x_calcfit(findgen(nrow)+row_shift, $
                            fitstr=all_arcfit[i+ordr_shift])

  x = findgen(norders) # replicate(1,nrow)
  xy2traceset, x, transpose(boxwave), wset, ncoeff=3, $
             invvar=transpose(boxwave NE 0), yfit=tbox

  final_wave = 10^(transpose(tbox))
  if keyword_set(reverse) then final_wave = reverse(final_wave)
  final_counts = final_wave*0.0
  ascii_out = strmid(bluefile, 0, strpos(bluefile,'.fit'))+'.asc'
  openw, ilun, ascii_out, /get_lun

for i=0, norders-1 do begin
;; 
;;   Multiply gain right here!
;;
  tempcts = 1.0d* gain * (boxcen[*,i] - (boxl[*,i]+boxr[*,i]))
  tempinv = 1.0d/(abs(tempcts)+100)
  tempinv[0:5] = 0.
  tempinv[nrow-5:nrow-1] = 0.
 
  if NOT keyword_set(bkspace) then bkspace = side EQ 1 ? 1.0 : 2.0

  fullbkpt = bspline_bkpts(final_wave[[5,nrow-5],i],bkspace=bkspace,no=4)
  bset = bspline_iterfit(final_wave[*,i], tempcts,fullbkpt=fullbkpt, $
             yfit=tempfit, invvar=tempinv, upper=20, lower=20)
  final_counts[*,i] = tempfit

  printf, ilun, string(final_wave[*,i],format='(f10.4)')+$
                string(tempfit,format='(f12.3)')+$
           string(ostep*i+ordr_shift+guess_ordr[0],format='(i4)'), format='(a)'
             
endfor
plot, final_wave, final_counts,ps=3
free_lun, ilun

return
end
