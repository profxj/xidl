;+ 
; NAME:
; z_arcpairs   
;     Version 1.1
;
; PURPOSE:
;  To identify the wavelength solution given arc lines and a list
;
; CALLING SEQUENCE:
;   
;  good = z_arcpairs(peaks, l, disp0, [dr1, NCOEFF=, INC=, TAN_BLAZE=])
;
; INPUTS:
;
;  peaks :	arc line positions in pixel space
;  l     :      possible arc line wavelengths
;  disp0 :      guess at central dispersion per pixel, needs to be within dr
;                of the correct dispersion
;
; RETURNS:
;   
;  good  :      An index list of the best meatched wavelengths based on the
;                number of unique hits
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
;  dr1   :      range of dispersion to search, 10% is default [0.9,1.1]
;  tan_blaze    For high dispersion gratings, the tangent of the blaze
;                angle, usually referred to a R. (like the R=2 grating in HIRES)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   May-2005 Written by SB
;------------------------------------------------------------------------------
function z_arcpairs, peaks, l, disp0, dr1 = dr1, ncoeff = ncoeff, inc = inc $
                    , tan_blaze = tan_blaze, CENTER_WAVE = CW_PL $
                    , DISP = RR_PL

   if N_PARAMS() LT 3 then begin
     print, 'Syntax:  lfit = arc_pairs(p, l, disp0)'
     return, 0
   endif

   if NOT keyword_set(ncoeff) then ncoeff=3L
   if NOT keyword_set(inc) then inc=80.0

   np =  n_elements(peaks)
   nl = n_elements(l)

   if NOT keyword_set(dr1) then dr = [0.8, 1.2]*disp0 $
   else dr = dr1
   
   dr = dr[sort(dr)]
   ds = dr
   

   central_pix = mean(peaks)
   dl = l # replicate(1,nl) - l ## replicate(1,nl)

   s = sort(peaks)
   p = peaks[s]
   for i=0, np-2 do begin
     if i EQ 0 then  begin
       dp = p[1:*] - p[0] 
       ilist = replicate(0,np-1-i)
       jlist = lindgen(np-i-1)+i+1
     endif else begin
       dp = [dp, p[i+1:*] - p[i]]
       ilist = [ilist, replicate(i,np-1-i)]
       jlist = [jlist, lindgen(np-i-1)+i+1]
     endelse
   endfor 

   sp = sort(dp)
   up = long(dp*0) -1L
   lo = long(dp*0) -1L
   sl = sort(dl)
   if disp0 LT 0 then sp = reverse(sp)

   i=0 
   ml = min(where(dl[sl] GE ds[0]*dp[sp[0]]))
   lo[0] = ml
   mh = max(where(dl[sl] LE ds[1]*dp[sp[0]]))
   up[0] = mh > 0  ;; Added by JXP 7/4/2007

   for i=1, n_elements(dp)-1 do begin
     ml = min(where(dl[sl[lo[i-1]:*]] GE ds[0]*dp[sp[i]]))
     lo[i] = ml+lo[i-1]
     mh = max(where(dl[sl[up[i-1]:*]] LE ds[1]*dp[sp[i]])) > 0 ;; Added by JXP
     up[i] = mh+up[i-1]
   endfor
   
   n = long(total(up - lo + 1))
   i1 = lonarr(n)
   i2 = i1
   hh = i1
   rr = fltarr(n)    
   cw = fltarr(n) 

   count = 0
   for i=0, n_elements(dp)-1 do begin
     nhere = up[i] - lo[i] + 1
     if nhere LE 0 then continue
     i1[count:count+nhere-1] = ilist[sp[i]]
     i2[count:count+nhere-1] = jlist[sp[i]]
     hh[count:count+nhere-1] = sl[lo[i]:up[i]]
     rr[count:count+nhere-1] = dl[hh[count:count+nhere-1]]/dp[sp[i]]
     cw[count:count+nhere-1] = rr[count:count+nhere-1] * $
           (central_pix-p[ilist[sp[i]]]) + l[hh[count:count+nhere-1] / nl]
     count = count + nhere
   endfor

  min_cw = min(cw)
  max_cw = max(cw)
  bin_cw = 2.0*abs(disp0)
  bin_rr = (dr[1] - dr[0])/200.
  his = hist_2d(cw, rr, min1=min_cw, max1=max_cw, bin1=bin_cw, $
                            min2=dr[0], max2=dr[1], bin2=bin_rr)

  hiscol = (size(his))[1]
  hisrow = (size(his))[2]
  smooth_his = convol(1.0*his,replicate(1,3,3))
  max_his = max(smooth_his, pl)
  cw_pl = min_cw + bin_cw * ((pl mod hiscol) + (hiscol mod 2)/2. )
  rr_pl = dr[0] + bin_rr * ((pl / hiscol) +  (hiscol mod 2)/2.)

  splog, 'Best fit center wave and dispersion are: ', cw_pl, rr_pl
  poss = where(cw GE cw_pl - inc*bin_cw AND cw LE cw_pl + inc*bin_cw $
           AND rr GE rr_pl - inc*bin_rr AND rr LE rr_pl + inc*bin_rr, nposs)

  pix_i = [i1[poss], i2[poss]]
  line_i = [hh[poss] / nl, hh[poss] mod nl]
  dp_i = p[i2[poss]] - p[i1[poss]]
  mp_i = (p[i2[poss]] + p[i1[poss]]) / 2.
  rr_i = dl[hh[poss]] / dp_i 
  cw_i = (central_pix - p[i1[poss]]) * rr_i + l[hh[poss] / nl]

  if keyword_set(tan_blaze) then begin
    rr_i = rr_i / $
      (1.0 - tan_blaze^2 *((l[hh[poss] / nl] + l[hh[poss] mod nl])/cw_i - 2.0))
  endif

  if nposs GT 1000 then begin

    min_mp = min(mp_i)
    max_mp = max(mp_i)
    bin_mp = (max_mp - min_mp)/inc
    small_his = hist_2d(rr_i, mp_i, min1=rr_pl - (inc+1)*bin_rr, $
                   max1 = rr_pl + (inc+1)*bin_rr, bin1=bin_rr, $
                   bin2 = bin_mp, min2=min_mp, max2=max_mp)
    small_sm = convol(small_his, fltarr(3,3)+1)
 
    rrs = (findgen((size(small_sm))[1]) +0.5) * bin_rr + rr_pl - (inc+1)*bin_rr 
    a = ladfit(rrs, total(small_sm,2)/total(small_sm))

    avg = (total(small_sm,1)) ## poly(rrs, a)
    small_sm = small_sm - avg
    ny = (size(small_sm))[2]
    xcen  = fltarr(ny)-1.
    xivar = fltarr(ny)-1.
    yhist = total(small_sm,1)/(size(small_sm))[1]
    for i=0,ny-1 do begin & $
      xcen[i] = find_npeaks(small_sm[*,i], nfind=1, ypeak=y) & $
      xivar[i] = ((y-yhist[i]) > 0)^2/(y+1) * (y GT 2) & $
    endfor

    yt = findgen(ny)
    xy2traceset, yt, xcen, tset, ncoeff=ncoeff, invvar=xivar, yfit=xfit, $
                upper=10, lower=10, maxrej=ny/3. , /silent
 
    mask = abs(xcen-xfit) LT inc/2.
    xy2traceset, yt, xcen, tset, ncoeff=ncoeff, invvar=xivar*mask, yfit=xfit, $
                upper=10, lower=10, maxrej=ny/3. , /silent
    

    yall = (mp_i - min_mp)/bin_mp - 0.5
    traceset2xy, tset, yall, xall
    x = (rr_i - (rr_pl - (inc+1)*bin_rr))/(bin_rr) - 0.5
 
    keep = where(x GE xall - 8 AND x LE xall + 8, nkeep)
  endif else begin
    nkeep = nposs
    keep = lindgen(nposs)
  endelse

    splog, 'Zeroing in on ', nkeep,' pairs out of ', n
  line_p = lonarr(np)-1
  num_p = lonarr(np)

  if nkeep LT np then begin
    splog, 'WARNING: Only found ',  nkeep, ' pairs'
    return, line_p
  endif

  total_hits = histogram(i1[poss[keep]], min=0, max=np-1) + $
               histogram(i2[poss[keep]], min=0, max=np-1) 

  order_hits = reverse(sort(total_hits))
  for ii=0,np-1 do begin
    i = order_hits[ii]

    a = where(i1[poss[keep]] EQ i)
    b = where(i2[poss[keep]] EQ i)
    line_hist = lonarr(nl)
    if a[0] NE -1 then $
      line_hist = histogram(line_i[keep[a]],min=0, max=nl-1)
    if b[0] NE -1 then $
      line_hist = line_hist + histogram(line_i[keep[b]+nposs],min=0, max=nl-1)

    am = max(line_hist,pl)
    if total(line_p EQ pl) eq 0 AND am GT 1 then begin
      line_p[i] = pl
      num_p[i] = am
    endif

  endfor

  splog, 'Finished: Found ', long(total(line_p GE 0)), $
         ' line matches out of ', np

; This is a hack. The logic is flawed somewhere and this nonmonotonic
; solutions should be fixed
  igood = WHERE(line_p GE 0)
  xpix = peaks[igood]
  wave = l[line_p[igood]]
  IF disp0 LT 0 THEN BEGIN
      dlambda = [-1e6, (wave - shift(wave, 1))[1:*]]
      nonmono = where(dlambda gt 0) 
  ENDIF ELSE BEGIN
      dlambda = [1e6, (wave - shift(wave, 1))[1:*]]
      nonmono = where(dlambda lt 0) 
  ENDELSE
  IF nonmono[0] NE -1 THEN line_p[igood[nonmono]] = -1
  
  return, line_p

end


