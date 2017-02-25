
function fit2std, s1, ncoeff=ncoeff, mask=mask, wmin=wmin, wmax=wmax, $
       corr=corr

   if NOT keyword_set(ncoeff) then ncoeff=7
   if NOT keyword_set(corr) then corr=1.0
   ; first mask pixels in both spectra
 
   norder = n_elements(s1)
   res = 0.

   wmin = s1.wave[0]
   wmax = wmin*0.0
   lower = lonarr(norder)
   upper = lonarr(norder)
   tpix = long(total(s1.npix))
   o = lonarr(tpix)
   l = fltarr(tpix)
   livar = fltarr(tpix)
   m = lonarr(tpix)
   w = lonarr(tpix)
   x = fltarr(tpix)
   
   
   for qq=0,norder-1 do begin
      i1 = qq EQ 0 ? 0 : i2 + 1
      npix = s1[qq].npix
      i2 = i1 + npix - 1
      w1 = s1[qq].wave[0:npix-1]
      nrow = s1[qq].nrow

      corr1 = w1*0.0 + 1.0
      if n_elements(corr) EQ n_elements(s1.box_wv[0:nrow-1]) then $
        linterp, s1[qq].box_wv[0:nrow-1], corr[*,qq], w1, corr1
      wmax[qq] = s1[qq].wave[npix-1]
      fx1 = s1[qq].fx[0:npix-1]  
      logfx = alog(fx1 > 1) + alog(corr1 > 1.0e-8)

      linterp, s1[qq].box_wv[0:nrow-1], s1[qq].box_fx[0:nrow-1], $
                  w1, box1
      logfx_ivar = fx1 * (fx1 GT 1) *(corr1 GT 0) * (fx1 LT 2.0*box1)
      mask1 = 1
      if keyword_set(mask) then begin
        linterp, s1[qq].box_wv[0:nrow-1], mask[*,qq], w1, mask1
        mask1 = mask1 EQ 1
      endif

      o[i1:i2] = qq
      l[i1:i2] = logfx
      livar[i1:i2] = logfx_ivar
      w[i1:i2] = long(w1*10000.0)
      m[i1:i2] = mask1
      x[i1:i2] = 2.0d*(w1 - wmin[qq])/(wmax[qq]-wmin[qq])-1.0
    endfor

;
;  enlarge mask to all common wavelengths
;

   ws = sort(w)
   wu = uniq(w[ws])
   ii=0L
   for i=0L,n_elements(wu)-1 do begin & $
     m[ws[ii:wu[i]]] = total(1-m[ws[ii:wu[i]]]) EQ 0 & $
     ii = wu[i]+1 & $
   endfor


    livar = livar* (livar LT 1.0e5)

    fleg = flegendre(x, ncoeff)
    np = ncoeff*norder
    chi = 0.0
    bw = 2* ncoeff
    bi = lindgen(bw)
    bo = lindgen(bw)
    for i=1L, bw-1 do bi = [bi, lindgen(bw-i)+(bw+1)*i]
    for i=1L, bw-1 do bo = [bo, lindgen(bw-i)+bw*i]


  ymodel = l
  for kk=1,3 do begin
    m = m * (chi^2 LT 400.0) 
    livar = livar * (abs(l-ymodel) LT 0.5)
    alpha = dblarr(bw,np+bw)
    beta = dblarr(np+bw)
    aull = 0


    for qq=0,norder-1 do begin
;
;    First put in normalized fit
;
      t = where(o EQ qq, nt)

      a = fltarr(nt, bw)

      a[*, 0:ncoeff-1] = fleg[t,*] 

      sqinv = sqrt(livar[t]*m[t])
      ct = 0
      if qq LT norder-1 then begin
        t2 = where(o EQ qq+1, nt2)
        if nt2 GT 1 then match, w[t], w[t2], sub1, sub2, count = ct
        if ct GT 1 then sqinv[sub1] = sqinv[sub1] * m[t2[sub2]]
      endif

      f = a * (sqinv # replicate(1,2*ncoeff))
      d = l[t] * sqinv

    
      if keyword_set(aull) then begin
        lower[qq] = (size(aull))[1]
        aull = [aull, a]
      endif else aull = a 
      upper[qq] = lower[qq] + nt -1

;
;    Now try to find overlaps
;
      if (qq LT norder-1 AND ct GT 1) then begin


        diff = l[t[sub1]] - l[t2[sub2]]
        good = where(livar[t[sub1]] GT 0 AND livar[t2[sub2]] GT 0, ngood)
;        print, qq, ct, ngood
        if ngood LT 2  then break

        diff_ivar = 1.0/(1.0/livar[t[sub1[good]]] + 1.0/livar[t2[sub2[good]]])
        
        adiff = [[fleg[t[sub1[good]], *]], [-1.0*fleg[t2[sub2[good]], *]]] 
        fdiff = adiff * (0.5*sqrt(diff_ivar) # replicate(1, 2*ncoeff))

        aull = [aull, adiff]
        f = [f, fdiff]
        d = [d, diff[good]*0.5*sqrt(diff_ivar)]
        upper[qq] = (size(aull))[1] - 1


      endif 

      work = f ## transpose(f)
      wb =   d # f
      itop = qq*ncoeff
      ibottom = (itop < np + bw) - 1

      alpha[bo+itop*bw] = alpha[bo+itop*bw] + work[bi]
      beta[itop:ibottom] = beta[itop:ibottom] + wb

     
    endfor

    errb = cholesky_band(alpha, mininf=min_influence)
    res = beta
    errs = cholesky_solve(alpha, res)

    res = reform(res[0:np-1], ncoeff, norder)
    ymodel = l*0.0
    for qq=0,norder-1 do begin
      t = where(o EQ qq)
      ymodel[t] = fleg[t,*] # res[*,qq] 
    endfor

    chi = (l-ymodel)*sqrt(livar*m)
  endfor

return, res
end
