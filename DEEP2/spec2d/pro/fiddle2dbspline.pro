
function shift_model, model,  sky 
; test if offsetting by a fractional pixel shift leads to improved fit
; 
      xx =  [-2., -1., 0.,  1.,  2.]
      cc = dblarr(5)
      sz = (size(sky, /dimen))[1]
      shift = fltarr(sz)
      for i=0, sz-1 do begin 
        for j=0, 4 do cc[j] = total(shift(model[*, i], j-2)*sky[*, i]/1.d6)
        fit = poly_fit(xx, cc, 2)
        shift[i] = -fit[1]/(2.*fit[2])  
      endfor
return,  shift
end


;pro junk
; silly code to study how to deal with 2d-bsplines
; md 28sep
; 
;restore, file='2dtest.sav', /verbose

print, 'start point of slitlets: ', slitlet_ptr
read, k, prompt = 'input desired slitlet number: '

;s1 = 737
;s2 = 795
;s1 = 781 +10
;s2 = 851 -10
;s1 =  997 +10
;s2 = 1077 -10
s1 = slitlet_ptr[k] +10
s2 = slitlet_ptr[k+1] -10

;select single slitlet, at end of chip1
arcwave = all_arcwave[*, s1:s2 ] 
arcflux = all_arcflux[*, s1:s2 ]
arcivar = all_arcivar[*, s1:s2 ]
syarray = yarray[*, s1:s2 ]


sz = size(arcwave, /dimen)
;index = lindgen(sz[1])


minlambda = min(arcwave, max=imax) ;*1.001 ;shade the bkpts to leave no gap
maxlambda = max(arcwave, min=imin) ;*.999
mminlambda = max(arcwave[0, *])
mmaxlambda =  min(arcwave[4095, *]) ;set min and max so all rows contribute

;xx = replicate(1.,  sz[0])
ymin = min(syarray)
ymax = max(syarray)

nbkpts = 4096*2
 bkpt = findgen(nbkpts)*(maxlambda-minlambda)/(nbkpts-1) +minlambda 


;set breakpoints to be uniform in this space, but why are we in log
;lambda??

;   Is = convol(arcflux,[1.,1.,1.],3.,/center) 
;   Is = convol(arcflux,[1.,2.,1.],3.,/center) 
; spectrum smoothed slightly -- done already in previous routine
;   Is = all_arcflux
   isort = bsort(arcwave)
   L1 = arcwave(isort)
   I1 = arcflux(isort)
   ivvar =  arcivar(isort)

 tstart = systime(1)
 sset = bspline_iterfit(L1, I1, $ ; invvar=ivvar, $
;         x2=syarray ,npoly=2 , xmin=ymin, xmax=ymax, $ 
;         upper = 20, lower=20, $
         nord = 4,  maxiter=0, bkpt=bkpt, yfit=skyfit)  
 time1 = systime(1)
 print, 'time for bspline estimation:', time1-tstart 
; wsky = (findgen(100000)+0.5)/100000*(max(wave)-min(wave))+min(wave)
  model= bspline_valu(arcwave, sset); , x2=syarray)
 print, 'time for bspline evaluation:', systime(1)-time1 


;  splot, 10^arcwave, arcflux, ps=3
;  soplot,10^wsky,sky, color=2

;  dfpsplot, 'arcfit.ps', /sq
;  plot, 10^arcwave, arcflux, ps=3, xr=[8120, 8170], yr=[1000, 2500], $
;    xtit='wavelength [Ang]', ytit='counts', title='Arcfit', chars=1.5
;  oplot,10^wsky,sky, color=3, thick=3
;  xyouts, 8142, 1400, 'Ghost', /data, chars=1.5
;  dfpsclose
 
   err = arcflux-model
   atv, transpose( [transpose(err), transpose(model),  $ 
      transpose(arcflux) ]),  max=500, min=-500

     nrr =  (size(err, /dimen))[1]
     mask = convol(float(arcflux gt 30000), $
           [1., 1., 1., 1., 1.], 5., /center)
      xpos = findgen(4096)#replicate(1,nrr) 
      ypos = replicate(1, 4096)#findgen(nrr)
 
      j = where(arcflux gt 500 and mask lt .2) ; unsaturated lines
      dd = histogram(err[j]/sqrt(arcflux[j]/1.2), bin=.2, max=10, min=-10)
      xposs = xpos[j] ;selected x,y position
      yposs = ypos[j]
      delta = err[j]/sqrt(arcflux[j]/1.2)
      jj = where(delta gt 3, njj)
      if njj gt 0 then atvplot, xposs[jj], yposs[jj]+nrr, psym=1, color=3
      jj = where(delta lt  -3, njj)
      if njj gt 0 then atvplot, xposs[jj], yposs[jj]+nrr, psym=1, color=2

      xx = findgen(101)/5.-10.
      window, 1
      wset, 1
      plot,xx, dd,  psym=10, xtitle='error (sigma)', ytitle='number', $
       title='residuals for points with 500<I<30000, slitlet: '+ $
           string(k,  format= '(i3)' )
      gauss = exp(-.5*xx^2)
      gauss = gauss*dd[50]
      oplot, xx, gauss, lines=1

      shift = shift_model(model, arcflux)
      window, 2
      plot, shift




  ;arcmodel = bspline_valu(wave, sset)

;  mask = arcflux LT 60000
;  block = [[(arcflux - model)*mask],[arcmodel-1800],[ext-1800]]
;  img = block[1010:1521, *]

; ------------------------------------------------------------------------
; Plot 6 - Arc subtraction
;  dfpsplot, 'arcsub.ps', /sq, /color, bits=8
;  bim = bytscl(img,min=-1000,max=5000, top=245)+8B
;  bim2 = bytscl((ext-arcmodel)[1010:1521, *],min=-250,max=250, top=245)+8B
;  display,[[bim2], [bim]], title='arc; arcmodel; residual; res stretched', $
;    xtit='lambda [pix]', ytit=$
;    '[pix]', chars=1.5, xmargin=[7, 2]
;  dfpsclose


  ;atv,block, max=5000, min=-1000



end
