pro flatfig
   loadct, 0
   rflat = mrdfits ('data/Flat_R_01_M.fits.gz')
   rflatavg = mrdfits ('data/rflat.fits')
   a = rflatavg[300:600,850:1150]
   b = rflat[300:600,850:1150]
   c = a/b
   b = 13930. + (b-1.)*300000.
   print, avg(a), avg(b), avg(c)
   print, stddev(a), stddev(b), stddev(c)
;   tvscl, c
;   stop

;   tvscl, b
;   stop

   img = [a,b,c]
;   tvscl, (img>7000)<18000  
;   stop

   imgb =  img 
   imgb[0,*]=0                  ; 0-300,301-601,602-902
   imgb[300,*] = 0
   imgb[300:301,*]=0
   imgb[601:602,*]=0
   imgb[901:902,*] = 0
   x_psopen, 'redflatdemo.ps', /portrait, /encap
   tvscl, ((imgb>7000)<22000  )
;   tvscl, 18000 - ((imgb>8000)<18000  )
   x_psclose
   ;spawn, 'gv redflatdemo.ps'
end



