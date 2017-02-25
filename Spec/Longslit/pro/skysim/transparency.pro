function transparency, waves

  openr, 10, getenv('XIDL_DIR') + '/Spec/Longslit/pro/skysim/atm_transmission_secz1.5_1.6mm.dat'

   tmp1=0.
   tmp2=0.

   ntrans=0

   while (NOT EOF(10))do begin
       readf, 10, tmp1, tmp2
       ntrans++
   endwhile

   close, 10

   wave = fltarr(ntrans)
   trans = fltarr(ntrans)

   openr, 10, getenv('XIDL_DIR') + '/Spec/Longslit/pro/skysim/atm_transmission_secz1.5_1.6mm.dat'
   i=0
   while (NOT EOF(10))do begin
       readf, 10, tmp1, tmp2
       wave[i] = tmp1
       trans[i] = tmp2
       i++
   endwhile

   close, 10

   trans = trans[where(wave GT 0.8 AND wave LT 2.6)]
   wave = wave[where(wave GT 0.8 AND wave LT 2.6)]

   sset = bspline_iterfit(wave, trans, everyn=1.2)

   ans = bspline_valu(waves, sset)
   ans[where(waves LT 0.9)] = 1.0

   return, ans

end
