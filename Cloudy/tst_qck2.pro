pro tst_qck2, grid

  rto = [0.12, 0.44, 0.8]
  sig = [0.05, 0.05, 0.05]
  ion1 = [ [14,3], [14,4], [6,4] ]    ; SiIII/SiII, SiIV/SiIII, CIV/CII
  ion2 = [ [14,2], [14,3], [6,2] ]
  model = lonarr(2)
  model[0] = where(grid.NHI EQ 16.0d AND grid.U EQ -6.d AND $
               grid.nH EQ -3.d AND grid.FeH EQ -1.d)
  model[1] = where(grid.NHI EQ 16.0d AND grid.U EQ -0.2d AND $
               grid.nH EQ -3.d AND grid.FeH EQ -1.d)


  cldy_qck2, grid, rto, sig, ion1, ion2, model


return
end
  


  
