pro tst_2phase, grid

; All Si
  rto = [0.12, 0.44, 0.56]
  sig = [0.07, 0.07, 0.07]
  ion1 = [ [14,3], [14,4], [14,4] ]    
  ion2 = [ [14,2], [14,3], [14,2] ]

;  All Si + CIV/CII
;  rto = [0.12, 0.44, 0.56, 0.8]
;  sig = [0.07, 0.07, 0.07, 0.07]
;  ion1 = [ [14,3], [14,4], [14,4], [6,4] ]    ; SiIII/SiII, SiIV/SiIII, CIV/CII
;  ion2 = [ [14,2], [14,3], [14,2], [6,2] ]


  cldy_2phas, grid, rto, sig, ion1, ion2, soltn, $
    NHILMT=[15.1, 19.1], FEHLMT=[-1.1, -0.9], UONLY=-1.d


  i = soltn.min[0]
  j = soltn.min[1]
  ii = soltn.sub[i]
  jj = soltn.sub[j]
  
  print, 'Values: '
  print, soltn.val, rto
  print, 'Answer: ', soltn.ans[i,j], soltn.sigans[i,j]
  print, 'Chisq: ', soltn.chisq[i,j]
  print, 'NHI: ', grid[ii].NHI, grid[jj].NHI
  print, 'U: ', grid[ii].U, grid[jj].U

return
end
  


  
