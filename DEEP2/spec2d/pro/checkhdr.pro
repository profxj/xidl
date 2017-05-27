pro checkhdr
  flist=findfile('d*.fits')
  for i=0,n_elements(flist)-1 do $
    print,flist[i],fix(total(strmid(headfits(flist[i]),0,1) ne ' '))

  return
end
