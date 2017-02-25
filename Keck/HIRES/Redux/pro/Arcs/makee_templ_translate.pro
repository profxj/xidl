pro makee_templ_translate,inlist,echnumlist,outfil

  readcol,inlist,n,name,ech,xd,col,l1,f1,tt,format='l,a,f,f,a,a,a,a'
  readcol,echnumlist,ewcen,eordr,format='f,l'

  openw,1,outfil
  
  ;read in each member of the list, grab header, return dimensions
  for i=0, n_elements(name)-1 do begin
     tmp=xmrdfits(name[i]+'.fits',0,hdr)
     hires_rdmakeewave,hdr,wave
     sz=size(wave,/dimensions)
    ;note, all of these are the *central* wavelegths
     wmin=wave[floor(sz[0]/2),0]
     wmax=wave[floor(sz[0]/2),sz[1]-1]
     wcen=wave[floor(sz[0]/2),floor(sz[1]/2)]

    ;translate central wavelengths to orders
     wmn=min(abs(wmin-ewcen),indx)
     owmn=eordr[indx]
     wmx=min(abs(wmax-ewcen),indx)
     owmx=eordr[indx]
     wmc=min(abs(wcen-ewcen),indx)
     owmc=eordr[indx]

     printf,1,name[i]+' '+string(ech[i])+'  '+string(xd[i])+'  '+string(owmn)+'  '+$
            string(owmc)+'  '+string(owmx)+'  '+string(wmin)+'  '+string(wcen)+'  '+string(wmax)

  endfor
  close,1
return
end
