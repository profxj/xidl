function makee_picktemplate,echang,xdang,col

;expects col=0 for oldest, col=1 for uv, col=2 for red (all pre
;mosaic)


  if (col eq 0) then reflist='adrd97_translate.list'
  if (col eq 1) then reflist='aduv_translate.list'
  if (col eq 2) then reflist='adrd_translate.list'

  print,'reading from '+reflist
  ;read in reference files
  readcol,reflist,template,rech,rxd,omn,ocen,omx,wmn,wcen,wmx,format='a,f,f,i,i,i,f,f,f'

  ;find minimum deviate in xdang first, then determine echelle

  xdmin=min(abs(xdang-rxd),xdindx)
  range=where(abs(rxd-rxd[xdindx]) LT 0.05)
  
  echmin=min(abs(echang-rech[range]),echindx)


out=template[range[echindx]]+'.fits'
print,'guessing template is '+out
return,out

end
