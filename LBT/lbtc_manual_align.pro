
;perform an interactive alignment with positions of nstar 
;selected by the user

;name list of file names---> multidim array with the fits to align

pro lbtc_manual_align, name, xref_star=xref_star,yref_star=yref_star,$
                       xshift=xshift, yshift=yshift, chipnum=chipnum

  ;some into stuff
 
  nsci=n_elements(name)
  

 ;-----------------------------------------------------
 ;Provide a reference star position and do first shift
 ;-----------------------------------------------------
      

  
     splog, 'Select reference stars for first alignemt. Then, press right-button to continue...'

     
     ctload, 0, /reverse
     window, 0, xsize=600, ysize=600
     fits=mrdfits(name[0])
     ;reform and quick skysub
    if keyword_set(chipnum) then  $
       imgzero=reform(fits[chipnum,*,*])-djs_median(fits[chipnum,*,*]) else $
          imgzero=fits-djs_median(fits)
    undefine, fits
     imdisp, bytscl(imgzero,min=0.,max=1.), /axis
    
     
     ;get first
     cursor, x, y, /wait, /data 
     ;centroid
     cntrd, imgzero, x, y, xcen, ycen, 15.
     xref_tmp=xcen
     yref_tmp=ycen
     splog, 'Star Position ', xcen, ycen
                

     ;then continue
     for i=0, 100 do begin
        cursor, x, y, /down, /data 
        if(!mouse.button eq 1) then begin
           cntrd, imgzero, x, y, xcen, ycen, 15.
           xref_tmp=[xref_tmp,xcen]
           yref_tmp=[yref_tmp,ycen]
           splog, 'Star Position ', xcen, ycen
        endif
        if(!mouse.button eq 4) then break
     endfor
     
 
 
     nstar=n_elements(xref_tmp)
     
     xref_star=fltarr(nsci,nstar)
     yref_star=fltarr(nsci,nstar)
     xshift_mat=fltarr(nsci,nstar)
     yshift_mat=fltarr(nsci,nstar)
     
     ;set reference
     xref_star[0,*]=xref_tmp[*]
     yref_star[0,*]=yref_tmp[*]
     

     ;display reference star
     erase
     display,  bytscl(imgzero,min=0.,max=1.)
     undefine, imgzero

     for sta=0, nstar-1 do begin
        x_oplotcirc, 40., x0=xref_star[0,sta], y0=yref_star[0,sta]
        xyouts, xref_star[0,sta], yref_star[0,sta], string(sta+1)
     endfor

     

     
     ;now get other star positions
     window, 1, xsize=600, ysize=600
     splog, 'Now mark other  reference star'

     for sci=1, nsci-1 do begin
        splog, 'Load new image...'
        ctload, 0, /reverse
        fits=mrdfits(name[sci])
        if keyword_set(chipnum) then  $
           thisimg=reform(fits[chipnum,*,*])-djs_median(fits[chipnum,*,*]) else $
              thisimg=fits-djs_median(fits)
        undefine, fits
        imdisp, bytscl(thisimg,min=0.,max=1.), /axis
        
        for sta=0, nstar-1 do begin
           cursor, x, y, /down, /data 
           cntrd, thisimg, x, y, xcen, ycen, 15.
           xref_star[sci,sta]=xcen
           yref_star[sci,sta]=ycen
           splog, 'Star Position ', xcen, ycen
        endfor
        
        undefine, thisimg
        
        xshift_mat[sci,*]=xref_star[0,*]-xref_star[sci,*]
        yshift_mat[sci,*]=yref_star[0,*]-yref_star[sci,*]

     endfor
    
     ctload, 39
 
     ;close windows
     wdelete, 0
     wdelete, 1

     ;get final shift
     xshift=djs_median(xshift_mat,2)
     yshift=djs_median(yshift_mat,2)
     
     splog, 'Shift X ', xshift
     splog, 'Shift Y ', yshift

end
