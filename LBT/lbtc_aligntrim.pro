;to run above the Sci folder

;test --> compute the shifts, but not apply them

pro lbtc_aligntrim, sci_list, test=test, xtrim=xtrim, ytrim=ytrim

  n_sci=n_elements(sci_list)
  
  filter=strarr(n_sci)
  target=strarr(n_sci)

  ;reconstruct names
  for i=0, n_sci-1 do begin
     p1=strpos(sci_list[i],'_')
     p2=strpos(sci_list[i],'_sci')
     target[i]=strmid(sci_list[i],0,p1)
     filter[i]=strmid(sci_list[i],p1+1,p2-p1-1)
  endfor
  root=strmid(sci_list[i-1],p2+4,6)
  

    ;now compute shifts
  splog, 'Align frames...'
 
  lbtc_manual_align,  sci_list, xref_star=xref_star,yref_star=yref_star,$
                       xshift=outx, yshift=outy
  
 
  if keyword_set (test) then return
    
 
 ;loop and apply shift
  
  if ~keyword_set(xtrim) then xtrim=100
  if ~keyword_set(ytrim) then ytrim=100

  for st=0, n_sci-1 do begin
     splog, 'Shift ', sci_list[st]
     
   
        ;science
     fits=mrdfits(sci_list[st],0,header,/silent)
     fits=shiftf(fits,outx[st],outy[st])
        ;trim
     sizefits=size(fits,/dimensions)
     fits=fits[xtrim:sizefits[0]-xtrim,ytrim:sizefits[1]-ytrim]
     mwrfits, fits, sci_list[st], header, /create
        
        ;median
     medname=strcompress(target[st]+'_'+filter[st]+'_median'+root,/remove_all)
     
     fits=mrdfits(medname,0,header,/silent)
     fits=shiftf(fits,outx[st],outy[st])
        ;trim
     fits=fits[xtrim:sizefits[0]-xtrim,ytrim:sizefits[1]-ytrim]
     mwrfits, fits, medname, header, /create
      
        ;mask
     maskname=strcompress(target[st]+'_'+filter[st]+'_mask'+root,/remove_all)
     fits=mrdfits(maskname,0,header,/silent)
     fits=shiftf(fits,outx[st],outy[st])
        ;trim
     fits=fits[xtrim:sizefits[0]-xtrim,ytrim:sizefits[1]-ytrim]
     mwrfits, fits, maskname, header, /create
        
  endfor


end
