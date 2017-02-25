;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_group.pro               
; Author: Kathy Cooksey                      Date: 29 Jan 2009
; Project: 
; Description: 
; Input: 
; Optional:
; Output: 
; Example:
; History:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_calcewn

pro civ_group,subciv_fil, civstrct_fil, rating=rating,flg_sys=flg_sys,$
              zlim=zlim,ewlim=ewlim,nlim=nlim, unsat=unsat, nciv=nciv, $
              siiv=siiv, dvgal=dvgal, dvqso=dvqso

  if not keyword_set(civstrct_fil) then begin
     if keyword_set(siiv) then $
        civstrct_fil = getenv('MLSS_DIR')+'/analysis/evalssiiv/siiv_rtgge5.fits' $
     else $
        civstrct_fil = getenv('MLSS_DIR')+'/analysis/evals/civ_rtgge5.fits'
  endif 

  ;; Read in file
  if size(civstrct_fil,/type) eq 7 then $
     civstrct = xmrdfits(civstrct_fil,1,/silent) $
  else civstrct = civstrct_fil
  nciv = n_elements(civstrct)

  ;; Assume first element define structure focus
  gd = where(strtrim(civstrct[0].ion,2) ne '')
  indx = gd[0]
  prs = strsplit(civstrct[0].ion[gd[0]],/extract,count=nprs) 
  dblt_name = prs[0]

  ;; Rating cut
  if keyword_set(rating) and nciv gt 0 then begin
     nrtg = n_elements(rating)
     for ii=0,nrtg-1 do begin
        gd = where(civstrct.rating_eye eq rating[ii],nciv)
        if nciv eq 0 then begin
           print,'civ_group: No '+dblt_name+' with rating = ',rating[ii]
        endif else begin
           if not keyword_set(substrct) then substrct = civstrct[gd] $
           else substrct = [substrct,civstrct[gd]]
        endelse 
     endfor                     ; loop rating
     if not keyword_set(substrct) then civstrct = -1 $
     else civstrct = substrct
     nciv = n_elements(civstrct)
  endif                         ; rating=

  ;; Flag cut
  if keyword_set(flg_sys) and nciv gt 0 then begin
     gd = where((civstrct.flg_sys[0] and flg_sys) ge flg_sys,nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' with flg_sys = ',flg_sys
        civstrct = -1
     endif else begin
        civstrct = civstrct[gd]
     endelse 
  endif                         ; flg_sys=

  ;; Galaxy buffer (> dvgal/c)
  if keyword_set(dvgal) and nciv gt 0 then begin
     gd = where(civstrct.zabs[indx] gt dvgal/2.998e5,nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' with z*c > ',dvgal
     endif else begin
        civstrct = civstrct[gd]  
     endelse 
  endif 

  ;; QSO buffer (< zqso-dvqso/c*(1+zqso))
  if keyword_set(dvqso) and nciv gt 0 then begin
     gd = where(civstrct.zabs[indx] lt civstrct.zqso-dvqso/2.998e5*$
                (1+civstrct.zqso),nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' with z*c > ',dvgal
     endif else begin
        civstrct = civstrct[gd]  
     endelse 
  endif 

  ;; Redshift cut
  if keyword_set(zlim) and nciv gt 0 then begin
     gd = where(civstrct.zabs[indx] ge zlim[0] and $
                civstrct.zabs[indx] lt zlim[1],nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' within z limits = ',zlim
        civstrct = -1
     endif else begin
        civstrct = civstrct[gd]
     endelse 
  endif                         ; zlim= 

  ;; EW cut
  if keyword_set(ewlim) and nciv gt 0 then begin
     gd = where(civstrct.ew[indx] ge ewlim[0] and $
                civstrct.ew[indx] lt ewlim[1],nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' within EW limits = ',ewlim
        civstrct = -1
     endif else begin
        civstrct = civstrct[gd]
     endelse 
  endif                         ; ewlim= 

  ;; Column density cut
  if keyword_set(nlim) and nciv gt 0 then begin
     ncolm = civ_calcewn_ndblt(civstrct,dblt_name,/log,$
                               flg_colm=flg_colm,/silent)
     gd = where(ncolm ge nlim[0] and ncolm lt nlim[1],nciv)
     if nciv eq 0 then begin
        print,'civ_group: No '+dblt_name+' within log N limits = ',nlim
        civstrct = -1
     endif else begin
        civstrct = civstrct[gd]
     endelse 
  endif                         ; nlim= 

  ;; Saturation cut
  if keyword_set(unsat) and nciv gt 0 then begin
     ncolm = civ_calcewn_ndblt(civstrct,dblt_name,/log,$
                               flg_colm=flg_colm,/silent)
     gd = where((flg_colm and 6) ne 2,nciv)
     if nciv eq 0 then begin
        print,'civ_group: No unsaturated '+dblt_name
        civstrct = -1
     endif else begin
        civstrct = civstrct[gd]
     endelse 
  endif                         ; nlim= 

  ;; Sort by RA
  unq = uniq(civstrct.qso,sort(civstrct.qso)) 
  nqso = n_elements(unq) 
  unq = unq[sort(civstrct[unq].ra)]
  
  for qq = 0,nqso-1 do begin 
     mtch = where(civstrct.qso eq civstrct[unq[qq]].qso) 
     srt = sort(civstrct[mtch].zabs[0]) 
     if qq eq 0 then tmp = civstrct[mtch[srt]] else $
        tmp = [tmp,civstrct[mtch[srt]]] 
  endfor
  if n_elements(tmp) ne nciv then $
     stop,'civ_group: sorting by RA has corrupted results'
  civstrct = tmp

  ;; Return results
  if keyword_set(subciv_fil) then begin
     if size(subciv_fil,/type) eq 7 then begin
        mwrfits,civstrct,subciv_fil,/create,/silent
        print,'civ_group: created ',subciv_fil
     endif else subciv_fil = civstrct
  endif else subciv_fil = civstrct

end
