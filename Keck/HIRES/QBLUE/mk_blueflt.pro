pro mk_blueflt, list, outfil

  if  N_params() LT 2  then begin
      print,'Syntax - ' + $
        'xcombine, img, [cmbimg, header], FCOMB=, MMEM=, MASKS=, SIGLO=' 
      return
  endif

 ;; Files
 readcol, list, flats, format='a'

 ;; Bias subtract
 nfil = n_elements(flats)
 for qq=0L,nfil-1 do hires_subbias, flats[qq]

 ;; Combine
 pos = strpos(flats[0],'/')
 pos = pos > 0
 ovflt = 'OV/'+strmid(flats,pos)

 xcombine, ovflt, finimg, fcomb=2  ; Median + rej

 ;; Output
 mwrfits, finimg, outfil, /create
 spawn, 'gzip -f '+outfil

 return
end

