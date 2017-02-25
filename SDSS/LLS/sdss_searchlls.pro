;+ 
; NAME:
; sdss_searchlls
;
; PURPOSE:
;    Search for LLS in the SDSS QSO distribution
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2005 Written by JO + JXP
;-
;------------------------------------------------------------------------------
;IDL> sdss_searchlls, 'dr3_qso.fits', 'foo.fits'
;sdss_searchlls,'dr3_qso.fits', 'hiztemplate.fits', $
;'foo.fits',szstart=3.5,szend=3.53, allmod=allmod
PRO sdss_searchlls, DR_fits, hiztemplate, $
                    outfil, szstart=SZSTART, szend=SZEND, $
                    drpath=DRPATH, ALLMOD=allmod

;set the DR path
   IF NOT keyword_set(drpath) THEN drpath=getenv('SDSSPATH')+'/DR3_QSO/'

;open the dr summary fits file
   dr=xmrdfits(drpath+DR_fits,1,/silent)

;determine the number of qsos to search
   IF NOT keyword_set (szstart) THEN szstart = 3.19
   IF NOT keyword_set (szend) THEN szend = max(dr.z)
   
   range = where(dr.z GE szstart AND dr.z LE szend,nqsos)

   
;make the big structure

   tmp = {sdssllsstrct}
   llssearch = replicate(tmp, nqsos)

;fill up the non-lls dependent arrays
   llssearch.plate=dr[range].plate
   llssearch.mjd=dr[range].mjd
   llssearch.fiber=dr[range].fiberid
   llssearch.zem=dr[range].z

;now for the big loop
   ;; Read in HI
    nnhi = 20
    nh1 = xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',0)
    npnh = n_elements(nh1)
    allnh = fltarr(npnh,nnhi)
    for i=0L,nnhi-1 do allnh[*,i] = xmrdfits(getenv('XIDL_DIR')+ $
                                             '/SDSS/LLS/nhi16_19b30.fits',i)


   ;; Create model image
   nztest = long((alog10(911.7633*6) - alog10(911.7833*4.2))/0.0001) + 11
   zllsarr = (10^(findgen(nztest)*0.0001 + alog10(4.2*911.7633)))/911.7633 - 1.0
   wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
   mn = min(abs(wv_mod-911.7633),mmn)

   wv_sdss = 10^(3.58 + dindgen(3900)*1e-4)
   mn = min(abs(4.2*911.7633-wv_sdss),imn)
   shft = imn-mmn

   if not keyword_set(ALLMOD) then begin
       allmod = fltarr(npnh, nztest, nnhi)
       for rr=0l,nnhi-1 do begin
           for qq=0L,nztest-1 do begin
               allmod[*,qq,rr] = shift(allnh[*,rr], shft+qq)
               allmod[(npnh+qq+shft)<(npnh-1):*,qq,rr] = 1.
           endfor
       endfor
   endif
   ;; Template
   template=xmrdfits(hiztemplate,/silent)

   ;; Loop
   FOR qq=0,nqsos-1 do begin
      ;make string verisons of the plate,mjd,fiber
       IF(llssearch[qq].plate LT 1000) THEN BEGIN
           splate=strcompress('0'+string(llssearch[qq].plate),/remove_all)
       ENDIF ELSE BEGIN
           splate=strcompress(string(llssearch[qq].plate))
       ENDELSE

       smjd=string(llssearch[qq].mjd, format='(i5.5)')
       
       sfiber=string(llssearch[qq].fiber, format='(i3.3)')

      ;create the file name 
       filename = strcompress('spSpec-'+smjd+'-'+splate+'-'+sfiber+'.fit.gz', $
                              /remove_all)
       filedir = strcompress(drpath+'spectro/1d_25/'+splate+'/1d/', $
                             /remove_all)
       file=filedir+filename
       

       ;now run the lls finder
       print,'Now doing qso '+filename+' (#'+ $
         strcompress(string(qq),/remove_all)+' of '+ $
         strcompress(string(nqsos-1),/remove_all)+')'
       gdmod = where(zllsarr LT llssearch[qq].zem+0.1)
       llsstate=sdss_findlls(template,2.8,file,/plotfit,/plotchi, $
                            allmod=allmod[*,gdmod,*])

       ;fill in the remaining variables into the structure
       llssearch[qq].llsflg=llsstate.llsflg
       llssearch[qq].zlls=llsstate.zlls
       llssearch[qq].taulls=llsstate.taulls
       llssearch[qq].blls=llsstate.blls
       llssearch[qq].zstart=llsstate.zstart
       llssearch[qq].zend=llsstate.zend
       llssearch[qq].uminusg=llsstate.uminusg
       llssearch[qq].svchi=llsstate.svchi
       llssearch[qq].svz=llsstate.svz
       ;stop
   ENDFOR
   mwrfits, llssearch, outfil, /create

return
END

