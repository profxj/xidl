;+ 
; NAME:
; sdss_functions.pro
;    Version 1.0
;
; PURPOSE:
;   Many functions that retrieve all little information for/from SDSS
;   CIV/SiIV/CaII stuff
;
; CALLING SEQUENCE:
;   
;   sdss_functions will print all the possible functions to screen
;   (and compile everything).
;   Otherwise, compile at beginning of codes with:
;   @sdss_functions or by calling sdss_functions,/compile
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-May-2011  created by KLC
;   30-Jul-2011  updated to handle new structure, KLC
;   06-Sep-2016  updated sdss_calcsigpoiss() to give actual c.l., DRM via KLC
;                updated sdss_calcsigbinom() to give actual c.l., DRM via KLC
;   20-Dec-2016  sdss_getcivstrct(/noBAL) bug fix, KLC
;-
;------------------------------------------------------------------------------
;; Conventions:
;; *Every function/procedure must have n_params() check or /help
;; keyword. 
;; *Every function/procedure must have description in cade and in
;; sdss_functions print-out.
;; *Comments!  Including associating final end with
;; function/procedure. 
;; *Function/procedure names:  sdss_get[], sdss_prnt[], sdss_plt[],
;; sdss_cp[]2[]. 

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getsdssdir,help=help
  ;; Set the default data structure (with appropriate /'s)
  if keyword_set(help) then begin
     print,'Syntax - sdss_setsdssdir([/help])'
     return,-1
  endif 
     
  return,getenv('SDSSPATH')+'/'+getenv('SDSSDR')+'/'
end                             ; sdss_getsdssdir()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getname,sdsstab,strct=strct,spec=spec,hdr=hdr,$
                      eig=eig, spl=spl, hyb=hyb, extrap=extrap, $
                      abslin=abslin, clean=clean, gz=gz, mc=mc, $
                      quick=quick, user=user, plate=plate, root=root, dir=dir
  ;; Return the spSpec-jjjjj-pppp-fff*.fits for given structure or
  ;; header 
  if n_params() ne 1 then begin
     print,'Syntax - sdss_getname(sdsstab,[/spec,/strct,/hdr,/eig,/spl,/abslin,plate='
     print,'                      /hyb,/clean,/gz,mc=,/quick,user=,root=, dir=])'
     return,-1
  endif 
  
  ;; Parse pieces of names
  nfil = (size(sdsstab,/dim))[0]
  if keyword_set(spec) or keyword_set(strct) then begin
     if keyword_set(spec) then begin
        ;; Assumes all names are of equal length and nature; avoids looping
        prs = strsplit(sdsstab[0],'/',count=nprs)
        name = strmid(sdsstab,prs[nprs-1])
        prs = strsplit(name[0],'.',count=nprs)
        name = strmid(name,0,prs[1]-1) ; remove .fit*
        prs = strsplit(name[0],'-',count=nprs)
        mjd = strmid(name,prs[1],prs[2]-prs[1]-1)
        plate = strmid(name,prs[2],prs[3]-prs[2]-1)
        fiber = strmid(name,prs[3])
     endif else begin
        ;; QSO Name
        if size(sdsstab,/type) eq 7 then $
           civstr = xmrdfits(sdsstab,1,/silent) $
        else civstr = sdsstab
        prs = strsplit(civstr[0].qso_name,'-',count=nprs)
        root = civstr.qso_name
        mjd = strmid(root,0,5)
        plate = strmid(root,prs[1],4)
        fiber = strmid(root,prs[2],4)
     endelse 
  endif else begin

     if keyword_set(hdr) then begin
        ;; sdsstab is actually a header; and only 1
        nm = sxpar(sdsstab,'NAME')
        fiber = strtrim(sxpar(sdsstab,'FIBERID'),2)
        plate = strtrim(sxpar(sdsstab,'PLATEID'),2)
        prs = strsplit(nm,'-',/extract,count=nprs)
        mjd = prs[1]
     endif else begin
        ;; Use structure information
        if size(sdsstab,/type) ne 8 then $
           stop,'sdss_getname(): input must be structure'
        mjd = strtrim(sdsstab.smjd,2)
        fiber = strtrim(sdsstab.fiber,2)
        plate = strtrim(sdsstab.plate,2)
     endelse 
     
     ;; Set appropriate lengths
     bd = where(strlen(fiber) lt 3,nbd)
     while nbd ne 0 do begin
        fiber[bd] = '0'+fiber[bd]
        bd = where(strlen(fiber) lt 3,nbd)
     endwhile
     bd = where(strlen(plate) lt 4,nbd)
     while nbd ne 0 do begin
        plate[bd] = '0'+plate[bd]
        bd = where(strlen(plate) lt 4,nbd)
     endwhile
  endelse                       ; header or structure
  
  ;; Concatenate
  if not keyword_set(root) then root = mjd + '-' + plate + '-' + fiber
  if not keyword_set(name) then name = 'spSpec-' + root

  ;; dir is relative to $SDSSPATH/$SDSSDR
  dir = 'spectro/'              ; may be overwritten

  ;; Suffix
  if keyword_set(eig) and keyword_set(spl) and keyword_set(abslin) and keyword_set(hyb) then $
     stop,'sdss_getname(): cannot have all keywords (/eig, /spl, /abslin, /hyb) set'
  if keyword_set(eig) then begin
     dir = 'conti/'
     if keyword_set(extrap) then tmp = '-eigxconti' $ ; extrapolated
     else tmp = '-eigconti'
     name = name+tmp
  endif 
  if keyword_set(spl) then begin
     dir = 'conti/'
     if keyword_set(extrap) then tmp = '-splxconti' $ ; extrapolated
     else tmp = '-splconti'
     name = name+tmp
  endif 
  if keyword_set(hyb) then begin
     dir = 'conti/'
     if keyword_set(extrap) then tmp = '-hybxconti' $ ; extrapolated
     else tmp = '-hybconti'
     name = name+tmp
  endif 
  if keyword_set(abslin) then begin
     dir = 'abslin/'
     if keyword_set(extrap) then tmp = '-abslinx' $ ; extrapolated
     else tmp = '-abslin'
     name = name+tmp
  endif 
  if keyword_set(clean) then begin
     dir = 'cleanspec/'
     name = name+'-clean'
  endif 
  if keyword_set(gz) then begin
     tmp = 'gz'
     if size(gz,/type) eq 7 then $
        tmpsfx = strtrim(strlowcase(gz),2) $
     else tmpsfx = ''
     dir = tmp+tmpsfx+'/'
     name = name+'-' +tmp+ tmpsfx
  endif 
;  if keyword_set(gz) then begin
;     dir = 'completeness/'
;     name = name+'-cmplt'
;  endif 
  if keyword_set(mc) then begin
     tmp = 'mc'
     if  keyword_set(quick) then $
        tmpnam = tmp $
     else begin
        if keyword_set(user) then tmpnam = 'user' $
        else tmpnam = 'conti'
     endelse 
     if size(mc,/type) eq 7 then $
        tmpsfx = strtrim(strlowcase(mc),2) $
     else tmpsfx = ''
     dir = tmp+tmpsfx+'/'
     name = name+'-'+tmpnam + tmpsfx
  endif 
  if keyword_set(user) then begin
     ;; User will be a number and can append it to any file name
     usersfx = strtrim(user,2)
     bd = where(strlen(usersfx) lt 3) ; 100's
     while bd[0] ne -1 do begin
        usersfx[bd] = '0' + usersfx[bd]
        bd = where(strlen(usersfx) lt 3)
     endwhile 
     name = name + '-' +usersfx
  endif 

  ;; All have same subdirectory
  dir = dir + '1d_26/'+plate+'/1d/' 

  ;; Default
  return,name+'.fit'
  
end                             ; sdss_getname()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getqsostrct,dr=dr, name=name, bal=bal, nobal=nobal, help=help
  ;; Return the appropriate data release catalog structure or file
  ;; name to be used 
  if keyword_set(help) then begin 
     print,'Syntax - sdss_getqsostrct([dr=, /name, /bal, /nobal, /help])'
     return, -1 
  endif 

  ;; Set Parameters
  if not keyword_set(dr) then dr = 7 ;data release
  
  case dr of 
     7: if keyword_set(bal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_BAL.fit' $
     else if keyword_set(nobal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_noBAL.fit' $
     else fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit'
     else: begin
        print,'sdss_getqsostrct(): DR='+strtrim(dr,2)+' not supported'
        return, -1
     end
  endcase
        
  ;; Check existence
  test = file_search(fil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_getqsostrct(): file DNE ',strtrim(fil,2)+'*'

  ;; Return either name or read
  if keyword_set(name) then return,fil $
  else begin
     strct = xmrdfits(fil,1,/silent)
     return, strct
  endelse 
  
end                             ; sdss_getqsostrct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getsnrstrct, dr=dr, name=name, bal=bal, nobal=nobal, help=help
  ;; Return the appropriate data release catalog structure or file
  ;; name to be used 
  if keyword_set(help) then begin 
     print,'Syntax - sdss_getsnrstrct([dr=, /name, /bal, /nobal, /help])'
     return, -1 
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Set Parameters
  if not keyword_set(dr) then dr = 7 ;data release
  
  case dr of 
     7: if keyword_set(bal) then $
        fil = sdssdir+'inputs/dr7qso_BAL_SNR.fit' $
     else if keyword_set(nobal) then $
        fil = sdssdir+'inputs/dr7qso_noBAL_SNR.fit' $
     else fil = sdssdir+'inputs/dr7qso_srt_SNR.fit'
     else: begin
        print,'sdss_getsnrstrct(): DR='+strtrim(dr,2)+' not supported'
        return, -1
     end
  endcase
        
  ;; Check existence
  test = file_search(fil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_getsnrstrct(): file DNE ',strtrim(fil,2)+'*'

  ;; Return either name or read
  if keyword_set(name) then return,fil $
  else begin
     strct = xmrdfits(fil,1,/silent)
     return, strct
  endelse 
  
end                             ; sdss_getsnrstrct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getqsolist,dr=dr, read=read, help=help, count=count, bal=bal, $
                         nobal=nobal, abslin_dir=abslin_dir 
  ;; Return the appropriate data release catalog spectra list to be used 
  if keyword_set(help) then begin 
     print,'Syntax - sdss_getqsolist([dr=, /read, /help, count=, /bal, '
     print,'                         /nobal, abslin_dir=])'
     return, -1 
  endif 

  ;; Set Parameters
  if not keyword_set(dr) then dr = 7 ;data release
  
  case dr of 
     7: if keyword_set(bal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_BAL.list' $
     else if keyword_set(nobal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_noBAL.list' $
     else  fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list'
     else: begin
        print,'sdss_getqsolist(): DR='+strtrim(dr,2)+' not supported'
        return, -1
     end
  endcase
        
  ;; Check existence
  test = file_search(fil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_getqsolist(): file DNE ',strtrim(fil,2)+'*'
  ;; Return either name or read
  if keyword_set(read) then begin
     readcol,fil,nam,format='a',/silent
     count = (size(nam,/dim))[0]
     abslin_dir = nam[0]
     count = count-1
     return,nam[1:count]
  endif else begin
     return, fil
  endelse 
  
end                             ; sdss_getqsolist()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_wrqsolist,sdsssum,list_fil,absdir=absdir,silent=silent
  ;; Write the specifically formatted files to work with
  ;; e.g. sdss_fndlin
  if n_params() ne 2 then begin
     print,'Syntax - sdss_wrqsolist, sdsssum, list_fil, [absdir=,/silent]'
     return
  endif
  if not keyword_set(absdir) then absdir = 'abslin/' 
  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)

  spec_fil = sdss_getname(sdsstab,plate=plate,dir=dir)

  openw,1,list_fil
  printf,1,absdir
  writecol,list_fil,dir+spec_fil,filnum=1
  close,1
  if not keyword_set(silent) then $
     print,'sdss_wrqsolist: created ',list_fil

end                             ; sdss_wrqsolist()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_getqsoinlist,list_fil, sdsssum, snrstrct_fil, snr=snr,$
                      silent=silent,_extra=extra
  ;; Read list and find the matching QSOs
  if n_params() lt 2 then begin
     print,'Syntax - sdss_getqsoinlist, list_fil, sdsssum, snrstrct_fil, [/silent,'
     print,'                            /snr, _extra=]'
     return
  endif 

  readcol,list_fil,spec,format='a',skip=1,/silent
  nspec = (size(spec,/dim))[0]
  tmp = sdss_getname(spec,/spec,root=subqso_name)

  ;; _extra= includes /noBAL, /BAL
  sdsstab = sdss_getqsostrct(_extra=extra)
  tmp = sdss_getname(sdsstab,root=qso_name)

  if keyword_set(snr) or keyword_set(snrstrct_fil) then begin
     if keyword_set(snr) then begin
        if size(snr,/type) eq 8 then snrstrct = snr $
        else snrstrct = sdss_getsnrstrct(_extra=extra)
     endif
  endif
  if keyword_set(snr) then begin
     sdsstab = snrstrct
     qso_name = snrstrct.qso_name
  endif 

  ;; Loop
  for ii=0L,nspec-1 do begin
     mtch = where(subqso_name[ii] eq qso_name,nmtch)
     if nmtch ne 1 then $
        stop,'sdss_getqsoinlist: multiple matches!'
     if ii eq 0 then subsdsstab = sdsstab[mtch[0]] $
     else subsdsstab = [subsdsstab,sdsstab[mtch[0]]]

     if keyword_set(snrstrct_fil) then begin
        mtch = where(subqso_name[ii] eq snrstrct.qso_name,nmtch)
        if nmtch ne 1 then $
           stop,'sdss_getqsoinlist: multiple matches in S/N structure!'
        if ii eq 0 then subsnrstrct = snrstrct[mtch[0]] $
        else subsnrstrct = [subsnrstrct,snrstrct[mtch[0]]]
     endif 
  endfor                        ; loop ii=nspec

  if size(sdsssum,/type) eq 7 then begin
     mwrfits,subsdsstab,sdsssum,/create,/silent
     spawn,'gzip -f '+sdsssum
     if not keyword_set(silent) then $
        print,'sdss_getqsoinlist: created ',sdsssum
  endif else sdsssum = subsdsstab

  if keyword_set(snrstrct_fil) then begin
     if size(snrstrct_fil,/type) eq 7 then begin
        mwrfits,subsnrstrct,snrstrct_fil,/create,/silent
        spawn,'gzip -f '+snrstrct_fil
        if not keyword_set(silent) then $
           print,'sdss_getqsoinlist: created ',snrstrct_fil
     endif else snrstrct_fil = subsnrstrct
  endif 

end                             ; sdss_getqsoinlist


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_exclbal, list_fil, ballist_fil, olist_fil, sdsssum=sdsssum, $
                  clobber=clobber,_extra=extra
  ;; Take ballist_fil spectra out of list_fil and write to olist_fil
  ;; and optionally, sdsssum
  if n_params() ne 3 then begin
     print,'Syntax - sdss_exclbal, list_fil, ballist_fil, olist_fil, '
     print,'                       [/clobber, sdsssum=, _extra=]'
  endif 
  
  test = file_search(olist_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_exclbal: output list exists; will not overwrite; exiting'
     return
  endif 
  
  readcol,list_fil,spec_fil0,format='a',/silent
  absdir = spec_fil0[0]         ; save for output
  spec_fil0 = spec_fil0[1:*]
  nfil0 = (size(spec_fil0,/dim))[0] > 1
  
  readcol,ballist_fil,bal_fil,format='a',/silent,skip=1
  nbal = (size(bal_fil,/dim))[0] > 1
  
  badmask = intarr(nfil0)       ; 1 = bad; 0 = good

  for bb=0L,nbal-1 do begin
     mtch = where(bal_fil[bb] eq spec_fil0,nmtch)
     if nmtch gt 1 then $
        stop,'sdss_exclbal stop: duplicate matches'
     if nmtch ne 0 then badmask[mtch] = 1
  endfor                        ; loop bb=nbal

  gd = where(badmask eq 0,nfil)
  if nfil eq 0 then begin
     print,'sdss_exclbal: no unmatched spectra left; exiting'
     return                     ; EXIT
  endif 

  print,'sdss_exclbal: Old and new number of spectra:',nfil0,nfil

  spec_fil = spec_fil0[gd]
  openw,1,olist_fil
  printf,1,absdir
  writecol,olist_fil,spec_fil,fmt='(a)',filnum=1
  close,1
  print,'sdss_exclbal: created ',olist_fil

  if keyword_set(sdsssum) then begin
     test = file_search(sdsssum+'*',count=ntest)
     if ntest eq 0 or keyword_set(clobber) then $
        ;; _extra includes /silent, snrstrct=
        sdss_getqsoinlist,olist_fil,sdsssum,_extra=extra
  endif 

end                             ; sdss_exclbal


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getcflg, eig=eig, spl=spl, hyb=hyb, custom=custom, help=help, $
                       index=index
  ;; Return convention for whether absorption system defined with
  ;; eigen continuum or spline continuum.
  ;; eigen will be default
  ;; To map to index, take fix(alog(cflg)/alog(2))
  if keyword_set(help) then begin
     print,'Syntax - sdss_getcflg([/eig,/spl,/custom,/help,/index]); default /eig'
     return,-1
  endif 
  cflg = 1                             ; /eig and default; index = 0
  if keyword_set(custom) then cflg = 8 ; index = 3
  if keyword_set(hyb) then cflg = 4    ; index = 2
  if keyword_set(spl) then cflg = 2    ; index = 1

  if keyword_set(index) then return,fix(alog(cflg)/alog(2)) $
  else return, cflg
end                             ; sdss_getcflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getrating, definite=definite, good=good, maybe=maybe, bad=bad, $
                         unrated=unrated, help=help
  ;; Return convention sdss_chkciv evaluation
  if keyword_set(help) then begin
     print,'Syntax - sdss_getrating([/definite,/good,/maybe,/bad,/unrated])'
     print,'                  Default: /bad'
     return,-1
  endif 
  flg = 0
  
  if keyword_set(definite) and (keyword_set(good) or keyword_set(maybe) or $
                                keyword_set(unrated)) then begin
     print,'sdss_getrating(): multiple keywords set; defaulting'
     return,-1
  endif

  if not (keyword_set(definite) or keyword_set(good) or keyword_set(maybe) or $
          keyword_set(unrated)) then return, 0
  
  if keyword_set(definite) then return, 3
  if keyword_set(good) then return, 2
  if keyword_set(maybe) then return, 1
  if keyword_set(unrated) then return, -1
  
end                             ; sdss_getrating()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_srtcivstrct, civstrct_fil, by_rating=by_rating, index=index
  ;; Sort entries by QSO name (jjjjj-pppp-fff) and zabs
  if n_params() ne 1 then begin
     print,'Syntax - sdss_srtcivstrct( civstrct_fil, [/by_rating, index=])'
     return,-1
  endif 
  
  if size(civstrct_fil,/type) eq 8 then civstrct = civstrct_fil $
  else civstrct = xmrdfits(civstrct_fil,1,/silent)

  if keyword_set(by_rating) then begin
     index = civstrct.rating[0] * 0 - 1 ; quick instantiation
     rtg = civstrct[uniq(civstrct.rating[0],sort(civstrct.rating[0]))].rating[0]
     nrtg = (size(rtg,/dim))[0] > 1 ; foils singularity
     for rr=0,nrtg-1 do begin
        sub = where(civstrct.rating[0] eq rtg[rr])
        tmp = sdss_srtcivstrct(civstrct[sub],index=isub)
        if rr eq 0 then newcivstr = tmp $ 
        else newcivstr = [newcivstr,tmp]
        index[sub] = sub[isub]
     endfor                     ; loop rr=nrtg
     civstrct = newcivstr
  endif else begin
     civstrct.qso_name = strtrim(civstrct.qso_name,2)
     index = sort(civstrct.qso_name)
     civstrct = civstrct[index]
     unq = uniq(civstrct.qso_name)
     nunq = (size(unq,/dim))[0] 
     
     rng = lindgen(unq[0]+1)
     srt = sort(civstrct[rng].zabs_orig[0])
     civstrct[rng] = civstrct[rng[srt]] ; sort by QSO and z
     index[rng] = index[rng[srt]]
     
     for qq=1L,nunq-1 do begin
        rng = unq[qq-1] + 1 + lindgen(unq[qq]-unq[qq-1])
        srt = sort(civstrct[rng].zabs_orig[0])
        civstrct[rng] = civstrct[rng[srt]] ; sort by QSO and z
        index[rng] = index[rng[srt]]
     endfor                                ; loop qq=nunq
  endelse

  return, civstrct
end                             ; sdss_srtcivstrct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_rmcivduplicates, civstrct_fil, count=count, index=index, $
                               complement=complement, icomplement=icomplement, $
                               ncomplement=ncomplement, check=check
  ;; Remove duplicate candidates based on qso_name, zabs_orig, and
  ;; cflg as in sdss_chkciv
  if n_params() ne 1 then begin
     print,'Syntax - sdss_rmcivduplicates( civstrct_fil, [count=, index=,'
     print,'                               complement=, ncomplement=, '
     print,'                               icomplement=, check=])'
  endif 

  if size(civstrct_fil,/type) eq 8 then civstrct = civstrct_fil $
  else civstrct = xmrdfits(civstrct_fil,1,/silent)

  civstrct = sdss_srtcivstrct(civstrct)
  nciv = (size(civstrct,/dim))[0] > 1 ; defeats singularity prob
  mask = intarr(nciv)                 ; 0 = keep; 1 = duplicate
  check = mask - 1L                    ; to hold indices of matches

  unq = uniq(civstrct.qso_name) ; defined by the last per block
  nqso = (size(unq,/dim))[0] > 1
  nciv_per_qso = [unq[0] + 1,(shift(unq,-1)-unq)[0:nqso-2]]
  
  for qq=0L,nqso-1 do begin
     if qq gt 0 then rng = unq[qq-1] + lindgen(nciv_per_qso[qq]) + 1 $
     else rng = lindgen(nciv_per_qso[qq])
     ;; printcol,rng,'  '+civstrct[rng].qso_name,civstrct[rng].zabs_orig[0],civstrct[rng].cflg 

     ;; redshifts are already sorted
     zunq = uniq(civstrct[rng].zabs_orig[0])
     nzunq = (size(zunq,/dim))[0] > 1 ; defeats singularity problem
     
     if nzunq ne nciv_per_qso[qq] then begin
        ;; Then there's duplicate to trace down so want to
        ;; modify rng and instantiate complement
        if nzunq eq 1 then nciv_per_z = [zunq[0]+1] $
        else $
           nciv_per_z = [zunq[0] + 1,(shift(zunq,-1)-zunq)[0:nzunq-2]]
        zunq = rng[zunq]
        bd = where(nciv_per_z ne 1,nbd)
        if nbd ne 0 then begin
           bd2 = where(nciv_per_z[bd] eq 2,nbd2,$ ; easier to handle
                      complement=bd3,ncomplement=nbd3)
           if nbd2 ne 0 then begin
              ;; Since just 2 then can do quick global where() grep
              ;; without looping; uniq() returns the *last* incidence
              ;; of a sorted list
              ;; qso_name and zabs_orig[0] already matched
              bd2 = zunq[bd[bd2]]
              subbd2 = where(civstrct[bd2].cflg eq $
                             civstrct[bd2-1].cflg,nsubbd2,$
                             complement=subgd,ncomplement=nsubgd) 
              if nsubbd2 ne 0 then begin
                 mask[bd2[subbd2]]++ ; duplicate

                 ;; Make sure to take most favorable rating for one
                 ;; being kept
;                 stop,'sdss_rmcivduplicates(): check this subbd2 is right'
                 ;; printcol,bd2,'  '+civstrct[bd2].qso_name,civstrct[bd2].zabs_orig[0],civstrct[bd2].cflg        
                 ;; printcol,bd2-1,'  '+civstrct[bd2-1].qso_name,civstrct[bd2-1].zabs_orig[0],civstrct[bd2-1].cflg
                 ;; print,bd2[subbd2],'  '+civstrct[bd2[subbd2]].qso_name,civstrct[bd2[subbd2]].zabs_orig[0],civstrct[bd2[subbd2]].cflg
                 ;; print,bd2[subbd2]-1,'  '+civstrct[bd2[subbd2]-1].qso_name,civstrct[bd2[subbd2]-1].zabs_orig[0],civstrct[bd2[subbd2]-1].cflg


                 test = where(civstrct[bd2[subbd2]-1].rating[0] ne $
                              civstrct[bd2[subbd2]].rating[0])
                 if test[0] ne -1 then $
                       printcol,'sdss_rmcivduplicates(): adopting higher rating[0] ',$
                             civstrct[bd2[subbd2[test]]-1].qso_name,$
                             civstrct[bd2[subbd2[test]]-1].zabs_orig[0],$
                             civstrct[bd2[subbd2[test]]-1].rating[0],$
                             civstrct[bd2[subbd2[test]]].rating[0]
                 civstrct[bd2[subbd2]-1].rating[0] = $
                    civstrct[bd2[subbd2]-1].rating[0] > $
                    civstrct[bd2[subbd2]].rating[0]

                 check[bd2[subbd2]-1] = bd2[subbd2]
              endif             ; nsubbd2 ne 0
           endif                 ; nbd2 ne 0

           bd3 = where(nciv_per_z[bd] eq 3,nbd3)
           if nbd3 ne 0 then begin
              print,'sdss_rmcivduplicates(): attempting to handle triplets',nciv_per_z[bd[bd3]]
              ;; Have to loop (slow!)
              bd3 = zunq[bd[bd3]]
;              stop,'sdss_rmcivduplicates(): do not know how to handle
;              triplets+'
              ;; Just assume it's a triplet
              subbd3 = where(civstrct[bd3].cflg eq civstrct[bd3-1].cflg $
                             and civstrct[bd3].cflg eq civstrct[bd3-2].cflg,$
                             nsubbd3,complement=subgd2,ncomplement=nsubgd2)

              if nsubbd3 ne 0 then begin
                 mask[bd3[subbd3]]++   ; duplicate
                 mask[bd3[subbd3]-2]++ ; duplicate

                 if civstrct[bd3[subbd3]-1].rating[0] ne $
                    civstrct[bd3[subbd3]].rating[0] or $
                    civstrct[bd3[subbd3]-1].rating[0] ne $
                    civstrct[bd3[subbd3]-2].rating[0] then $
                       print,'sdss_rmcivduplicates(): adopting higher rating[0] ',$
                             civstrct[bd3[subbd3]-1].qso_name,$
                             civstrct[bd3[subbd3]-1].zabs_orig[0],$
                             civstrct[bd3[subbd3]-1].rating[0],$
                             civstrct[bd2[subbd3]].rating[0],$
                             civstrct[bd2[subbd3]-2].rating[0]
                 civstrct[bd3[subbd3]-1].rating[0] = $
                    civstrct[bd3[subbd3]].rating[0] > $
                    civstrct[bd3[subbd3]-1].rating[0] > $
                    civstrct[bd3[subbd3]-2].rating[0]
                 civstrct[bd3[subbd3]-2].rating[0] = civstrct[bd3[subbd3]-1].rating[0] 

                 check[bd3[subbd3]] = bd3[subbd3] - 2
                 check[bd3[subbd3]-1] = bd3[subbd3]
              endif             ; nsubbd3 ne 0
           endif                ; nbd3 ne 0
        endif                   ; nbd ne 0

     endif                      ; nzunq ne nciv_per_qso[qq]

  endfor                        ; loop qq=nqso

  index = where(mask eq 0,count,complement=icomplement,ncomplement=ncomplement)
  subcivstr = civstrct[index]
  if ncomplement eq 0 then complement = -1 $
  else complement = civstrct[icomplement]
  
  return, subcivstr
end                             ; sdss_rmcivduplicates()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getrateddblt, strct_fil, count=count, _extra=extra
  ;; Return the subsample with the desired rating[0]
  if n_params() ne 1 then begin
     print,'Syntax - sdss_getrateddblt(strct_fil, [count=, _extra=])'
     return,-1
  endif 

  if size(strct_fil,/type) eq 8 then strct = strct_fil $
  else strct = xmrdfits(strct_fil,1,/silent)

  gd = where(strct.rating[0] eq sdss_getrating(_extra=extra),count)
  if count eq 0 then return, -1 $
  else return, strct[gd]

end                             ; sdss_getrateddblt()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getlimflg, upper=upper, lower=lower, help=help
  ;; Report upper (4) or lower (2) limit flags
  if keyword_set(help) then begin
     print,'Syntax - sdss_getlimflg([/upper,/lower,/help])'
     return, -1
  endif 
  if keyword_set(upper) then return,4
  if keyword_set(lower) then return,2
  if not keyword_set(upper) and not keyword_set(lower) then return,1 ; good!
end                             ; sdss_getlimflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_setlimflg, inflg, upper=upper, lower=lower
  ;; Binary combine upper (4) and/or lower (2) limit flags to inflg
  ;; For example, a good (flg=1) line that is a lower limit (2) has
  ;; binary combination 3. 
  if n_params() lt 1 then begin
     print,'Syntax - sdss_setlimflg(inflg,[/upper,/lower])'
     return, -1
  endif 
  outflg = inflg
  if keyword_set(upper) then outflg = outflg or sdss_getlimflg(/upper)
  if keyword_set(lower) then outflg = outflg or sdss_getlimflg(/lower)

  return, outflg
end                             ; sdss_setlimflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getewflg, boxcar=boxcar, custom=custom, $
                        orig=orig, help=help
  ;; Report EW flags
  ;; 1: orig; 8: tau; 16: boxcar; 32: gauss; 64: custom
  if keyword_set(help) then begin
     print,'Syntax - sdss_getewflg([/boxcar,/custom,/orig,/help])'
     return, -1
  endif 
  if keyword_set(orig) then return,1
  if keyword_set(boxcar) then return,16
  if keyword_set(custom) then return,64
end                             ; sdss_getewflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_setewflg, inflg, boxcar=boxcar, custom=custom, $
                        orig=orig, bad=bad
  ;; Binary combine EW flags to inflg
  ;; For example, a orig (flg=1) line that is a lower limit (2) has
  ;; binary combination 3. 
  if n_params() lt 1 then begin
     print,'Syntax - sdss_setewflg(inflg,[/boxcar,/custom,/orig,/bad])'
     return, -1
  endif 
  outflg = inflg
  if keyword_set(boxcar) then outflg = outflg or sdss_getewflg(/boxcar)
  if keyword_set(custom) then outflg = outflg or sdss_getewflg(/custom)
  if keyword_set(orig) then outflg = outflg or sdss_getewflg(/orig)
  if keyword_set(bad) then outflg = 0 ; don't do anything!

  return, outflg
end                             ; sdss_setewflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getblendflg, Iwv=Iwv, IIwv=IIwv, both=both, self=self, $
                           unblend=unblend, help=help
  ;; Report profile blending flag (and IDL stinks to say wvI, wvII
  ;; would be ambiguous so I have to do this backward stuff)
  ;; 1 = generic; 2: Iwv (e.g. 1548); 4: IIwv (e.g. 1550); both: 6;
  ;; self: 8; unblend: 0
  if keyword_set(help) then begin
     print,'Syntax - sdss_getblendflg([/Iwv,/IIwv,/both,/self,/help])'
     return, -1
  endif 

  outflg = 1                    ; default
  if keyword_set(Iwv) then outflg = 2
  if keyword_set(IIwv) then outflg = 4
  if keyword_set(both) then outflg = 6 ; binary combine of 2 and 4
  if keyword_set(self) then outflg = 8
  if keyword_set(unblend) then outflg = 0 

  return,outflg
end                             ; sdss_getblendflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_setblendflg, inflg, Iwv=Iwv, IIwv=IIwv, both=both, self=self, $
                           unblend=unblend
  ;; Binary combination of flags
  if n_params() ne 1 then begin
     print,'Syntax - sdss_setblendflg( inflg, [/Iwv,/IIwv,/both,/self])'
     return,-1
  endif                         

  outflg = inflg or sdss_getblendflg() ; generic
  if keyword_set(Iwv) then outflg = outflg or sdss_getblendflg(/Iwv)
  if keyword_set(IIwv) then outflg = outflg or sdss_getblendflg(/IIwv)
  if keyword_set(both) then outflg = outflg or sdss_getblendflg(/both)
  if keyword_set(self) then outflg = outflg or sdss_getblendflg(/self)
  if keyword_set(unblend) then outflg = sdss_getblendflg(/unblend)

  return,outflg
end                             ; sdss_setblendflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getbalflg, civbal=civbal, mgiibal=mgiibal, eig=eig, $
                         visual=visual, help=help
  ;; Report Shen et al. (2011) catalog (1: CIV BAL, 2: MgII BAL, 3:
  ;; both) or eig-selected (16) BAL flags
  ;; Shen et al. (2011), ApJS, in press, arXiv:1006.5178
  ;; http://das.sdss.org/va/qso_properties_dr7/dr7.htm
  if keyword_set(help) then begin
     print,'Syntax - sdss_getbalflg([/civbal or /mgiibal or /eig or /visual,/help]) [default: 4]'
     return, -1
  endif 
  if keyword_set(civbal) then return,1
  if keyword_set(mgiibal) then return,2
  if keyword_set(visual) then return,8
  if keyword_set(eig) then return,16
  return,4                      ; unknown
end                             ; sdss_getbalflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_fndbalqso, sdsssum, count=count, balflg=balflg, full=full, $
                         flgonly=flgonly
  ;; Match BAL QSOs from Shen et al (2011) to input QSO structure by
  ;; returning the where() and optionally, the full list of flags
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fndbalqso(sdsssum, [count=, balflg=, /full, /flgonly]'
     return, -1
  endif 
  sdssdir = sdss_getsdssdir()

  if keyword_set(full) then begin
     ;; Restore save file (ibal and balflg); faster!
     restore,getenv('XIDL_DIR')+'/SDSS/CIV/sdss_balqso_map.sav'
     count = (size(ibal,/dim))[0]
     if keyword_set(flgonly) then return, balflg $
     else return,ibal
  endif 

  ;; Otherwise must search
  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]
  balflg = replicate(0,nsdss)
  
  ;; Load BAL data
  baltab = sdss_getqsostrct(/bal)
  nbal = (size(baltab,/dim))[0]
  
  ;; Compare with minimal looping
  nfil = nbal < nsdss
  for qq=0L,nfil-1 do begin
     ;; Search two ways
     if nbal lt nsdss then begin
        ibal = qq
        mtch = where(baltab[qq].mjd eq sdsstab.smjd and $
                     baltab[qq].plate eq sdsstab.plate and $
                     baltab[qq].fiber eq sdsstab.fiber,nmtch)
        if nmtch gt 1 then $
           stop,'sdss_instantbalflg: multiple matches!'
        if nmtch ne 0 then mtch = mtch[0]
     endif else begin
        mtch = where(baltab.mjd eq sdsstab[qq].smjd and $
                     baltab.plate eq sdsstab[qq].plate and $
                     baltab.fiber eq sdsstab[qq].fiber,nmtch)
        if nmtch gt 1 then $
           stop,'sdss_instantbalflg: multiple matches!'
        if nmtch ne 0 then begin
           ibal = mtch[0]
           mtch = qq            ; want sdsstab reference
        endif 
     endelse 

     ;; If no match then move on
     if nmtch eq 0 then continue
     balflg[mtch] = baltab[ibal].bal_flag
  endfor                        ; loop qq=nfil

  if keyword_set(flgonly) then return, balflg $
  else return, where(balflg gt 0,count)
end                             ; sdss_fndbalqso()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_instantbalflg, sdsssum, debug=debug
  ;; Instantiate sdsscontistrct.balflg from Shen et al. (2011) BAL QSO
  ;; catalog 
  if n_params() ne 1 then begin
     print,'Syntax - sdss_instantbalflg, sdsssum, [/debug]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]
  
  ibal = sdss_fndbalqso(sdsstab,count=nbal,balflg=balflg)
  if nbal eq 0 then return

  for bb=0L,nbal-1 do begin 
     cfil = sdss_getname(sdsstab[ibal[bb]],/abslin,dir=cdir)
     cfil = cfil[0]
     cstrct = xmrdfits(sdssdir+cdir+cfil,1,/silent)
     cstrct.balflg = balflg[ibal[bb]]
     mwrfits,cstrct,sdssdir+cdir+cfil,/create,/silent
     spawn,'gz -f '+sdssdir+cdir+cfil
  endfor                        ; loop bb=nbal

  print,'sdss_instantbalflg: Number of BALs flagged',nbal

end                             ; sdss_instantbalflg


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getspecwave
  ;; Retun SDSS wavelength range
  return, [3820.d,9200.d]
end                             ; sdss_getspecwave


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getskylinwave,dpix=dpix,dwv=dwv
  ;; Retun strong skyline wavelengths that sdss_fndlin and sdss_fndciv
  ;; work around
  dpix = [5,2]                  ; +/-5 pix, +/-2 pix in sdss_fndlin
  dwv = [5d,5d]                ; +/-5 Angstrom in sdss_fndciv
  return, [5579d,6302d]
end                             ; sdss_getskylinwave


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_setcosmology, cosmology=cosmology, songaila=songaila, $
                            dodorico=dodorico, c737=c737, age=age, help=help
  if keyword_set(help) then begin
     print,'Sytax - sdss_setcosmology([cosmology=, /songaila, /dodorico, /c737, /help])'
     return,-1
  endif 
  
  if not keyword_set(cosmology) then begin
     if keyword_set(songaila) then begin
        cosmology = [65., 1.0, 0.0]  ; q0=0.5 (Lambda=0)
        age = 10.028                 ; Gyr
     endif else begin
        if keyword_set(dodorico) then begin
           cosmology = [72., 0.26, 0.74]  ; unclear origin
           age = 13.616                   ; Gyr
        endif else begin
           if keyword_set(c737) then begin
              cosmology = [70., 0.3, 0.7]  ; common cosmology
              age = 13.462                 ; Gyr
           endif else begin
              cosmology = [71.9, 0.258, 0.742] ; WMAP5
              age = 13.664                     ; Gyr
           endelse                             ; = WMAP5                           
        endelse                                ; != D'Odorico et al. (2010)
     endelse                                   ; != Songaila (2001)
  endif                                        ; cosmology not set

  cosm_common, H0=cosmology[0], Omegavac=cosmology[2], $
               OmegaDM=cosmology[1],/silent 

  ;; Cosmology can be input or returned default
  return, cosmology

end                             ; sdss_setcosmology()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcdxdz, z, cosmology
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcdxdz(z,cosmology)'
     return,-1
  endif

  dxdz = (1.+z)^2 / sqrt(cosmology[2] + cosmology[1]*(1.+z)^3)
  return,dxdz

end


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getspecpixscale, loglam=loglam
  ;; Retun SDSS pixel scale
  val = alog(10) * 0.0001 
  if keyword_set(loglam) then return, val $ ; dwv = val * wave
  else return, val * 299792.458 ; 69 km/s
end                             ; sdss_getspecpixscale


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_getrandqsosmpl, nrand, ostrct_fil, list_fil=list_fil, $
                         zlim=zlim, rmaglim=rmaglim, bal=bal, _extra=extra
  ;; Draw random sampe of QSOs with distribution of whole
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getrandqsosmpl, nrand, ostrct_fil, [list_fil=,'
     print,'                  /bal, zlim=,rmaglim=,_extra=]'
     return
  endif 
  
  sdsstab = sdss_getqsostrct()
  if keyword_set(bal) then begin
     ibal = sdss_fndbalqso(sdsstab,count=nsdss,/full) ; full is short cut
     sdsstab = sdsstab[ibal]
  endif else nsdss = (size(sdsstab,/dim))[0]

  ;; Subsets
  if keyword_set(zlim) then begin
     gd = where(sdsstab.z ge zlim[0] and sdsstab.z le zlim[1],nsdss)
     if nsdss eq 0 then begin
        print,'sdss_getrandqsosmpl: no QSOs in z cut'
        return
     endif 
     sdsstab = sdsstab[gd]
  endif 
  if keyword_set(rmaglim) then begin
     gd = where(sdsstab.psf_r ge rmaglim[0] and sdsstab.psf_r le rmaglim[1],nsdss)
     if nsdss eq 0 then begin
        print,'sdss_getrandqsosmpl: no QSOs in R-magnitude cut'
        return
     endif 
     sdsstab = sdsstab[gd]
  endif 


  ;; Get a random sample
  ;; See http://www.idlcoyote.com/code_tips/randperm.html
  x = lindgen(nsdss)
  y = randomu(dseed, nsdss)     ; uniformly scrambled #s
  z = x[sort(y)]
  if nrand lt nsdss then $
     substrct = sdsstab[z[0:nrand-1]] $; take first bit 
  else begin
     stop,'sdss_getrandqsosmpl: With cuts, not enough QSOs ',nsdss
     substrct = sdsstab
  endelse 

  ;; Sort by plate name
  spec_fil = sdss_getname(substrct,plate=plate)
  srt = sort(plate)
  plate = plate[srt]
  spec_fil = spec_fil[srt]
  substrct = substrct[srt]

  if keyword_set(list_fil) then $
     sdss_wrqsolist,substrct,list_fil,_extra=extra ; absdir=, /silent

  if size(ostrct_fil,/type) eq 7 then begin
     mwrfits,substrct,ostrct_fil,/create,/silent
     print,'sdss_getrandqsosmpl: created ',ostrct_fil
  endif else ostrct_fil = substrct
  
end                             ; sdss_getrandqsosmpl



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcnormerr, flux, error, conti_fil, cflg=cflg, unnorm=unnorm, $
                           baderrval=baderrval
  ;; Compute the total error for spectrum, including continuum and
  ;; flux errors.
  ;; Can parse conti_fil as file name, sdsscontistrct format, or as
  ;; normal [npix,3] array
  if n_params() ne 3 then begin
     print,'Syntax - sdss_calcnormerr(flux,error,conti_fil,[cflg=,/unnorm,baderrval=])'
     return,-1
  endif 
  sdssdir = sdss_getsdssdir()
  if keyword_set(unnorm) then power = 1. else power = 2.
  if not keyword_set(baderrval) then baderrval = 0. ; for badpix
  npix = (size(flux,/dim))[0]
  bdpix = where(error eq 0.)

  ;; Read in many options
  case size(conti_fil,/type) of 
     7: begin                   ; file name
        conti = xmrdfits(sdssdir+conti_fil,0,/silent)
        if size(conti,/n_dimen) ne 2 then begin
           cstrct = xmrdfits(sdssdir+conti_fil,1,/silent) ; next extension
           if size(cstrct,/type) ne 8 then $
              stop,'sdss_calcnormerr(): non-sensical input file ',conti_fil
           if not keyword_set(cflg) then $
              stop,'sdss_calcnormerr(): must set cflg'
           ;; Use structure
           conti = dblarr(cstrct.npix,3,/nozero)
           cindx = fix(alog(cflg)/alog(2))
           conti[*,0] = cstrct.conti[0:cstrct.npix-1,cindx]
           conti[*,2] = cstrct.sigconti[0:cstrct.npix-1,cindx]
           ipix0 = cstrct.ipix0
        endif 
     end 
     8: begin                   ; structure
        cstrct = conti_fil
        if keyword_set(cflg) then cstrct.cflg = cflg 
        if not keyword_set(cstrct.cflg) then $
           stop,'sdss_calcnormerr(): must set abslin cflg or cflg keyword'
        ;; Use structure
        conti = dblarr(cstrct.npix,3,/nozero)
        cindx = fix(alog(cstrct.cflg)/alog(2))
        conti[*,0] = cstrct.conti[0:cstrct.npix-1,cindx]
        conti[*,2] = cstrct.sigconti[0:cstrct.npix-1,cindx]
        ipix0 = cstrct.ipix0
     end
     else: conti = conti_fil    ; better be right already
  endcase
  if not keyword_set(ipix0) then begin
     ;; Figure out ipix0
     bd = where(conti[*,0] eq 0. or finite(conti[*,0],/nan))
     if bd[0] eq -1 then ipix0 = 0 $
     else begin
        istrt = where(bd ne shift(bd,1)+1,ngap)
        istop = where(bd ne shift(bd,-1)-1,ntest)
        ipix0 = istop[0]        ; just the beginning
     endelse 
  endif 

  if npix ne (size(conti,/dim))[0] then $
     stop,'sdss_calcnormerr(): input flux and continuum not same size'
  normerr = error * 0.
  gd = where(conti[*,0] ne 0. and finite(conti[*,0]),complement=bd)
  normerr[gd] = sqrt(error[gd]^2*conti[gd,0]^2 + $
                     conti[gd,2]^2*flux[gd]^2)/conti[gd,0]^power
  if ipix0 ne 0 then begin
     normerr[0:ipix0] = error[0:ipix0] 
     if not keyword_set(unnorm) then begin
        gd = where(conti[0:ipix0,0] ne 0. and finite(conti[0:ipix0]),complement=bd)
        if gd[0] ne -1 then normerr[gd] = normerr[gd]/conti[gd,0]
        if bd[0] ne -1 then normerr[bd] = baderrval
     endif 
  endif 
  if bdpix[0] ne -1 then normerr[bdpix] = baderrval

  return, normerr
end                             ;  sdss_calcnormerr()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_normspec, spec_fil, conti_fil, out_fil, cflg=cflg, clobber=clobber, $
                   stack=stack
  if n_params() ne 3 then begin
     print,'Syntax - sdss_normspec, spec_fil, conti_fil, out_fil, '
     print,'                        [cflg=, /clobber, /stack]'
     return
  endif
  
  ;; parse_sdss handles different SDSS formats
  parse_sdss, spec_fil, flux, wave, sig=error, head=head, npix=npix

  if keyword_set(stack) then begin
     ;; Assume spec_fil is output of sdss_stackciv.pro so that the
     ;; continuum structure is in 1st extension of spec_fil; this just
     ;; ignores conti_fil, whatever it's set to 
     cstrct = xmrdfits(spec_fil,1,/silent)
     if not keyword_set(cflg) then cflg = cstrct.cflg
  endif else begin
     ;; Figure out the type of continuum: structure from Precious Metals
     ;; survey (eigconti, hybconti, splconti) or x_continuum
     if size(conti_fil,/type) eq 7 then begin
        cstrct = xmrdfits(conti_fil,0,/silent) 
        if not (size(cstrct,/n_dimen) ge 2 or $
                n_elements(cstrct) eq npix) then begin ; not x_continuum output
           cstrct = xmrdfits(conti_fil,1,/silent)      ; next extension
           if not keyword_set(cflg) then cflg = cstrct.cflg
        endif
     endif else cstrct = conti_fil ; structure or multi-dim array
  endelse 

  if size(cstrct,/type) eq 8 then begin
     if not keyword_set(cflg) then $
        stop,'sdss_normspec: must set cflg'
     ;; Use structure
     cindx = fix(alog(cflg)/alog(2))
     conti0 = cstrct ; for sdss_calcnormerr()
     conti = cstrct.conti[0:cstrct.npix-1,cindx]
     sigconti = cstrct.sigconti[0:cstrct.npix-1,cindx]
     ipix0 = cstrct.ipix0
  endif else begin
     conti = cstrct
     if size(conti,/n_dim) gt 1 then begin
        conti0 = conti ; for sdss_calcnormerr()
        sigconti = conti[*,2]
        conti = conti[*,0]
        gd = where(conti gt 0,npix_conti)
        ipix0 = gd[0]           ; or zero?
     endif else begin
        ipix0 = 0
        sigconti = 0
     endelse 
  endelse

  spec_norm = dblarr(npix,5)  ; output will be traditional SDSS

  gdpix = where(conti[ipix0:*] gt 0.,complement=bdpix)
  
  if gdpix[0] ne -1 then begin
     spec_norm[gdpix,0] = flux[gdpix]/conti[gdpix]
     if keyword_set(sigconti) then $
        spec_norm[*,2] = sdss_calcnormerr(flux,error,conti0,_extra=extra) $
     else $
        spec_norm[gdpix,2] = error[gdpix]/conti[gdpix]
  endif

  ;; Write out 
  test = file_search(out_fil+'*',count=ntest)
  if ntest eq 0 or keyword_set(clobber) then begin
     mwrfits,spec_norm,out_fil,head,/create,/silent
     spawn,'gzip -f '+out_fil
     print,'sdss_normspec: created ',out_fil
  endif else $
     stop,'sdss_normspec stop: will not clobber file ',out_fil

end ; sdss_normspec


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltconti, spec_fil, all=all, norm=norm, zabs=zabs, $
                   stack=stack, _extra=extra
  ;; Read in and plot the spectrum with continuum
  ;; if /all set, still have to put /eig or /spl
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltconti, spec_fil, [/all, /norm, zabs=, /stack, _extra=]'
     return
  end 
  sdssdir = sdss_getsdssdir()
  if keyword_set(all) and keyword_set(norm) then begin
     print,'sdss_pltconti: cannot display all spline and eigen-normalized spectra.'
     return
  endif 

  ;; Find file
  if not keyword_set(stack) then begin
     sfil = sdss_getname(spec_fil,/spec,dir=dir)
     test = file_search(sdssdir+dir+sfil+'*',count=ntest)
     if ntest eq 0 then $
        stop,'sdss_pltconti: spec_fil DNE as given or in full path'

     ;; Find continuum (_extra includes /spl or /eig or /abslin)
     cfil = sdss_getname(spec_fil,/spec,dir=dir,_extra=extra)
     
     ;; Find continuum file
     test2 = file_search(sdssdir+dir+cfil+'*',count=ntest2)
     if ntest2 eq 0 then $
        stop,'sdss_pltconti: continuum file DNE ',dir+cfil

     ;; Read in from explicitly stored continuum
     ;; Test kind of file
     cstrct = xmrdfits(test2,0,/silent) ; array or header
  endif else begin
     test = spec_fil 
     prs = strsplit(spec_fil,'/',/extract,count=nprs)
     cfil = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'.',/reverse_search))
     cstrct = xmrdfits(spec_fil,1,/silent)
  endelse

  ;; Read in files (use test b/c has full path and any .gz)
  parse_sdss,test,fx,wv,sig=er0,zqso=zqso ; zqso may be bad!
  
  if size(cstrct,/n_dimen) eq 2 then begin
     conti = cstrct[*,0]
     ;; er = sqrt(er0^2*cstrct[*,0]^2 + cstrct[*,2]^2*fx^2)/cstrct[*,0]
     er = sdss_calcnormerr(fx,er0,cstrct,/unnorm)

     if keyword_set(all) then begin
          cfil2 = sdss_getname(spec_fil,/spec,/spl)
          if cfil2 eq cfil then $
             cfil2 = sdss_getname(spec_fil,/spec,/eig) ;  try other
          tmp = xmrdfits(sdssdir+dir+cfil2,0,/silent)
          altconti = tmp[*,0]
     endif 
     print,'sdss_pltconti: Using spectrum header redshift, which may be wrong'
  endif else begin
     ;; Figure out which version
     if size(cstrct,/type) ne 8 then $
        cstrct = xmrdfits(test,1,/silent) ; structure
     zqso = cstrct.z_qso               ; better than the header!
     cindx = fix(alog(cstrct.cflg)/alog(2))
     conti = cstrct.conti[0:cstrct.npix-1,cindx]
     ;; er = sqrt(er0^2*cstrct.conti[0:cstrct.npix-1,cindx]^2 + $
     ;;           cstrct.sigconti[0:cstrct.npix-1,cindx]^2*fx^2)/$
     ;;           cstrct.conti[0:cstrct.npix-1,cindx]
     er = sdss_calcnormerr(fx,er0,cstrct,/unnorm,cflg=cstrct.cflg)
     if cstrct.cflg eq sdss_getcflg(/eig) or $
        cstrct.cflg eq sdss_getcflg(/hyb) then begin
        altconti = cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)]
     endif else begin
        conti = cstrct.conti[0:cstrct.npix-1,cindx]
        altconti = cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/eig,/index)]
     endelse 
  endelse
  if keyword_set(all) then begin
     psym3 = -3 
     print,'sdss_pltconti: blue is other conti'
  endif else begin 
     altconti = er0
     psym3 = 10
     print,'sdss_pltconti: blue is original error'
  endelse 

  if keyword_set(norm) then begin
     fx = fx/conti
     er = er/conti
     conti = 0
  endif 

  ;; Plot; _extra= include xrange=, etc
  if keyword_set(zabs) then zin = zabs else zin = zqso
  x_specplot,fx,er,wav=wv,inflg=4,ytwo=conti,psym2=-3,$
             ythr=altconti,psym3=psym3,zin=zin,/lls,$
             title=cfil+' zqso = '+string(zqso,format='(f8.5)'),$
             _extra=extra, /block
end                             ; sdss_pltconti



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltobj, input, qsostrct=qsostrct, radec=radec, jname=jname, $
                 _extra=extra
  ;; Find and plot the object's spectrum, assuming standard
  ;; directory structure. Default is to search by obj id
  ;; (MJD-plate-fiber) 
  sdssdir = sdss_getsdssdir()

  ;; _extra= inccludes: dr=, name=, /BAL, /noBAL
  if not keyword_set(qsostrct) then $
     qsostrct = sdss_getqsostrct(_extra=extra)

  if keyword_set(radec) then begin
     if size(input,/type) eq 7 then begin
        ;; String RA, Dec
        x_radec,radec[0],radec[1],ra,dec
     endif else begin
        ;; in degrees
        ra = input[0]
        dec = input[1]
     endelse

     gcirc,2,ra,dec,qsostrct.ra,qsostrct.dec,das
     mn = min(das,iqso)
     print,'sdss_pltobj: matched on RA, Dec to '+$
           string(das,format='(f6.2)')+' arcsec'
  endif else begin
     
     if keyword_set(jname) then $
        ;; Matching by coordinate in name
        iqso = where(stregex(input,'J'+qsostrct.sdssj,/boolean),nmtch) $
     else begin
        ;; Matching by MJD-plate-fiber
        spec_fil = sdss_getname(qsostrct,root=qso_name)
        iqso = where(qso_name eq input,nmtch)
     endelse
     
     case nmtch of
        0: begin
           print,'sdss_pltobj: no matches for ',input 
           return
        end
        1: iqso = iqso[0]
        else: begin
           print,'sdss_pltobj: multiple matches for ',input
           if keyword_set(jname) then $
              printcol,indgen(nmtch),qsostrct[iqso].sdssj,format='(i3,1x,a)' $
           else $
              printcol,indgen(nmtch),qso_name[iqso],format='(i3,1x,a)'
           print,'  Which number or q(uit):'
           done = 0
           while not done do begin
              ch = get_kbrd(1)  ; w/(1) wait for character
              case ch of
                 113: begin     ; q
                    print,'sdss_pltobj: exiting'
                    done = 1
                    return
                 end
                 else: begin
                    gd = where(ch eq indgen(ntmch),ngd)
                    if ngd eq 0 then begin
                       print,'  Invalid entry: which number or q(uit):'
                    endif else begin
                       iqso = iqso[ch] ; set number
                       done = 1
                    endelse
                 end
              endcase
           endwhile             ; wait for selection
        end
     endcase

  endelse                       ; jname = 0

     ;; ploti 

end ; sdss_pltobj

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltsnr, spec_fil, lsnr=lsnr, cfil2=cfil2, debug=debug,$
                 _extra=extra
  ;; Plot convolved S/N used to find absorption lines
  ;; _extra= passed to x_specplot and includes zin=, /lls, etc.
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltsnr, spec_fil, [lsnr=, _extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  if not keyword_set(lsnr) then lsnr = 3.5 ; to match sdss_fndlin
  
  cfil = sdss_getname(spec_fil,/spec,/abslin,dir=cdir)    
  
  cstrct = xmrdfits(sdssdir+cdir[0]+cfil[0],1,/silent)        
  
  parse_sdss,sdssdir+cstrct.sdss_obs[0],fx,wv
  
  cindx = fix(alog(cstrct.cflg)/alog(2))
  gd = where(finite(cstrct.snr_conv[0:cstrct.npix-1,cindx]),ngd)      
  title = strtrim(spec_fil,2)+$
          ' zqso = '+string(cstrct.z_qso,format='(f8.5)')

  if keyword_set(cfil2) then begin
     if size(cfil2,/type) eq 8 then cstrct2 = cfil2 $
     else cstrct2 = xmrdfits(sdssdir+cfil2,1,/silent)
     cindx2 = fix(alog(cstrct2.cflg)/alog(2))
     gd2 = where(finite(cstrct2.snr_conv[0:cstrct2.npix-1,cindx2]),ngd)
     x_specplot,cstrct.snr_conv[gd,cindx],replicate(lsnr,ngd),wave=wv[gd],inflg=4,$
                /block,title=title,_extra=extra,$
                two_wave=wv[gd2], ytwo=cstrct2.snr_conv[gd2,cindx2]
  endif else $
     x_specplot,cstrct.snr_conv[gd,cindx],replicate(lsnr,ngd),wave=wv[gd],inflg=4,$
                /block,title=title,_extra=extra

  if keyword_set(debug) then $
     stop,'sdss_pltsnr debug: stopping before exiting'
end                             ; sdss_pltsnr


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltabslin, spec_fil, abslin_fil, cflg=cflg, noplt=noplt, $
                    fx_lin=fx_lin, wv_lin=wv_lin, wave=wave, $
                    sigma=sigma, _extra=extra
  ;; Plot the spectrum with the sdss_chkciv indication of lines
  ;; _Extra includes functions to x_specplot
  if n_params() ne 2 then begin
     print,'Syntax - sdss_pltabslin, spec_fil, abslin_fil, [cflg=, /noplt=,'
     print,'                         fx_lin=, wv_lin=, wave=, sigma=, _extra=]'
     return
  endif 

  if size(spec_fil,/type) eq 7 then $
     parse_sdss, spec_fil, flux, wave, sig=sigma, npix=npix $
  else begin
     flux = spec_fil
     if not keyword_set(wave) then $
        stop,'sdss_pltabslin: must set wave if input spec_fil is flux'
     npix = (size(flux,/dim))[0]
     if npix ne (size(wave,/dim))[0] then $
        stop,'sdss_pltabslin: size of flux and wave arrays do not match'
  endelse 

  if size(abslin_fil,/type) eq 8 then cstrct = abslin_fil $
  else cstrct = xmrdfits(abslin_fil,1,/silent)

  if not keyword_set(cflg) then cflg = cstrct.cflg
  cindx = fix(alog(cflg)/alog(2))

  if cstrct.ncent[cindx] gt 0. then begin
     fx_lin= fltarr(2*cstrct.ncent[cindx],/nozero)
     wv_lin = fx_lin
     tmp = cstrct.centroid[$
           sort(cstrct.centroid[0:cstrct.ncent[cindx]-1,$
                                cindx]),cindx]

     ncentarr = lindgen(cstrct.ncent[cindx]) ; instantiate once
     wv_lin[2*ncentarr] = tmp
     if cstrct.ncent[cindx] gt 1 then begin
        nhalfarr = lindgen(cstrct.ncent[cindx]/2.)

        ;; Prevent x_spepclot's sorting from screwing this up
        ;; by toggling a little
        wv_lin[2*ncentarr+1] = tmp       
        
        fx_lin[4*nhalfarr] = 500. ; just something large
        wv_lin[4*nhalfarr] = $
           wv_lin[4*nhalfarr] - 5.e-4
        
        fx_lin[4*nhalfarr+1] = -500. ; just something small       
        wv_lin[4*nhalfarr+1] = $
           wv_lin[4*nhalfarr+1] - 4.e-4
        
        fx_lin[4*nhalfarr+2] = -500. 
        wv_lin[4*nhalfarr+2] = $
           wv_lin[4*nhalfarr+2] - 3.e-4
        
        fx_lin[4*nhalfarr+3] = 500.       
        wv_lin[4*nhalfarr+3] = $
           wv_lin[4*nhalfarr+3] - 2.e-4
        
        fx_lin[2*cstrct.ncent[cindx]-2+[0,1]] = $
           -1.*fx_lin[2*cstrct.ncent[cindx]-4+[0,1]]
        
     endif else begin
        ;; Hand instantiate
        fx_lin = [500.,-500.]
     endelse
  endif else begin
     ;; cstrct.ncent[cstrct.cflg] = 0
     wv_lin = [0.,0.]
     fx_lin = [0.,0.]
     print,'sdss_pltabslin: no centroids detected for ',cstrct.qso_name
  endelse 
  

  if not keyword_set(noplt) then begin
     ;;
     x_specplot, flux, sigma, ytwo=cstrct.conti[0:cstrct.npix-1,cindx],$
                 psym2=-3, wave=wave, /block, $
                 inflg=4, ythree=fx_lin, three_wave=wv_lin, $
                 title=cstrct.qso_name+$
                 string(cstrct.z_qso,format='(1x,"zqso = ",f7.5)')
  endif                         ; plot

end                             ; sdss_pltabslin


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_contiqual, spec_fil, wvlim=wvlim, _extra=extra
  ;; Measure the continuum chi^2
  if n_params() ne 1 then begin
     print,'Syntax - sdss_contiqual(spec_fil, [wvlim=,_extra=])'
     return,[-1,-1]
  end 
  sdssdir = sdss_getsdssdir()

  ;; Find file
  sfil = sdss_getname(spec_fil,/spec,dir=dir)
  test = file_search(sdssdir+dir+sfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_contiqual(): spec_fil DNE as given or in full path'

  ;; Read in files (use test b/c has full path and any .gz)
  parse_sdss,test,fx,wv,sig=er,npix=npix

  ;; Find continuum (_extra includes /spl or /eig or /abslin)
  cfil = sdss_getname(spec_fil,/spec,dir=dir,_extra=extra)

  ;; Find continuum file
  test = file_search(sdssdir+dir+cfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_contiqual(): continuum file DNE ',dir+cfil
  
  ;; Read in from explicitly stored continuum
  ;; Test kind of file
  cstrct = xmrdfits(test,0,/silent) ; array or header
  if size(cstrct,/n_dimen) eq 2 then begin
     conti = cstrct[*,0]
     ;; Test what error should be becasue eigconti and splconti stored
     ;; differently; original er = 0. indicates gap
     test = where(er lt cstrct[*,2] and er ne 0.,ntest)
     if ntest eq 0 then $
        ;; er = sqrt(er^2*cstrct[*,0]^2 + cstrct[*,2]^2*fx^2)/conti^2 $
                er = sdss_calcnormerr(fx,er,cstrct) $
     else er = cstrct[*,2]
  endif else begin
     ;; Figure out which version
     cstrct = xmrdfits(test,1,/silent) ; structure
     cindx = fix(alog(cstrct.cflg)/alog(2))
     conti = cstrct.conti[0:cstrct.npix-1,cindx]
     er = cstrct.sigconti[0:cstrct.npix-1,cindx]
  endelse
  ;; Normalize 
  fx = fx/conti
  
  ;; Trim
  if keyword_set(wvlim) then begin
     gd = where(wv ge wvlim[0] and wv le wvlim[1],npix)
     if npix eq 0 then begin
        print,'sdss_contiqual(): no spectrum in wavelength limits'
        return,[-1,-1]
     endif 
     fx = fx[gd]
     er = er[gd]
  endif 

  ;; RMS and ...
  rms = sqrt(total(fx^2)/npix)
;  sigrms = sqrt( total((fx*er)^2) / (npix * total(fx^2)) ) ; weighted?
  sigrms = sqrt(total(er^2)/npix)

  return,[rms,sigrms]
end                             ; sdss_contqual()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_measuresnr, spec_fil, list=list, wvlim_obs=wvlim_obs, $
                          dvgal=dvgal, dvqso=dvqso, dblt_name=dblt_name, $
                          zqso=zqso, sdsssum=sdsssum, no_snr=no_snr, $
                          gz_fil=gz_fil, nsig=nsig, limz=limz, clobber=clobber, $
                          silent=silent
  ;; Measure S/N for unnormalized input spectrum
  if n_params() ne 1 then begin
     print,'Syntax - sdss_measuresnr(spec_fil, [/list, wvlim_obs=, dvgal='
     print,'                         dvqso=, zqso=, sdsssum=, dblt_name=, /no_snr'
     print,'                         gz_fil=, nsig=, limz=, /clobber, /silent])'
     return,-1
  endif 

  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if not keyword_set(dvgal) then dvgal = 5000.  ; km/s
  if not keyword_set(dvqso) then dvqso = -3000. ; km/s

  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  wvspecmnx = sdss_getspecwave()
  cinv = 1./2.998e5             ; km^-1 s
  wvr_lim = [1250.,1.e5]        ; matches sdss_fndciv

  if dblt.ion eq 'CIV' or dblt.ion eq 'FkIV' then $
     wvr_lim[0] = 1310.         ; exclude OI/SiII forest

  sdssdir = sdss_getsdssdir()
  if keyword_set(list) then begin
     readcol,spec_fil,file,format='(a)',skip=1,/silent
     sv_spec_fil = spec_fil     ; will revert later
     spec_fil = file
  endif  
  nfil = (size(spec_fil,/dim))[0] > 1 ; foils singularity

  pixscale = sdss_getspecpixscale(/loglam) 

  if keyword_set(nsig) then begin
     if keyword_set(gz_fil) then begin
        if ((size(gz_fil,/dim))[0] > 1) ne nfil then $
           stop,'sdss_measuresnr(): gz_fil= size does not match spec_fil size'
     endif 
  endif else if keyword_set(gz_fil) then $
     stop,'sdss_measuresnr(): must sent nsig if setting gz_fil='

  if not keyword_set(sdsssum) and not keyword_set(zqso) then begin
     print,'sdss_measuresnr(): sdsssum and zqso not set so measure S/N for full spectrum.'
     zqso = replicate(100.,nfil) ; something large
  endif else begin
     if keyword_set(sdsssum) then begin
        if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
        else sdsstab = xmrdfits(sdsssum,1,/silent)
        if (size(sdsstab,/dim))[0] ne nfil then $
           stop,'sdss_measuresnr(): sdsstab not same size as input file list'
        ;; Assume everything in right order
        zqso = sdsstab.z
     endif
  endelse 

  ;; Set wavelength bounds for comparison
  ;; and it's ALLOWED to have wvlim_obs[*,1] < wvlim_obs[*,0]
  ;; but S/N will equal -9.99
  wvlim_obs = fltarr(nfil,2,/nozero)
  wvlim_obs[*,0] = (wvspecmnx[0] > (wvr_lim[0]*(1.+zqso))) > $
                   (dblt.wvI*(1.+dvgal*cinv))
  wvlim_obs[*,1] = wvspecmnx[1] < (dblt.wvII*(1.+zqso)*(1.+dvqso*cinv))

  if keyword_set(limz) then begin
     ;; For limited redshift range, measure <S/N>
     wvlim_obs[*,0] = wvlim_obs[*,0] > dblt.wvI*(1+limz[0])
     wvlim_obs[*,1] = wvlim_obs[*,1] < dblt.wvI*(1+limz[1])
  endif                         ; limz=
  
  ;; Exit point; very fast
  if keyword_set(no_snr) then begin
     if nfil eq 1 then $
        ;; Return normal looking arrays
        wvlim_obs = reform(wvlim_obs)
     if keyword_set(sv_spec_fil) then spec_fil = sv_spec_fil
     return,-1
  endif 


  ;; ;;;;;;;
  snr = fltarr(nfil,3,/nozero)  ; median flux, error, median(flux/error)
  for ss=0L,nfil-1 do begin 
     parse_sdss, sdssdir+spec_fil[ss], flux, wave, npix=npix, sig=error, $
                 head=hdr

     sub = where(wave ge wvlim_obs[ss,0] and wave le wvlim_obs[ss,1] and $
                 error ne 0.,nsub)
     if nsub ne 0 then begin
        snr_arr = flux[sub]/error[sub]
        snr[ss,*] = [median(flux[sub],/even),$
                     median(error[sub],/even),$
                     median(snr_arr,/even)] 

     endif else $
        snr[ss,*] = -9.99


     if keyword_set(nsig) then begin
        ;; Estimate g(z) using Danforth & Shull (2008) EW_lim
        ;; EW_lim = nsig * lambda / (R_lambda * (S/N)_lambda)

        ;; gz array will be [npix, [zabs, dz, EWlim, mask, nsig]]
        gz_los = fltarr(npix,5,/nozero) ; instantiate efficiently
        gz_los[*,0] = wave/dblt.wvI - 1.
        gz_los[*,1] = pixscale * wave ; dlambda (later dz)
        gz_los[*,2] = !values.f_infinity   
        gz_los[*,3] = 0         ; mask (0 = bad; 1 = useful)
        gz_los[*,4] = nsig      ; useful number to preserve
        
        ;; Before overwrite this:
        if nsub ne 0 then $
           gz_los[sub,3] = 1    ; useful, qualified pixels

        dispstrct = xmrdfits(sdssdir+spec_fil[ss],6,/silent)
        if size(dispstrct,/type) eq 8 then $
           dispersion = dispstrct.dispersion $
        else begin
           print,'sdss_measuresnr(): dispersion structure DNE ',spec_fil[ss]
           dispersion = replicate(1,npix)
        endelse

        ;; calculate limit for every pixel even outside of bounds
        sub = where(error ne 0 and dispersion ne 0,nsub)
        if nsub ne 0 then begin
           snr_arr = flux[sub]/error[sub]
           gz_los[sub,2] = nsig * gz_los[sub,1] / $ ; why gz_los[*,1] = dlambda
                           (dispersion[sub] * snr_arr)
           
           bd = where(dispersion[sub] eq 0,complement=gd)
           if bd[0] ne -1 then begin
              ;; Don't leave these as infinite though don't know why
              ;; would ever get here
              ;; likely never will
              gz_los[sub[bd],2] = interpol(gz_los[sub[gd],2],$
                                           gz_los[sub[gd],0],$
                                           gz_los[sub[bd],0])
              print,'sdss_measuresnr(): interpolated some pixels ',spec_fil[ss]
           endif                ; interpolate dispersion
           
           bd = where(gz_los[sub,2] le 0.)
           if bd[0] ne -1 then $
              gz_los[sub[bd],3] = 0 ; back to bad, not useful
        endif                       ; nsub ne 0

        gz_los[*,1] = gz_los[*,1]/dblt.wvI ; now dz
        
        if size(gz_fil,/type) eq 7 then begin
           ;; Write file
           test = file_search(sdssdir+gz_fil[ss]+'*',count=ntest)
           if ntest eq 0 or keyword_set(clobber) then begin
              if ntest eq 0 then begin
                 ;; Check if need to make directory
                 dir = strmid(gz_fil[ss],0,strpos(gz_fil[ss],'/',$
                                                  /reverse_search))
                 dtest = file_search(sdssdir+dir,count=ndtest)
                 if ndtest eq 0 then $
                    spawn,'mkdir -p '+dir
              endif 
              mwrfits,gz_los,sdssdir+gz_fil[ss],hdr,/create,/silent
              spawn,'gzip -f '+sdssdir+gz_fil[ss]
              if not keyword_set(silent) then $
                 print,'sdss_measuresnr(): created ',gz_fil[ss]
           endif                ; created file
           
        endif else gz_fil = gz_los ; return array
     endif                         ; keyword_set(nsig)
  
  endfor                        ; loop ss=nfil

  if nfil eq 1 then begin
     ;; Return normal looking arrays by collapsing single element dimensions
     snr = reform(snr)
     wvlim_obs = reform(wvlim_obs)
  endif 
  
  if keyword_set(sv_spec_fil) then spec_fil = sv_spec_fil
  return, snr
end                             ; sdss_measuresnr()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_mksnrtab, list_fil, sdsssum, out_fil, _extra=extra
  ;; Create structure of all S/N per doublet
  if n_params() ne 3 then begin
     print,'Syntax - sdss_mksnrtab, list_fil, sdsssum, out_fil, [_extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  dblt_arr = ['SiIV','CIV','MgII','CaII']
  ndblt = (size(dblt_arr,/dim))[0]

  readcol,list_fil,spec_fil,skip=1,format='a',/silent
  nfil = (size(spec_fil,/dim))[0]
  abs_fil = sdss_getname(spec_fil,/spec,/abslin,dir=absdir)
  abs_fil = absdir+abs_fil

  if size(sdsssum,/type) eq 7 then sdsstab = xmrdfits(sdsssum,1,/silent) $
  else sdsstab = sdsssum
  if nfil ne (size(sdsstab,/dim))[0] then $
     stop,'sdss_mksnrtab: list_fil and sdsssum unequal sizes.'

  strct_tmp = {qso_name:'', $
               z_qso:0d, $
               ncent:lonarr(3)-1, $ ; match sdsscontistrct
               SNR_SiIV:fltarr(3,/nozero), $
               wvobs_SiIV:fltarr(2,/nozero), $
               SNR_CIV:fltarr(3,/nozero), $
               wvobs_CIV:fltarr(2,/nozero), $
               SNR_MgII:fltarr(3,/nozero), $
               wvobs_MgII:fltarr(2,/nozero), $
               SNR_CaII:fltarr(3,/nozero), $
               wvobs_CaII:fltarr(2,/nozero) $
              }
  tags = tag_names(strct_tmp)

  snr_strct = replicate(strct_tmp,nfil)
  istrt = strpos(spec_fil[0],'-')
  istop = strpos(spec_fil[0],'.',/reverse_search) 
  snr_strct.qso_name = strmid(spec_fil,istrt+1,istop-istrt-1)
  snr_strct.z_qso = sdsstab.z

  for ii=0,ndblt-1 do begin
     ;; _extra includes dvgal=, dvqso= 
     snr_arr = sdss_measuresnr(spec_fil,wvlim_obs=wvlim_obs,$
                               dblt_name=dblt_arr[ii],zqso=snr_strct.z_qso, $
                              _extra=extra)
     snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
     wvtag = where(tags eq 'WVOBS_'+strupcase(dblt_arr[ii]))
     snr_strct.(snrtag)[0] = snr_arr[*,0]
     snr_strct.(snrtag)[1] = snr_arr[*,1]
     snr_strct.(snrtag)[2] = snr_arr[*,2]
     snr_strct.(wvtag)[0] = wvlim_obs[*,0]
     snr_strct.(wvtag)[1] = wvlim_obs[*,1]
  endfor                        ; for ii=ndblt

  for ii=0L,nfil-1 do begin
     test = file_search(sdssdir+abs_fil[ii]+'*',count=ntest)
     if ntest gt 0 then begin
        cstrct = xmrdfits(sdssdir+abs_fil[ii],1,/silent)
        snr_strct[ii].ncent = cstrct.ncent
     endif ; else will stay -1 in all
  endfor                        ; loop ii=nfil

  ;; Write out
  if size(out_fil,/type) eq 7 then begin
     mwrfits,snr_strct,out_fil,/create,/silent
     spawn,'gzip -f '+out_fil
     print,'sdss_mksnrtab: created ',out_fil
  endif else out_fil = snr_strct

end                             ; sdss_mksnrtab


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntsnrtab, snrstrct_fil, out_fil, snrsumm=snrsumm, zbin=zbin, $
                     dblt_name=dblt_name, zmin=zmin
  ;; Print out a variety of cuts based on S/N
  ;; if dblt_name not set, input all 4 usuals
  ;; if zmin not set, then use all (zmin = 0)
  ;; Default table is QSO ID, zqso, S/N for each dblt_name
  ;; if snrsumm = 1, then print out a by-redshift breakdown for all
  ;; dblt_names for S/N >= 3, 4, and 5
  ;; if snrsumm > 1, then print out list of spectra with S/N >= snrsumm
  if n_params() ne 2 then begin
     print,'Syntax - sdss_prntsnrtab, snrstrct_fil, out_fil, [snrsumm=, zbin=,'
     print,'                          dblt_name=, zmin=]'
     return
  endif 
  if not keyword_set(zbin) then zbin = 0.2
  if not keyword_set(zmin) then zmin = 0.0

  if size(snrstrct_fil,/type) eq 7 then $
     snr_strct = xmrdfits(snrstrct_fil,1,/silent) $
  else snr_strct = snrstrct_fil
  tags = tag_names(snr_strct[0])

  close,1
  openw,1,out_fil

  if keyword_set(snrsumm) then begin
     if keyword_set(dblt_name) then dblt_arr = dblt_name $
     else dblt_arr = ['SiIV','CIV','MgII','CaII'] 
     ndblt = n_elements(dblt_arr) ; do not use size(/dim)


     if snrsumm gt 1 then begin
        ;; Print list
        for ii=0,ndblt-1 do begin
           snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
           gd = where(snr_strct.(snrtag)[2] ge snrsumm and $
                      snr_strct.z_qso ge zmin,ngd)
           prs = strsplit(snr_strct[0].qso_name,'-')
           spec_fil = 'spectro/1d_26/'+strmid(snr_strct[gd].qso_name,prs[1],prs[2]-prs[1]-1)+$
                      '/1d/spSpec-'+snr_strct[gd].qso_name+'.fit'

           printf,1,'abslin/'
           
           writecol,out_fil,spec_fil,snr_strct[gd].(snrtag)[2],fmt='(a,5x,f7.2)',$
                    filnum=1
        endfor                  ; loop ii=ndblt

        close,1
        print,'sdss_prntsnrtab: created spectra list for S/N >= '+$
                 string(snrsumm,format='(f4.1," ")')+out_fil
        
     endif else begin
        ;; Print S/N cut summary
        for ii=0,ndblt-1 do begin

           snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
           
           sub = where(snr_strct.(snrtag)[2] gt 0. and $
                       snr_strct.z_qso ge zmin,nsub) ; has coverage
           zmn = min(snr_strct[sub].z_qso,max=zmx)
           

           gd3 = sub[where(snr_strct[sub].(snrtag)[2] ge 3.,ngd3)]
           hist3 = histogram(snr_strct[gd3].z_qso,loc=loc3,min=zmn,max=zmx,binsize=zbin)

           gd4 = sub[where(snr_strct[sub].(snrtag)[2] ge 4.,ngd4)]
           hist4 = histogram(snr_strct[gd4].z_qso,loc=loc4,min=zmn,max=zmx,binsize=zbin)

           gd5 = sub[where(snr_strct[sub].(snrtag)[2] ge 5.,ngd5)]
           hist5 = histogram(snr_strct[gd5].z_qso,loc=loc5,min=zmn,max=zmx,binsize=zbin)

           printf,1,'' 
           printf,1,dblt_arr[ii]," Total: ",string(nsub,format='(i5)') 
           printf,1,'zbin','S/N=3','4','4/3','5','5/3',$
                  format='(a6,2x,a5,2x,2(a5,1x,"(",a4,")",2x))'
           writecol,out_fil,loc3,hist3,hist4,hist4/float(hist3),hist5,hist5/float(hist3),$
                    fmt='(f6.2,2x,i5,2x,2(i5,1x,"(",f4.2,")",2x))' ,filnum=1
           printf,1,"Total",total(hist3),total(hist4),total(hist5), $
                  format='(a6,2x,i5,2x,2(i5,9x))' 
           
        endfor                  ; loop ii=ndblt
        close,1
        print,'sdss_prntsnrtab: created S/N-z summary tables ',out_fil

     endelse 

  endif else begin
     ;; Print all S/N information
     printf,1,"DR7 QSO S/N Summary (for regions where ion's detectable)"
     printf,1,'QSO ID','zqso','SiIV','CIV','MgII','CaII',$
            format='(a15,2x,a7,4(1x,a6))'

     writecol,out_fil,snr_strct.qso_name,snr_strct.z_qso,$
              snr_strct.SNR_SiIV[2],snr_strct.SNR_CIV[2],$
              snr_strct.SNR_MgII[2],snr_strct.SNR_CaII[2],$
              fmt='(a15,2x,f7.5,4(1x,f6.2))',filnum=1
     close,1
     print,'sdss_prntsnrtab: created all-S/N summary table ',out_fil
  endelse 

end                             ; sdss_prntsnrtab


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_reorganize, help=help
  ;; Move old EIGCONTI/, SPLCONTI/, and ABSLIN/ stuff to new directory structure
  ;; Be in the right directory (also tests that $SDSSDR set)
  if keyword_set(help) then begin
     print,'Syntax - sdss_reorganize [/help]'
     return
  endif 
  pwd, cur_dir
  cd,sdss_getsdssdir()

  ;; Do it for all
  sdsstab = sdss_getqsostrct() 
  nsdss = (size(sdsstab,/dim))[0] 
  
  ;; Generate appropriate names
  spl_fil = sdss_getname(sdsstab,/spl,plate=plate,dir=cdir) 
  eig_fil = sdss_getname(sdsstab,/eig) 
  fnd_fil = sdss_getname(sdsstab,/abslin,dir=absdir) 

  ;; Loop by plates (since have to make directories)
  unq = uniq(plate,sort(plate))
  nunq = (size(unq,/dim))[0]

  for pp=0L,nunq-1 do begin 
     mtch = where(plate eq plate[unq[pp]],nmtch) 
     
     if nmtch eq 0 then begin 
        print,'sdss_reorganizeo: no match!' 
        return 
     endif 

     spawn,'mkdir -p '+cdir[unq[pp]]
     odir = '.gz '+cdir[unq[pp]] + '/.' 

     for ff=0,nmtch-1 do spawn,'rsync -ravu EIGCONTI/'+spl_fil[mtch[ff]]+odir 
     for ff=0,nmtch-1 do spawn,'rsync -ravu SPLCONTI/'+eig_fil[mtch[ff]]+odir 

     spawn,'mkdir -p '+absdir[unq[pp]]
     odir = '.gz '+absdir[unq[pp]] + '/.' 

     for ff=0L,nmtch-1 do spawn,'rsync -ravu ../ABSLIN/'+fnd_fil[mtch[ff]]+odir 

  endfor                        ; loop pp=nunq
  cd, cur_dir
end                             ; sdss_reorganize


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_srtbalqso, strct_fil, bal_fil, qso_fil
  ;; Take sdsscivstrct and sort between candidates in BAL QSO
  ;; sightlines and those in regular QSOs
  if n_params() ne 3 then begin
     print,'Syntax - sdss_srtbalqso, strct_fil, bal_fil, qso_fil'
     return
  endif 

  if size(strct_fil,/type) eq 8 then begin
     ;; Input is structure; sort and either return as structure or
     ;; write out to files
     strct = strct_fil 

     if tag_exist(strct,'BALFLG') then $
        qso = where(strct.balflg eq 0,complement=bal) $
     else qso = where(strct.flg_bal eq 0,complement=bal) 

     if qso[0] ne -1 then begin
        if size(qso_fil,/type) eq 7 then begin
           mwrfits,strct[qso],qso_fil,/create,/silent
           print,'sdss_srtbalqso: sorted QSOs into ',qso_fil
        endif else qso_fil = strct[qso]
     endif 
     if bal[0] ne -1 then begin
        if size(bal_fil,/type) eq 7 then begin
           mwrfits,strct[bal],bal_fil,/create,/silent
           print,'sdss_srtbalqso: sorted BAL QSOs into ',bal_fil
        endif else bal_fil = strct[bal]
     endif 
  endif else begin
     ;; Input is file which may have two extensions 
     done = 0
     ext = 1
     while not done do begin
        strct = xmrdfits(strct_fil,ext,/silent)
        if size(strct,/type) ne 8 then begin
           done = 1             ; EOF
           continue
        endif 

        if tag_exist(strct,'BALFLG') then $
           ;; Either {sdsscivstrct}
           qso = where(strct.balflg eq 0,complement=bal) $
           ;; Else {qalcharstrct}
        else qso = where(strct.flg_bal eq 0,complement=bal)  
        
        if qso[0] ne -1 then begin
           mwrfits,strct[qso],qso_fil,create=(ext mod 2),/silent
           print,'sdss_srtbalqso: sorted QSOs into ',qso_fil,$
                 ' ext = '+strtrim(ext,2)
        endif 
        if bal[0] ne -1 then begin
           mwrfits,strct[bal],bal_fil,create=(ext mod 2),/silent
           print,'sdss_srtbalqso: sorted BAL QSOs into ',bal_fil,$
                 ' ext = '+strtrim(ext,2)
        endif 

        ext = ext + 1
     endwhile
  endelse 

end                             ; sdss_srtbalqso



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntdr7summ, help=help
  ;; Print out some facts
  if keyword_set(help) then begin
     print,'Syntax - sdss_prntdr7summ [/help]'
     return
  endif 

  sdsstab = sdss_getqsostrct()
  file = sdss_getname(sdsstab,/abslin) ; dummy name
  wvspecmnx = sdss_getspecwave()
  nsdss = (size(sdsstab,/dim))[0]
  print,strtrim(nsdss,2)+' (100%) -- Total number of DR7 QSOs '
  nsdss = float(nsdss); for fractions
  
  balflg = sdss_fndbalqso(sdsstab,/full,count=nbal,/flgonly)
  print,' '+strtrim(nbal,2)+' ('+string(nbal/nsdss*100.,format='(f4.1)')+'%)'+$
        ' -- Number of BAL QSOs '
  
  ;; various redshift cuts reflect sdss_civsearch/sdss_fndciv 
  ;; ;;;;;;;;;
  ;; CIV
  ;; ;;;;;;;;;
  civ = dblt_retrieve('CIV')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_civ,dblt_name=civ,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_civ[*,0] ge wvspecmnx[0] and $
             wvobs_civ[*,1] le wvspecmnx[1] and $
             wvobs_civ[*,0] lt wvobs_civ[*,1],ngd)
  zmin = min(wvobs_civ[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+$
        ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+$
        ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CIV with Lya
  ;; ;;;;;;;;;
  lya = dblt_retrieve('Lya')  
  lya.wvII = civ.wvI
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_civ_lya,dblt_name=lya,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_civ_lya[*,0] ge wvspecmnx[0] and $
             wvobs_civ_lya[*,1] le wvspecmnx[1] and $
             wvobs_civ_lya[*,0] lt wvobs_civ_lya[*,1],ngd)
  zmin = min(wvobs_civ_lya[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ_lya[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and Lya coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and Lya coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CIV & SiIV: this seems wrong...
  ;; ;;;;;;;;;
  siiv = dblt_retrieve('SiIV')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_siiv,dblt_name=siiv,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_siiv[*,0] ge wvspecmnx[0] and $
             wvobs_siiv[*,1] le wvspecmnx[1] and $
             wvobs_siiv[*,0] lt wvobs_siiv[*,1] and $
             wvobs_civ[*,0] ge wvspecmnx[0] and $
             wvobs_civ[*,1] le wvspecmnx[1] and $
             wvobs_civ[*,0] lt wvobs_civ[*,1],ngd)
  zmin = min(wvobs_civ[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and SiIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and SiIV coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; SiIV
  ;; ;;;;;;;;;
  gd = where(wvobs_siiv[*,0] ge wvspecmnx[0] and $
             wvobs_siiv[*,1] le wvspecmnx[1] and $
             wvobs_siiv[*,0] lt wvobs_siiv[*,1],ngd)
  zmin = min(wvobs_siiv[gd,0]/siiv.wvI-1.)
  zmax = max(wvobs_siiv[gd,1]/siiv.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SiIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zsiiv < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SiIV coverage (excl. BAL QSOs) '


  ;; ;;;;;;;;;
  ;; SiV with Lya
  ;; ;;;;;;;;;
  lya.wvII = siiv.wvI
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_siiv_lya,dblt_name=lya,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_siiv_lya[*,0] ge wvspecmnx[0] and $
             wvobs_siiv_lya[*,1] le wvspecmnx[1] and $
             wvobs_siiv_lya[*,0] lt wvobs_siiv_lya[*,1],ngd)
  zmin = min(wvobs_siiv_lya[gd,0]/siiv.wvI-1.)
  zmax = max(wvobs_siiv_lya[gd,1]/siiv.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SIIV and Lya coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zsiiv < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)',+$
        ' -- Number with SIIV and Lya coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CaII
  ;; ;;;;;;;;;
  caii = dblt_retrieve('CaII')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_caii,dblt_name=caii,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_caii[*,0] ge wvspecmnx[0] and $
             wvobs_caii[*,1] le wvspecmnx[1] and $
             wvobs_caii[*,0] lt wvobs_caii[*,1],ngd)
  zmin = min(wvobs_caii[gd,0]/caii.wvI-1.)
  zmax = max(wvobs_caii[gd,1]/caii.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CaII coverage '+$
        string(zmin,zmax,format='(f5.3," < zcaii < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+ ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CaII coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; MgII
  ;; ;;;;;;;;;
  mgii = dblt_retrieve('MgII')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_mgii,dblt_name=mgii,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_mgii[*,0] ge wvspecmnx[0] and $
             wvobs_mgii[*,1] le wvspecmnx[1] and $
             wvobs_mgii[*,0] lt wvobs_mgii[*,1],ngd)
  zmin = min(wvobs_mgii[gd,0]/mgii.wvI-1.)
  zmax = max(wvobs_mgii[gd,1]/mgii.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with MgII coverage '+$
        string(zmin,zmax,format='(f5.3," < zmgii < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+ ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with MgII coverage (excl. BAL QSOs) '


  print,'CIV cut excludes OI/SiII "forest."'
  print,'SiIV cut excludes Lya forest.'
  print,'CaII cut exlcudes HVCs (up to 5000 km/s).'
  print,'MgII for z > 0.'

end                             ; sdss_prntdr7summ


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcparalleljob, sdsssum, processor, count=count
  ;; Paired with sdss_runparallelsrch.sh can divide big structure and
  ;; run sdss_fndlin and sdss_civsearch on multiple processors
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcparalleljob(sdsssum, processor, [count=]'
     return, -1
  endif 

  if size(sdsssum,/type) eq 8 or (size(sdsssum,/dim))[0] gt 0 then $
     sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]   ; if full DR7, 105783
  nsub = round(nsdss/float(processor[1]))  ; if processor[1]=4, nsub = 26445
  istrt = (processor[0]-1)*nsub 
  ;; Stuff the remainder into the last processor job
  if processor[0] eq processor[1] then iend = nsdss - 1 $
  else iend = processor[0]*nsub - 1
  count = iend - istrt + 1L
  
  return,[istrt,iend]                        ; make sure to exit for script
end                             ; sdss_calcparalleljob()


pro sdss_catparalleljob, cand_fil, processor
  ;; Paired with sdss_runparallelsrch.sh can combine sub-structures  
  ;; from sdss_civsearch into one
  if n_params() ne 2 then begin
     print,'Syntax - sdss_catparalleljob, cand_fil, processor'
     return
  endif
  
  ;; Load conventional name from sdss_runparallelsrch.sh
  sub_fil = cand_fil+'.'+strtrim((lindgen(processor[1])+1),2)
  for pp=1,processor[1] do begin
     strct = xmrdfits(sub_fil[pp-1],1,/silent)

     if pp eq 1 then begin
        if size(strct,/type) eq 8 then allciv = strct ; may be nothing
     endif else begin
        if size(strct,/type) eq 8 then begin
           if keyword_set(allciv) then allciv = [allciv,strct] $
           else allciv = strct
        endif 
     endelse 

     strct = xmrdfits(sub_fil[pp-1],2,/silent)

     if pp eq 1 then begin
        if size(strct,/type) eq 8 then allbroadciv = strct ; may be nothing
     endif else begin
        if size(strct,/type) eq 8 then begin
           if keyword_set(allbroadciv) then allbroadciv = [allbroadciv,strct] $
           else allbroadciv = strct
        endif 
     endelse 
  endfor                        ; loop pp=processor[1]
  
  if keyword_set(allciv) then begin
     allciv = sdss_srtcivstrct(allciv)
     mwrfits,allciv,cand_fil,/create,/silent
  endif else mwrfits,-1,cand_fil,/create,/silent
  if keyword_set(allbroadciv) then begin
     allbroadciv = sdss_srtcivstrct(allbroadciv)
     mwrfits,allbroadciv,cand_fil,/silent ; ext=2
  endif else mwrfits,-1,cand_fil,/silent ; ext=2
  spawn,'gzip -f '+cand_fil
  print,'sdss_catparalleljob: created ',cand_fil
  
  return                        ; make sure to exit for script
end                             ; sdss_catparalleljob



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_histogram, data, n_per_bin, min=min, max=max, scale=scale, $
                         location=location, szloc=szloc, nbin=nbin, $
                         force=force, debug=debug
  ;; Histogram data with un-equal bins in order to get (roughly) the
  ;; same number per bin
  if n_params() ne 2 then begin
     print,'Syntax - sdss_histogram(data, n_per_bin, [min=, max=, scale=, '
     print,'                        location=, szloc=, nbin=, /force, /debug])'
     return,-1
  endif 

  srtdata = data[sort(data)]
  if not keyword_set(min) then begin
     min = min(data) 
     if keyword_set(scale) then $
        min = floor(min * scale)/scale
  endif else begin
     ;; Must trim
     mn = min(srtdata-min,imn,/abs)
     if srtdata[imn] lt min then $
        imn++                   ; go up one
     srtdata = srtdata[imn:*]
  endelse 
  if not keyword_set(max) then begin
     max = max(data)
     if keyword_set(scale) then $
        max = ceil(max * scale)/scale
  endif else begin
     ;; Must trim
     mx = min(srtdata-max,imx,/abs)
     if srtdata[imx] gt max then $
        imx--                   ; go down one
     srtdata = srtdata[0:imx]
  endelse 
  ndata = (size(srtdata,/dim))[0]

  ;; set up output arrays
  nbin = ceil(ndata/float(n_per_bin))
  histogram = lonarr(nbin,/nozero)
  location = dblarr(nbin,/nozero)
  szloc = dblarr(nbin,/nozero)

  ;; lay down a first pass
  indx = n_per_bin * lindgen(nbin)
  location = srtdata[indx]
  szloc = shift(location,-1) - location
  szloc[nbin-1] = srtdata[ndata-1] - location[nbin-1] 
  histogram = [replicate(n_per_bin,nbin-1),0]
  histogram[nbin-1] = ndata - total(histogram)
  
  ;; Check the last bin
  percent_off = histogram[nbin-1]/float(n_per_bin)
  if percent_off lt 0.1 then begin ; arbitrary number
     ;; Shove last bin into second to last 
     if keyword_set(debug) then $
        print,'sdss_histogram() debug: last bin has < 10% of previous bin; lump.',$
              n_per_bin
     histogram[nbin-2] = total(histogram[nbin-2:*])
     szloc[nbin-2] = total(szloc[nbin-2:*])
     nbin--
     ;; Trim
     histogram = histogram[0:nbin-1]
     location = location[0:nbin-1]
     szloc = szloc[0:nbin-1]
  endif else $
     if keyword_set(debug) then $
        print,'sdss_histogram() debug: last bin has >= 10% of previous bin; leave.',$
              n_per_bin
  
  if not keyword_set(force) then begin
     ;; use the mean binsize in one more iteration
     nmean = mean(histogram)
     histogram = sdss_histogram(srtdata,round(nmean),min=min,max=max,nbin=nbin,$
                                location=location,szloc=szloc,scale=0,$
                                debug=debug,/force)
     if keyword_set(debug) then $
        print,'sdss_histogram() debug: changing n_per_bin from, to',n_per_bin,$
              round(nmean)
  endif
  
  if keyword_set(scale) then begin
     ;; Want scale to give roughly equal binning but to have truncation
     ;; in significant figures
     ;; Do some rounding with the best-case scenario given
     nwloc = round(location*scale)/double(scale)
     nwloc[0] = floor(location[0]*scale)/double(scale) ; include all
     nwloc[nbin-1] = ceil(location[nbin-1]*scale)/double(scale) ; include all
     nwszloc = shift(nwloc,-1) - nwloc
     nwszloc[nbin-1] = ceil(srtdata[ndata-1]*scale)/double(scale) - $
                       nwloc[nbin-1]
     nwhist = lonarr(nbin,/nozero)
     for ii=0L,nbin-1 do begin
        gd = where(srtdata ge nwloc[ii] and srtdata lt nwloc[ii]+nwszloc[ii],ngd)
        nwhist[ii] = ngd
     endfor                     ; loop ii=nbin
     if keyword_set(debug) then begin
        print,'Prior to scaling:'
        print,'Loc','Size','Number',format='(3(a12,1x))'
        printcol,location,szloc,histogram
     endif
     ;; Overwrite
     location = nwloc
     szloc = nwszloc
     histogram = nwhist

  endif else begin
     szloc[nbin-1] = szloc[nbin-1] + 1.e-5 ; little extra to handle floating point comparison
     if keyword_set(debug) then $
     print,'sdss_histogram() debug: forcing bins as is.'
  endelse


  if keyword_set(debug) then begin
     print,'Loc','Size','Number',format='(3(a12,1x))'
     printcol,location,szloc,histogram
  endif
  
  return, histogram
end                             ;  sdss_histogram()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_bintoloc, data, loc, szloc
  ;; complements sdss_histogram()
  if n_params() ne 3 then begin      
     print,'Syntax - sdss_bintoloc(data, loc, szloc)'
     return,-1
  endif 

  nbin = (size(loc,/dim))[0] > 1
  histogram = lonarr(nbin,/nozero)
  for ii=0L,nbin-1 do begin
     gd = where(data ge loc[ii] and data lt loc[ii]+szloc[ii],ngd)
     histogram[ii] = ngd
  endfor                        ; loop ii=nbin

  return, histogram
end                             ; sdss_bintoloc()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_mtchcivstrct, strct1_fil, strct2_fil, mtch_fil, diff_fil, $
                       log_fil=log_fil,final=final, force=force, dvtol=dvtol
  ;; Take sdsscivstrct and match on qso_name and wavelength limits...
  if n_params() ne 4 then begin
     print,'Syntax - sdss_mtchcivstrct, strct1_fil, strct2_fil, mtch_fil, diff_fil'
     print,'                            [log_fil=,/final,/force,dvtol]'
     return
  endif 
  ;; Defaults
  cinv = 1./3.e5                             ; km^-1 s
  if not keyword_set(log_fil) then $
     log_fil = strmid(mtch_fil,0,strpos(mtch_fil,'.',/reverse_search))+'.log'
  close,1
  openw,1,log_fil

  ;; Read in
  if size(strct1_fil,/type) eq 7 then $
     strct1 = xmrdfits(strct1_fil,1,/silent) $
  else strct1 = strct1_fil
  strct1.qso_name = strtrim(strct1.qso_name,2)

  if size(strct2_fil,/type) eq 7 then $
     strct2 = xmrdfits(strct2_fil,1,/silent) $
  else strct2 = strct2_fil
  strct2.qso_name = strtrim(strct2.qso_name,2)

  ;; Plates
;  plate1 = strmid(strct1.qso_name,strpos(strct1[0].qso_name,'-')+1,4)
;  strct1 = strct1[sort(plate1)]
;  plate2 = strmid(strct2.qso_name,strpos(strct2[0].qso_name,'-')+1,4)
;  strct2 = strct2[sort(plate2)]

  nstrct1 = (size(strct1,/dim))[0]
  nstrct2 = (size(strct2,/dim))[0]
  if nstrct1 gt nstrct2 then begin
     ;; Swap
     swap = 1                   ; true
     tmp = strct1
     strct1 = strct2
     strct2 = tmp
     print,'sdss_mtchcivstrct: swapping order of input files.'
     printf,1,'sdss_mtchcivstrct: swapping order of input files.'
     if size(srct2_fil,/type) eq 7 then $
        printf,1,'strct1 = ',strct2_fil
     if size(srct1_fil,/type) eq 7 then $
        printf,1,'strct2 = ',strct1_fil
  endif else begin
     swap = 0                   ; false
     if size(srct1_fil,/type) eq 7 then $
        printf,1,'strct1 = ',strct1_fil
     if size(srct2_fil,/type) eq 7 then $
        printf,1,'strct2 = ',strct2_fil
  endelse 
  nstrct1 = (size(strct1,/dim))[0]
  nstrct2 = (size(strct2,/dim))[0]

  tags = tag_names(strct1)
  if keyword_set(final) then begin
     ;; Use final zabs and wvlim
     ztag = where(tags eq 'ZABS_FINAL')
     wvlimtag = where(tags eq 'WVLIM_FINAL')
  endif else begin
     ;; Use original zabs and wvlim
     ztag = where(tags eq 'ZABS_ORIG')
     wvlimtag = where(tags eq 'WVLIM_ORIG')
  endelse 

  ;; Instantiate useful arrays
  mask1 = replicate(-1L,nstrct1)
  mask2 = replicate(-1L,nstrct2)
  wvobs1 = strct1.wrest[0]*(1+strct1.(ztag)[0])
  wvobs2 = strct2.wrest[0]*(1+strct2.(ztag)[0])

  ;; Match on wavelength bounds
  for ii=0L,nstrct1-1 do begin
     ;; Use wavelength bounds to maximize inclusion
     ;; Allow boundaries to be "fuzzy" if think bounds aren't
     ;; naturally big enough (default: dvtol = 0.)
     ;; Bug: Always adding 1 Ang to the wavelength bounds
     if keyword_set(dvtol) then $
        mtch = where(strct1[ii].qso_name eq strct2.qso_name and $
                     ((wvobs1[ii] ge strct2.(wvlimtag)[0,0]+(1-dvtol*cinv) and $
                       wvobs1[ii] le strct2.(wvlimtag)[0,1]+(1+dvtol*cinv)) or $
                      (wvobs2 ge strct1[ii].(wvlimtag)[0,0]+(1-dvtol*cinv) and $
                       wvobs2 le strct1[ii].(wvlimtag)[0,1]+(1+dvtol*cinv))) and $
                     mask2 eq -1,nmtch) $
     else $
        mtch = where(strct1[ii].qso_name eq strct2.qso_name and $
                     ((wvobs1[ii] ge strct2.(wvlimtag)[0,0] and $
                       wvobs1[ii] le strct2.(wvlimtag)[0,1]) or $
                      (wvobs2 ge strct1[ii].(wvlimtag)[0,0] and $
                       wvobs2 le strct1[ii].(wvlimtag)[0,1])) and $
                     mask2 eq -1,nmtch)         
     if nmtch gt 1 then begin
        printf,1,strct1[ii].qso_name,wvobs1[ii],$
               format='("Multiple matches to ",a15,2x,"wvobs = ",f7.2)'
        ;; Take closest
        mn = min(wvobs1[ii]-wvobs2[mtch],imn,/abs)
        name = strct2[mtch].qso_name
        name[imn] = '*'+name[imn] ; indidcate
        writecol,log_fil,name,wvobs2[mtch],filnum=1,$
                 fmt='(5x,a15,2x,"wvobs = ",f7.2)'
        ;; Set
        mtch = mtch[imn]
     endif else mtch = mtch[0]  ; which may be -1

     if mtch ne -1 then begin
        mask1[ii] = mtch
        mask2[mtch] = ii
     endif 
  endfor                        ; loop ii=nstrct1

  mtch1 = where(mask1 ge 0,nmtch1,complement=unmtch1,ncomplement=nunmtch1)
  mtch2 = where(mask2 ge 0,nmtch2,complement=unmtch2,ncomplement=nunmtch2)
  if keyword_set(force) and (nunmtch1 ne 0 or nunmtch2 ne 0) then begin
     ;; Secondary forced match by closest in sightline
     ;; First status
     printf,1,''
     printf,1,'sdss_mtchcivstrct: initial matches ',nmtch1
     printf,1,'sdss_mtchcivstrct: initial non-matches ',nunmtch1,nunmtch2

     if nunmtch1 ne 0 then begin
        printf,1,''
        printf,1,'Forced pairing between unmatched strct1 and all strct2.'
        frc_mask1 = replicate(' ',nstrct1)
        for ii=0L,nunmtch1-1 do begin
           los = where(strct1[unmtch1[ii]].qso_name eq strct2.qso_name)
           if los[0] eq -1 then begin
              printf,1,'sdss_mtchcivstrct: LOS not in strct2: ',$
                    strct1[unmtch1[ii]].qso_name

              if not keyword_set(qso1notin2) then $
                 qso1notin2 = strct1[unmtch1[ii]].qso_name $
              else qso1notin2 = [qso1notin2,strct1[unmtch1[ii]].qso_name]
              continue
           endif 
           mn = min(wvobs2[los]-wvobs1[unmtch1[ii]],imn,/abs)
           dv = 3.e5*mn/strct1[unmtch1[ii]].(ztag)[0]
           if abs(dv) lt 1.e3 then begin
              frc_mask1[unmtch1[ii]] = frc_mask1[unmtch1[ii]] + ' ' + $
                                       string(los[imn])
              printf,1,strct2[los[imn]].qso_name,wvobs2[los[imn]],$
                 round(dv),format='(5x,a15,2x,"wvobs2 = ",f7.5,2x,"dv21 = ",i6)'
           endif 
        endfor                  ; loop ii=nstrct1
        gd = where(frc_mask1 ne ' ')
        if gd[0] ne -1 then $
           writecol,log_fil,gd,' '+strct1[gd].qso_name,frc_mask1[gd],filnum=1
     endif                      ; nunmtch != 0

     if nunmtch2 ne 0 then begin
        printf,1,''
        printf,1,'Forced pairing between unmatched strct2 and all strct1.'
        frc_mask2 = replicate(' ',nstrct2)
        for ii=0L,nunmtch2-1 do begin
           los = where(strct2[unmtch2[ii]].qso_name eq strct1.qso_name)
           if los[0] eq -1 then begin
              printf,1,'sdss_mtchcivstrct: LOS not in strct1: ',$
                    strct2[unmtch2[ii]].qso_name
              if not keyword_set(qso2notin1) then $
                 qso2notin1 = strct2[unmtch2[ii]].qso_name $
              else qso2notin1 = [qso2notin1,strct2[unmtch2[ii]].qso_name]
              continue
           endif 
           mn = min(wvobs1[los]-wvobs2[unmtch2[ii]],imn,/abs)
           dv = 3.e5*mn/strct2[unmtch2[ii]].(ztag)[0]
           if abs(dv) lt 1.e3 then begin
              frc_mask2[unmtch2[ii]] = frc_mask2[unmtch2[ii]] + ' ' + $
                                       string(los[imn])
              printf,1,strct1[los[imn]].qso_name,wvobs1[los[imn]],$
                    round(dv),format='(5x,a15,2x,"wvobs1 = ",f7.5,2x,"dv12 = ",i6)'
           endif 
        endfor                  ; loop ii=nstrct1
        gd = where(frc_mask2 ne ' ')
        if gd[0] ne -1 then $
           writecol,log_fil,gd,' '+strct2[gd].qso_name,frc_mask2[gd],filnum=1
     endif                      ; nunmtch != 0

     ;; New conditionals
     mtch1 = where(mask1 ge 0 or frc_mask1 ne ' ',nmtch1,$
                   complement=unmtch1,ncomplement=nunmtch1)
     mtch2 = where(mask2 ge 0 or frc_mask2 ne ' ',nmtch2,$
                   complement=unmtch2,ncomplement=nunmtch2)

  endif                         ; /force
  
  ;; Switch back order if necessary so can automatically use results
  ;; without knowing "the answer"
  if keyword_set(swap) then begin
     tmp = strct1
     strct1 = strct2
     strct2 = tmp

     tmp = mask1
     mask1 = mask2
     mask2 = tmp

     tmp = mtch1
     mtch1 = mtch2
     mtch2 = tmp
     tmp = nmtch1
     nmtch1 = nmtch2
     nmtch2 = tmp

     tmp = unmtch1
     unmtch1 = unmtch2
     unmtch2 = tmp
     tmp = nunmtch1
     nunmtch1 = nunmtch2
     nunmtch2 = tmp

     print,'sdss_mtchcivstrct: swapping order of input files back.'
     printf,1,'sdss_mtchcivstrct: swapping order of input files back.'
     if size(srct2_fil,/type) eq 7 then $
        printf,1,'strct1 = ',strct1_fil
     if size(srct1_fil,/type) eq 7 then $
        printf,1,'strct2 = ',strct2_fil
  endif  

  
  ;; Output
  printf,1,''
  if nmtch1 ne nmtch2 and not keyword_set(force) then $
     stop,'sdss_mtchcivstrct: matched indices not aligned.'

  if nmtch1 ne 0 then begin
     mwrfits,strct1[mtch1],mtch_fil,/create,/silent ; ext = 1
     ;; Keep it in the same order
     mwrfits,strct2[mask1[mtch1]],mtch_fil,/silent ; ext = 2
     mwrfits,mtch1,mtch_fil,/silent ; ext = 3
     mwrfits,mask1[mtch1],mtch_fil,/silent ; ext = 4
     spawn,'gzip -f '+mtch_fil
     print,'sdss_mtchcivstrct: created 4-extensions of ',mtch_fil,nmtch1
     printf,1,'sdss_mtchcivstrct: created 4-extensions of ',mtch_fil,nmtch1
  endif

  if nunmtch1 ne 0 then begin
;     ;; Sort
;     tmp = sdss_srtcivstrct(strct1[unmtch1])
     mwrfits,strct1[unmtch1],diff_fil,/create,/silent ; ext = 1
     diff = sdss_srtcivstrct(strct1[unmtch1])
  endif else mwrfits,[0],diff_fil,/create,/silent ; ext = 1
  if nunmtch2 ne 0 then begin
;     tmp = sdss_srtcivstrct(strct2[unmtch2])
     mwrfits,strct2[unmtch2],diff_fil,/silent            ; ext = 2
     diff2 = sdss_srtcivstrct(strct2[unmtch2])
     
     ;; Other info
     dvabs2 = (diff2.zabs_orig[1]-diff2.zabs_orig[0])/(1+diff2.zabs_orig[0])*3.e5
     gd = where(abs(dvabs2) le 150.,ngd)
     print,'sdss_mtchcivstrct: diff2 with dvabs <= 150 km/s',ngd
     printf,1,'sdss_mtchcivstrct: diff2 with dvabs <= 150 km/s',ngd
  endif else mwrfits,[0],diff_fil,/silent ; ext = 2

  if nunmtch1 ne 0 then mwrfits,unmtch1,diff_fil,/silent $ ; ext = 3
  else mwrfits,[0],diff_fil,/silent
  if nunmtch2 ne 0 then mwrfits,unmtch2,diff_fil,/silent $ ; ext = 4
  else mwrfits,[0],diff_fil,/silent

  spawn,'gzip -f '+diff_fil
  print,'sdss_mtchcivstrct: created 4-extensions of ',$
        diff_fil,nunmtch1,nunmtch2
  printf,1,'sdss_mtchcivstrct: created 4-extensions of ',$
         diff_fil,nunmtch1,nunmtch2

  save,/all,filename=strmid(mtch_fil,0,strpos(log_fil,'.',/reverse_search))+'.sav'

  close,1
  print,'sdss_mtchcivstrct: created ',log_fil
  

end                             ; sdss_mtchcivstrct



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntcandsumm, cand_fil, final=final, outfil=outfil, logarr=logarr, $
                       by_rtg=by_rtg, silent=silent
  ;; Print # of candidates, # LOS, z and EW median and range (for *_orig or
  ;; *_final) and division by rating
  if n_params() ne 1 then begin
     print,'Syntax - sdss_prntcandsumm, cand_fil, [/final, outfil=, logarr=, /by_rtg, /silent'
     return
  endif 

  logarr = strarr(100)
  ii = 0

  if size(cand_fil,/type) eq 7 then cand = xmrdfits(cand_fil,1,/silent) $
  else cand = cand_fil
  ncand = (size(cand,/dim))[0]
  dblt = dblt_retrieve(cand[0].wrest[0])
  ion = dblt.ion

  tags = tag_names(cand)
  if keyword_set(final) then begin
     ztag = where(tags eq 'ZABS_FINAL')
     ewtag = where(tags eq 'EW_FINAL')
  endif else begin
     ztag = where(tags eq 'ZABS_ORIG')
     ewtag = where(tags eq 'EW_ORIG')
  endelse 
  
  unq = uniq(cand.qso_name,sort(cand.qso_name))
  nqso = (size(unq,/dim))[0]

  logarr[ii] = '# '+ion+' candidates = '+string(ncand,format='(i6)')+$
               '; # LOS = '+string(nqso,format='(i5)')
  ii++

  logarr[ii] = 'zabs median, min, max: '+$
               string(median(cand.(ztag)[0],/even),$
                      min(cand.(ztag)[0],max=mx),mx,$
                      format='(3(f7.5,2x))')
  ii++

  logarr[ii] = 'EW median, min, max: '+$
               string(median(cand.(ewtag)[0],/even),$
                      min(cand.(ewtag)[0],max=mx),mx,$
                      format='(3(f8.3,2x))')
  ii++  

  if keyword_set(by_rtg) then begin
     ;; Will this program recursively
     unq = uniq(cand.rating[0],sort(cand.rating[0]))
     nrtg = (size(unq,/dim))[0]
     for rr=0,nrtg-1 do begin
        logarr[ii] = ion+' with rating = '+strtrim(cand[unq[rr]].rating[0],2)
        ii++
        sub = where(cand.rating[0] eq cand[unq[rr]].rating[0])
        sdss_prntcandsumm,cand[sub],final=final,logarr=subarr,/silent
        gd = where(subarr ne '',ngd)
        logarr[ii:(ii+ngd-1)] = subarr[gd]
        ii = ii+ngd
     endfor                     ; loop rr=nrtg
  endif                         ; /by_rtg


  if keyword_set(outfile) then begin
     openw,1,outfile,append=ss
     printf,1,logarr[0:ii-1],format='(a)'
     close,1
     print,'sdss_prntcandsumm: created ',outfile
  endif else if not keyword_set(silent) then $
     print,logarr[0:ii-1],format='(a)' 

end                             ;  sdss_prntcandsumm


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_mtchcandlist, cand_fil, list_fil, silent=silent, outfil=outfil
  ;; Given a candidate structure and a list of sightlines, find the
  ;; candidates in the list
  if n_params() ne 2 then begin
     print,'Syntax - sdss_mtchcandlist(cand_fil, list_fil, [/silent,outfil=])'
     return,-1
  endif 

  if size(cand_fil,/type) eq 7 then cand = xmrdfits(cand_fil,1,/silent) $
  else cand = cand_fil
  cand.qso_name = strtrim(cand.qso_name,2)
  cand = cand[sort(cand.qso_name)]

  ncand = (size(cand,/dim))[0]
  eqsowcand = uniq(cand.qso_name) ; ending index
  nqsowcand = (size(eqsowcand,/dim))[0]
  bqsowcand = [0,eqsowcand[0:nqsowcand-2]+1]; beginning index


  readcol,list_fil,spec_fil,skip=1,/silent,format='a'
  nspec = (size(spec_fil,/dim))[0]
  dum = sdss_getname(spec_fil,/spec,root=qsoinlist)

  nloop = nqsowcand < nspec     ; minimal looping
  for qq=0L,nloop-1 do begin
     ;; Search two ways
     if nqsowcand lt nspec then begin
        mtch = where(cand[eqsowcand[qq]].qso_name eq qsoinlist,nmtch)
        ;; Save candidates
        if nmtch eq 0 then begin
           if not keyword_set(silent) then $
              print,'sdss_mtchcandlist(): QSO w/ candidate not in list: ',$
                    cand[eqsowcand[qq]].qso_name
        endif else begin
           if keyword_set(subcand) then $
              subcand = [subcand,cand[bqsowcand[qq]:eqsowcand[qq]]] $
           else subcand = cand[bqsowcand[qq]:eqsowcand[qq]]
        endelse 
     endif else begin
        mtch = where(cand[eqsowcand].qso_name eq qsoinlist[qq],nmtch)
        if nmtch eq 0 then begin
           if not keyword_set(silent) then $
              print,'sdss_mtchcandlist(): QSO in list_fil has no candidate: ',$
                    qsoinlist[qq]
        endif else begin
           if nmtch ne 1 then $
              stop,'sdss_mtchcandlist(): should not be here...'

           if keyword_set(subcand) then $
              subcand = [subcand,cand[bqsowcand[mtch[0]]:eqsowcand[mtch[0]]]] $
           else subcand = cand[bqsowcand[mtch[0]]:eqsowcand[mtch[0]]]
        endelse 
     endelse 
  endfor                        ; loop qq=nloop
  
  subcand = sdss_srtcivstrct(subcand)
  if keyword_set(outfil) then begin
     mwrfits,subcand,outfil,/create,/silent
     spawn,'gzip -f '+outfil
     if not keyword_set(silent) then $
        print,'sdss_mtchcandlist(): created ',outfil
  endif 

  return,subcand

end                             ; sdss_mtchcandlist()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_wrnotefil, civstrct_fil, note_fil, notes=notes, silent=silent
  if n_params() ne 2 then begin
     print,'Syntax - sdss_wrnotefil, civstrct_fil, note_fil, [notes=]'
     return
  endif 

  if size(civstrct_fil,/type) eq 7 then $
     civstr = xmrdfits(civstrct_fil,1,/silent) $
  else civstr = civstrct_fil
  ncivstrct = (size(civstr,/dim))[0]

  if keyword_set(notes) then begin
     ;; Add it on
     num = (size(notes,/dim))[0]
     if num eq 0 then $
        ;; One message for all
        tmp = replicate(notes,ncivstrct) $
     else if num ne ncivstrct then $
        stop,'sdss_wrnotefil: notes must have same size as structure' $
     else tmp = notes

     blank = where(strtrim(civstr.notes,2) eq '',complement=filled)
     if blank[0] ne -1 then civstr[blank].notes = tmp[blank]
     if filled[0] ne -1 then $
        civstr[filled].notes = strtrim(civstr[filled].notes,2)+'; '+tmp[filled]
  endif  
  civstr.notes = strtrim(civstr.notes,2)

  gd = where(civstr.zabs_final[0] gt 0.)
  zabs = civstr.zabs_orig[0]
  if gd[0] ne -1 then zabs[gd] = civstr[gd].zabs_final[0]
  writecol,note_fil,lindgen(ncivstrct),civstr.qso_name,$
           zabs,civstr.rating[0],civstr.notes,$
           fmt='(i7,1x,a16,1x,f7.5,1x,i2,":",2x,a)'
  if not keyword_set(silent) then print,'sdss_wrnotefil: saved ',note_fil

end                             ; sdss_wrnotefil


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_rdnotefil, note_fil, index, qso_name, zabs, rating, notes, $
                    quick=quick
  if n_params() ne 6 and not keyword_set(quick) then begin
     print,'Syntx - sdss_rdnotefil, note_fil, index, qso_name, zabs, '
     print,'                        rating, notes'
     print,' or     sdss_rdnotefil, note_fil, notes, /quick'
     return 
  endif 

  readcol, note_fil, dum1, notes, format='a,a', delimiter=':', /silent
  nqso = (size(dum1,/dim))[0]

  if keyword_set(quick) then begin
     ;; Just read the notes
     index = notes
     return
  endif 

  tmpnote_fil = 'sdss_rdnotefil_tmp.notes'
  close,1
  openw,1,tmpnote_fil
  printf,1,dum1,format='(a)'
  close,1

  readcol, tmpnote_fil, index, qso_name, zabs, rating, $
           format='i,a,f,i', /silent
  spawn,'\rm '+tmpnote_fil
 
end                             ; sdss_rdnotefil


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpstrct, oldstrct_fil, newstrct_tmplt, excl_tag=excl_tag, $
                       nshift=nshift, verbose=verbose
  ;; Copy all elements of one civstrct to another (primarily used to
  ;; add new NOTES tag)
  if n_params() ne 2 then begin
     print,'Syntax - sdss_cpstrct(oldstrct_fil, newstrct_tmplt, '
     print,'                     [excl_tag=, nshift=, /verbose])'
     return,-1
  endif 

  if not keyword_set(nshift) then nshift = 0

  if size(oldstrct_fil,/type) eq 8 then oldstrct = oldstrct_fil $
  else oldstrct = xmrdfits(oldstrct_fil,1,/silent)
  nstrct = (size(oldstrct,/dim))[0] ; won't work if oldstrct is 1
  oldtags = tag_names(oldstrct[0])
  noldtags = (size(oldtags,/dim))[0]

  if size(newstrct_tmplt,/type) ne 8 then $
        stop,'sdss_cpstct() stop: newstrct_tmplt must be structure'
  newstrct = replicate(newstrct_tmplt,nstrct)

  newtags = tag_names(newstrct[0])
  nnewtags = (size(newtags,/dim))[0]

  for tt=0L,nnewtags-1 do begin
     mtch = where(newtags[tt] eq oldtags)
     if mtch[0] eq -1 then begin
        if keyword_set(verbose) then $
           print,'sdss_cpstrct(): tag '+newtags[tt]+' DNE in old structure'
        continue
     endif else begin
        if keyword_set(excl_tag) then $
           test = where(newtags[tt] eq strupcase(excl_tag)) $
        else test = [-1]

        if test[0] eq -1 then begin
           ;; Must monitor array size changes
           ndim_new = (size(newstrct[0].(tt),/dim))[0]
           ndim_old = (size(oldstrct[0].(mtch[0]),/dim))[0]

           if ndim_new eq ndim_old and not keyword_set(nshift) then begin
              ;; dimensionality is [ntagarr, nstrct]
              newstrct.(tt) = oldstrct.(mtch[0])

           endif else begin 

              if ndim_new ge ndim_old then begin
                 ndim2_new = (size(newstrct[0].(tt),/n_dim))[0] ; e.g., wvlim_* would == 2
                 ndim2_old = (size(oldstrct[0].(mtch[0]),/n_dim))[0]
                 if ndim2_old gt ndim2_new then $
                    stop,'sdss_cpstrct() stop: new number of dimensions smaller than old for tag ',$
                         newtags[tt]
                 
                 for ff=0,ndim2_old-1 do begin
                    ;; faster to have less looping outside
                    for ee=0,ndim_old-1 do begin

                       if keyword_set(nshift) then begin
                          ;; IDL convention is nshift < 0, then new
                          ;; 0th element becomes old 1st element
                          ;; and since ndim_new >= ndim_old, refer to
                          ;; ndim_old always
                          if ee-nshift lt 0 then begin
                             ;; nshift > 0 
                             ;; ee --> ndim_old-(ee-nshift)
                             ;; e.g., wrap beginning to end
                             if ndim2_old eq 1 then $
                                newstrct.(tt)[ee] = $
                                oldstrct.(mtch[0])[ndim_old-nshift-ee] $
                             else $
                                newstrct.(tt)[ee] = $
                                oldstrct.(mtch[0])[ndim_old-nshift-ee,ff] 

                          endif else begin

                             if ee-nshift ge ndim_old then begin
                                ;; nshift < 0
                                ;; ndim_old-(ee-nshift) --> ee
                                ;; e.g., wrap end to beginning
                                if ndim2_old eq 1 then $
                                   newstrct.(tt)[ee] = $
                                   oldstrct.(mtch[0])[(ee-nshift) mod ndim_old] $
                                else $
                                   newstrct.(tt)[ee,ff] = $
                                   oldstrct.(mtch[0])[(ee-nshift) mod ndim_old,ff] 
                             
                             endif else begin
                                ;; 0 <= ee-nshift <= ndim_new-1
                                ;; e.g., ee is old ee-nshift

                                if ndim2_old eq 1 then $
                                   newstrct.(tt)[ee] = $
                                   oldstrct.(mtch[0])[ee-nshift] $
                                else $
                                   newstrct.(tt)[ee,ff] = $
                                   oldstrct.(mtch[0])[ee-nshift,ff]
                             endelse

                          endelse 

                       endif else begin ; nshift = 0 
                          ;; element-by-element copy
                          if ndim2_old eq 1 then $
                             newstrct.(tt)[ee] = oldstrct.(mtch[0])[ee] $ 
                          else $
                             newstrct.(tt)[ee,ff] = $
                             oldstrct.(mtch[0])[ee,ff]
                       endelse

                    endfor      ; loop ee=ndim_old

                 endfor         ; loop ff=ndim2_old

              endif else begin
                 stop,'sdss_cpstrct() stop: new dimensions smaller than old dimensions for tag ',$
                      newtags[tt]
              endelse           ; ndim_new < ndim_old

           endelse              ; ndim_new != ndim_old

        endif else if keyword_set(verbose) then $
           print,'sdss_cpstrct(): excluding tag '+excl_tag[test[0]]
     endelse

  endfor                        ; loop tt=ntags

  return, newstrct
end                             ; sdss_cpstrct()
                    

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_srtvpstrct, vpstrct, recurs=recurs
  if n_params() ne 1 then begin
     print,'sdss_srtvpstrct(vpstrct)'
     return,-1
  endif 
  ;; May have duplicate id_sys
  id_sys = vpstrct[uniq(vpstrct.id_sys,sort(vpstrct.id_sys))].id_sys
  id_sys = id_sys[sort(id_sys)]
  nsys = (size(id_sys,/dim))[0] > 1 ; foil singularity
  mxid_sys = max(id_sys) + 1    ; fix duplicate system IDs

  for ss=0L,nsys-1 do begin
     ;; Find all with system ID
     sub = where(id_sys[ss] eq vpstrct.id_sys,nsub)
     tmpstrct = vpstrct[sub]
     ndblt = tmpstrct[0].ncomp+1 ; # of doublets
     
     if nsub ne 2*ndblt then begin
        if keyword_set(recurs) then $ ; emergency stop
           stop,'sdss_srtvpstrct() stop: WARNING! should never get here'

        ;; Test if more than two systems with same number
        test = where(tmpstrct.id_comp eq 0 and tmpstrct.id_lin eq 1,ntest)
        print,'sdss_srtvpstrct: duplicate id_sys = ',id_sys[ss],ntest
        
        for tt=1,ntest-1 do begin ; leave first with orig id_Sys
           ;; There's been a problem; so increment systems in id_sys 
           ;; but have to grep on ncomp and hope they aren't the same!
           gd = where(tmpstrct[test[tt]].ncomp eq tmpstrct.ncomp)
           tmpstrct[gd].id_sys = mxid_sys
           mxid_sys++           ; for next case
        endfor                  ; loop tt=ntest

        tmpstrct = sdss_srtvpstrct(tmpstrct,/recurs) ; now sort this subset

     endif else begin           ; called recursively
        ;; Pull out wvI and wvII indices and sort them by component ID 
        gdI = where(tmpstrct.id_lin eq 1,ngdI,complement=gdII)
        srtI = gdI[sort(tmpstrct[gdI].id_comp)]
        srtII = gdII[sort(tmpstrct[gdII].id_comp)]
        ;; Set the wvI index and put everything back
        rng = 2 * lindgen(tmpstrct[0].ncomp+1)
        tmpstrct[rng] = tmpstrct[srtI]
        tmpstrct[rng+1] = tmpstrct[srtII]
     endelse                    ; no recursive call needed
     
     ;; Append
     if ss eq 0 then srtstrct = tmpstrct $
     else srtstrct = [srtstrct,tmpstrct]
  endfor                        ; loop ss=nsys
   
  return, srtstrct
end                             ; sdss_srtvpstrct


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpvp2mcstrct, vpstrct, mcstrct
  ;; Instantiates the (input) system and component values from the
  ;; sdss_genprof Voigt Profile structure to the Monte Carlo structure
  ;; Instantiates tags: ID_SYS, NSYS, WREST, WVLIM_INPUT, ZABS_INPUT, 
  ;;  NCOLM_INPUT, B_INPUT, EW_INPUT
  if n_params() lt 1 then begin
     print,'Syntax - sdss_cpvp2mcstrct(vpstrct,[mcstrct])'
     return,-1
  endif 
  vpstr = sdss_srtvpstrct(vpstrct)         ; passing in ad hoc array
  index_sys = where(vpstr.id_comp eq 0 and $
                    vpstr.id_lin eq 1,nsys)

  if not keyword_set(mcstrct) then $
     mcstrct = {sdssmcstrct}
  nsysmax = (size(mcstrct.id_sys,/dim))[0]
  ncompmax = (size(mcstrct.zabs_input,/dim))[1] - 1 ; really, comp.
  
  ;; Checks
  if nsys gt nsysmax then $
     stop,'sdss_cpvp2mcstrct() stop: number of systems too large for sdssmcstrct',nsys
  if max(vpstr.ncomp+1,imx) gt ncompmax then $
     stop,'sdss_cpvp2mcstrct() stop: number of components too large for sdssmcstrct',vpstr[imx].ncomp+1

  mcstrct.id_sys[0:nsys-1] = vpstr[index_sys].id_sys
  mcstrct.nsys = nsys
  for ss=0,nsys-1 do begin
     mcstrct.ncomp[ss] = vpstr[index_sys[ss]].ncomp
     rng = index_sys[ss] + 2*indgen(mcstrct.ncomp[ss]+1) 
     iwvI = rng[0]
     
     mcstrct.wrest[ss,0] = vpstr[iwvI].wrest   ; wvI
     mcstrct.wrest[ss,1] = vpstr[iwvI+1].wrest ; wvII

     mcstrct.wvlim_input[ss,0,0] = vpstr[iwvI].wvlim_sys[0] ; wvI
     mcstrct.wvlim_input[ss,0,1] = vpstr[iwvI].wvlim_sys[1] ; wvI
     mcstrct.wvlim_input[ss,1,0] = vpstr[iwvI+1].wvlim_sys[0] ; wvII
     mcstrct.wvlim_input[ss,1,1] = vpstr[iwvI+1].wvlim_sys[1] ; wvII

     mcstrct.zabs_input[ss,0] = vpstr[iwvI].z_sys
     mcstrct.zabs_input[ss,1:mcstrct.ncomp[ss]+1] = vpstr[rng].zabs

     mcstrct.ncolm_input[ss,0] = vpstr[iwvI].n_sys
     mcstrct.ncolm_input[ss,1:mcstrct.ncomp[ss]+1] = vpstr[rng].n

     mcstrct.b_input[ss,0] = sqrt(total(vpstr[rng].b^2))
     mcstrct.b_input[ss,1:mcstrct.ncomp[ss]+1] = vpstr[rng].b

     mcstrct.ew_input[ss,0] = vpstr[iwvI].ew_obs[0] ; wvI
     mcstrct.ew_input[ss,1] = vpstr[iwvI+1].ew_obs[0] ; wvII

  endfor                        ; loop ss=nsys
  
  return, mcstrct
end                             ; sdss_cpvp2mcstrct()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpmc2civstrct, mcstrct, mcindx, instrct=instrct
  ;; Primarily for use with sdss_calcncolmdblt() in
  ;; sdss_completeness_czw
  if n_params() ne 2 then begin
     print,'Syntax - sdss_cpmc2civstrct(mcstrct,mcindx,[indstrct=])'
     return,-1
  endif 

  nstrct = (size(mcstrct,/dim))[0] > 1 ; foil singularity 
  if keyword_set(instrct) then civstrct = instrct $
  else begin
     civstrct = replicate({sdsscivstrct},nstrct)
     ;; Instantiate all the easy stuff 
     civstrct.qso_name = mcstrct.qso_name
     civstrct.mjd = mcstrct.mjd
     civstrct.plate = mcstrct.plate
     civstrct.fiber = mcstrct.fiber
     civstrct.z_qso = mcstrct.z_qso
     civstrct.cflg = mcstrct.cflg
     civstrct.balflg = mcstrct.balflg
  endelse 

  ;; For now have to check that this field exists
  tags = tag_names(mcstrct)
  flgtag = (where(tags eq 'NCOLMFLG_REC'))[0]
  for ii=0,1 do begin
     civstrct.wrest[ii] = mcstrct.wrest[mcindx,ii]
     ;; FINAL = INPUT
     civstrct.wvlim_final[ii,0] = mcstrct.wvlim_input[mcindx,ii,0]
     civstrct.wvlim_final[ii,1] = mcstrct.wvlim_input[mcindx,ii,1]
     civstrct.zabs_final[ii] = mcstrct.zabs_input[mcindx,ii]
     civstrct.ncolm_final[ii] = mcstrct.ncolm_input[mcindx,ii]
     civstrct.ew_final[ii] = mcstrct.ew_input[mcindx,ii]

     ;; ORIG = REC
     civstrct.wvlim_orig[ii,0] = mcstrct.wvlim_rec[mcindx,ii,0]
     civstrct.wvlim_orig[ii,1] = mcstrct.wvlim_rec[mcindx,ii,1]
     civstrct.sigzabs_orig[ii] = mcstrct.sigzabs_rec[mcindx,ii]
     civstrct.ncolm_orig[ii] = mcstrct.ncolm_rec[mcindx,ii]
     civstrct.signcolm_orig[ii] = mcstrct.signcolm_rec[mcindx,ii]
     if flgtag ne -1 then $
        civstrct.ncolmflg_orig[ii] = mcstrct.ncolmflg_rec[mcindx,ii]
     civstrct.ew_orig[ii] = mcstrct.ew_rec[mcindx,ii]
     civstrct.sigew_orig[ii] = mcstrct.sigew_rec[mcindx,ii]
     civstrct.gwidth[ii] = mcstrct.lsnr[mcindx,ii]

     ;; Store flag
     civstrct.ewflg[ii] = mcstrct.flg_rec[mcindx]
  endfor                        ; loop ii=0,1
  ;; Store rating
  civstrct.rating[0] = sdss_getrating(/definite)

  return,civstrct
end                             ; sdss_cpmc2civstrct



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_getqsosubset, snrstrct_fil, out_fil, zbinsize=zbinsize, $
                       list=list, sdsssum=sdsssum, no_write=no_write, $
                       only_list=only_list, index=index,  $
                       snrbinsize=snrbinsize, zmin=zmin, snrmin=snrmin, $
                       dblt_name=dblt_name,mask_fil=mask_fil, clobber=clobber, $
                       _extra=extra
  ;;
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getqsosubset, snrstrct_fil, out_fil, [zbinsize=,'
     print,'                       /list, sdsssum=, /no_write, /only_list, index=,'
     print,'                       snrbinsize=, zmin=, snrmin=, '
     print,'                       dblt_name=, mask_fil=, /clobber, _extra=]'
     return
  endif
  sdssdir = sdss_getsdssdir()
  
  if not keyword_set(zbinsize) then zbinsize = 0.25
  if not keyword_set(snrbinsize) then snrbinsize = 0.25
  if not keyword_set(snrmin) then snrmin = 4
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  if keyword_set(list) then begin
     ;; snrstrct_fil is actually spectra list and should create structure
     if not keyword_set(sdsssum) then $
        stop,'sdss_getqsosubset stop: if input is spectra list, must set sdsssum'

     ;; _extra goes to sdss_measuresnr() and includes dvgal=, dvqso=
     sdss_mksnrtab, snrstrct_fil, sdsssum, snrstrct, _extra=extra
  endif else begin
     if size(snrstrct_fil,/type) eq 7 then $
        snrstrct = xmrdfits(snrstrct_fil,1,/silent) $
     else snrstrct = snrstrct_fil
  endelse 
  tags = tag_names(snrstrct)
  snrtag = (where(tags eq 'SNR_'+strupcase(dblt.ion)))[0]
  
  if keyword_set(zmin) then begin
     ;; Float because that's how the list was selected originally
     gdsnr = where(snrstrct.(snrtag)[2] ge snrmin and $
                   float(snrstrct.z_qso) ge zmin*1d,ngdsnr)
  endif else begin
     gdsnr = where(snrstrct.(snrtag)[2] ge snrmin,ngdsnr)
     if not keyword_set(zmin) then $
        zmin = min(snrstrct.z_qso*100) / 100.
  endelse
  snrstrct = snrstrct[gdsnr]
  plate = strmid(snrstrct.qso_name,6,4)
  snrspec_fil = 'spectro/1d_26/'+plate+'/1d/spSpec-'+snrstrct.qso_name+ '.fit'

  test = file_search(out_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_getqsosubset: out_fil exists and clobber not set; exiting'
     return
  endif 
  
  if keyword_set(mask_fil) then begin
     readcol,mask_fil,mask_spec,format='a',/silent
     mask_spec = mask_spec[1:*]
     nmask_spec = (size(mask_spec,/dim))[0]

     mask4snrspec = intarr(ngdsnr)
     mask4maskspec = intarr(nmask_spec)

     nloop = ngdsnr < nmask_spec ; minimal looping
     for ll=0L,nloop-1 do begin
        if ngdsnr lt nmask_spec then begin
           cur_spec_fil = snrspec_fil[ll]
           isnr = ll
           mtch = where(snrspec_fil[ll] eq mask_spec,nmtch)
           if nmtch gt 1 then $
              stop,'sdss_getqsosubset stop: multiple matches in mask_fil for',$
                   snrspec_fil[ll],nmtch
           if nmtch eq 1 then mtch = mtch[0]
        endif else begin
           cur_spec_fil = mask_spec[ll]
           mtch = where(snrspec_fil eq mask_spec[ll],nmtch)
           if nmtch gt 1 then $
              stop,'sdss_getqsosubset stop: multiple matches in SNR structure for ',$
                   mask_spec[ll],nmtch
           if nmtch eq 1 then begin
              isnr = mtch[0]
              mtch = ll
           endif
        endelse 

        if nmtch eq 0 then $
           print,'sdss_getqsosubset: no matches for ',cur_spec_fil $
        else begin
           mask4snrspec[isnr] = mask4snrspec[isnr]+1   ; exclude
           mask4maskspec[mtch] = mask4maskspec[mtch]+1 ; excluded
        endelse 
     endfor                     ; loop ll=nloop

     ;; Cull
     gdsnr = where(mask4snrspec eq 0,ngdsnr)
     if ngdsnr eq 0 then $
        stop,'sdss_getqsosubset stop: WARNING! Masked out all spectra!'

     snrstrct = snrstrct[gdsnr]
     plate = plate[gdsnr]
     snrspec_fil = snrspec_fil[gdsnr]
  endif                         ; mask_fil set

  ;; Bins will be defined on the left-hand side
  iz = floor((snrstrct.z_qso - zmin) / zbinsize) > 0
  isnr = floor((snrstrct.(snrtag)[2] - snrmin) / snrbinsize)

  ;; Going to make nz by nsnr array accessed by [iz,isnr]
  ;; but linear rank is iz + isnr * nz
  nz = max(iz)                  ; minor rank; rank of unit stride

  index = iz + isnr * nz

  srt = sort(index)
  unq = uniq(index,srt)
  nunq = (size(unq,/dim))[0] > 1

  iz_gd = index[unq] mod nz 
  isnr_gd = (index[unq] - iz_gd) / nz

  srtunq = uniq(index[srt])
  num = (size(srtunq,/dim))[0] > 1
  if num eq 1 then nqso_per_bin = [srtunq[0] + 1] $
  else $
     nqso_per_bin = [srtunq[0]+1,(shift(srtunq,-1)-srtunq)[0:nunq-2]] 
  
  openw,1,out_fil
  printf,1,'abslin/',zmin,snrmin,format='(a7,8x,f7.5,2x,f5.2)'
  writecol,out_fil,snrspec_fil[unq],snrstrct[unq].z_qso,$
           snrstrct[unq].(snrtag)[2],filnum=1,fmt='(a,8x,f7.5,2x,f5.2)'
  close,1
  print,'sdss_getqsosubset: created ',out_fil

  if not keyword_set(only_list) then begin
     ;; _extra= includes /noBAL, /BAL
     root = strmid(out_fil,0,strpos(out_fil,'.',/reverse_searcH))
     sdsssum = root+'.fit'
     snrstrct_fil = root+'_SNR.fit'
     sdss_getqsoinlist,out_fil,sdsssum,snrstrct_fil,_extra=extra
  endif                         ; /only_list

  if not keyword_set(no_write) then begin
     ;; Random information
     ofil = 'sdss_getqsosubset.tab'
     openw,1,ofil
     printf,1,'z','S/N','No.',format='(a6,2x,a5,2x,a6)'
     writecol,ofil,iz_gd*zbinsize+zmin,isnr_gd*snrbinsize+snrmin,nqso_per_bin,$
              fmt='(f6.4,2x,f5.2,2x,i6)',filnum=1
     close,1
     print,'sdss_getqsosubset: created ',ofil

;    x_splot,snrstrct[unq].z_qso,snrstrct[unq].(snrtag)[2],psym1=4,/block
  endif                         ; /no_write

end                             ; sdss_getqsosubset



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_iterqsosubset, snrstrct_fil, out_fil, niter, only_list=only_list,$
                        mask_fil=mask_fil, clobber=clobber, _extra=extra
  if n_params() ne 3 then begin
     print,'Syntax - sdss_iterqsosubset, snrstrct_fil, out_fil, niter, [/only_list,'
     print,'                        mask_fil=, /clobber, _extra=]'
     return
  endif

  test = file_search(out_fil+'*', count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_iterqsosubset: file exists; will not clobber ',out_fil
     return
  endif
  
  ;; First pass
  ;; _extra= includes zbinsize=, /list, sdsssum=, index=, snrbinsize=,
  ;; zmin=, snrmin=, dblt_name=, [_extra=]
  print,'sdss_iterqsosubset: iter = ',1
  sdss_getqsosubset, snrstrct_fil, out_fil, /only_list, mask_fil=mask_fil, $
                     /no_write, /clobber, _extra=extra
  openr,17,out_fil
  dum = ''
  readf,17,dum   
  all_fil = dum                 ; first line
  while not eof(17) do begin
     readf,17,dum
     all_fil = [all_fil, dum] 
  endwhile
  close,17

  for ii=1,niter-1 do begin
     print,'sdss_iterqsosubset: iter = ',ii+1

     if keyword_set(mask_fil) then begin
        ;; Read in previously made file and append to mask_fil
        if ii eq 1 then begin
           spawn,'cp '+mask_fil+' '+mask_fil+'.orig'
           print,'sdss_iterqsosubset: created backup ',mask_fil+'.orig'
        endif

        readcol, out_fil, file, skip=1, format='a', /silent
        openw, 1, mask_fil, /append
        writecol, mask_fil, file, fmt='(a)', filnum=1
        close,1
        print,'sdss_iterqsosubset: updated ',mask_fil
     endif

     sdss_getqsosubset, snrstrct_fil, out_fil, /only_list, $
                        mask_fil=mask_fil, /no_write, /clobber, _extra=extra
  
     ;; Would be better to just track numbers instead of reading twice
     openr,17,out_fil
     dum = ''
     readf,17,dum               ; skip first line
     while not eof(17) do begin
        readf,17,dum
        all_fil = [all_fil, dum] 
     endwhile
     close,17
     
  endfor                        ; loop ii=niter
  
  if keyword_set(mask_fil) then begin
     ;; Read in previously made file and append to mask_fil
     readcol, out_fil, file, skip=1, format='a', /silent
     openw, 1, mask_fil, /append
     writecol, mask_fil, file, fmt='(a)', filnum=1
     close,1
     print,'sdss_iterqsosubset: updated ',mask_fil
  endif

  ;; Write files
  writecol, out_fil, all_fil, fmt='(a)'
  print,'sdss_iterqsosubset: created concatenated ',out_fil
  
  ;; _extra= includes /noBAL, /BAL
  root = strmid(out_fil,0,strpos(out_fil,'.',/reverse_searcH))
  sdsssum = root+'.fit'
  snrstrct_fil = root+'_SNR.fit'
  sdss_getqsoinlist,out_fil,sdsssum,snrstrct_fil,_extra=extra

  
end                             ; sdssS_iterqsosubset



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcfracpix, limits, global_arr, binsize=binsize, npix=npix, $
                           index_lim=index_lim, limits2=limits2, $
                           array2=array2, global_arr2=global_arr2, $
                           dglobal_arr2=dglobal_arr2
  ;; Largely needed for completeness tests
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcfracpix(limits, global_arr, [binsize=, npix=,'
     print,'                          index_lim=, limits2=, array2=, '
     print,'                          global_arr2=, dglobal_arr2=])'
     return,-1
  endif 

  if not keyword_set(binsize) then binsize = global_arr[1] - global_arr[0]
  if keyword_set(global_arr2) then begin
     if not keyword_set(limits2) then $
        stop,'sdss_calcfracpix() stop: global_arr2 set but not limits2'
     if not keyword_set(dglobal_arr2) then begin
        ;; this will be a MONSTER fail for co-moving pathlength X, for
        ;; which this secondary option is likely to work with
        dglobal_arr2 = shift(global_arr2,-1) - global_arr2
        nbin = (size(global_arr2,/dim))[0] > 1
        dglobal_arr2[nbin-1] = dglobal_arr2[nbin-2] ; better this way
     endif 
  endif 
  
  index_lim = floor((limits - global_arr[0]) / binsize) 
  index_lim = reform(index_lim) ; collapse dimensions of 1
  npix = index_lim[1] - index_lim[0] + 1
  case npix of
     1: begin
        ;; Fit all inside one z grid
        array = limits[1] - limits[0]
        if keyword_set(global_arr2) then $
           array2 = limits2[1] - limits2[0]
     end                        ; 1
     2: begin
        ;; Spans a border
        array = [global_arr[index_lim[1]]-limits[0], $
                 limits[1]-global_arr[index_lim[1]]]
        if keyword_set(global_arr2) then $
           array2 = [global_arr2[index_lim[1]]-limits2[0], $
                     limits2[1]-global_arr2[index_lim[1]]]
     end                        ; 2
     else: begin
        ;; Handle two sides and full bins in the middle
        array = replicate(binsize, npix)
        array[[0,npix-1]] = [global_arr[index_lim[0]+1]-limits[0], $
                             limits[1]-global_arr[index_lim[1]]]
        if keyword_set(global_arr2) then begin
           array2 = dglobal_arr2[index_lim[0]:index_lim[1]]
           array2[[0,npix-1]] = [global_arr2[index_lim[0]+1]-limits2[0], $
                                 limits2[1]-global_arr2[index_lim[1]]]
        endif 
     end
  endcase                       ; npix

  return, array                 ; dimensions of index_lim
end                             ; sdss_calcfracpix()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcdzlos, z_qso, dblt_name, zlim, z_global, spec_fil, rz_los=rz_los,$
                         zlim_los=zlim_los, zbinsize=zbinsize, dvem=dvem, $
                         nzpix=nzpix, iz_los=iz_los, dx_los=dx_los, $
                         x_global=x_global, dx_global=dx_global, $
                         xlim_los=xlim_los, wave=wave, flux=flux, error=error, $
                         debug=debug, _extra=extra
  ;; Returns: rz_los=, zlim_los=, nzpix=, iz_los=, dx_los=, xlim_los=
  ;; _extra includes cosmology=
  if n_params() lt 4 then begin
     print,'Syntax - sdss_calcdzlos( z_qso, dblt_name, zlim, z_global, [spec_fil, rz_los=,'
     print,'                         zlim_los=, zbinsize=, nzpix=, iz_los=, dx_los=, '
     print,'                         x_global=, dx_global=, xlim_los=, wave=, flux=, error=,'
     print,'                         /debug, _extra=])'
     return, -1
  endif 
  
  cosmology = sdss_setcosmology(_extra=extra) 
  skylinwv = sdss_getskylinwave(dwv=dwv) ; Angstrom cut 
  cinv = 1./299792.458d

  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  ;; 5579 +/- 5 Ang for CIV 1548 is 2.60032 to 2.60678 in z
  ;;                    MgII 2796 is 0.99331 to 0.99689 in z
  ;; 6302 +/- 5 Ang for CIV 1548 is 3.06732 to 3.07378 in z
  ;;                    MgII 2796 is 1.25186 to 1.25545 in z
  zsky_lim1 = (skylinwv[0] + [-1.,1.]*dwv[0])/ dblt.wvI - 1.
  zsky_lim2 = (skylinwv[1] + [-1.,1.]*dwv[1])/ dblt.wvI - 1.
  zsky_lim = transpose([[zsky_lim1],[zsky_lim2]])

  ;; Care about fractional pixels
  ;; _extra includes dvqso=, dvgal=, nsig=, and more
  snr = sdss_measuresnr('blank',wvlim_obs=wvlim_los,/no_snr,$ 
                        dblt_name=dblt,zqso=z_qso,_extra=extra)

  zlim_los = wvlim_los / dblt.wvI - 1. ; returned value

  if zlim_los[1] lt zlim[0] or zlim_los[0] gt zlim[1] or $
     zlim_los[1] lt zlim_los[0] then begin
     nzpix = 0
     if zlim_los[1] lt zlim_los[0] then zlim_los[1] = zlim_los[0]
     return, -1                 ; EXIT
  endif 

  ;; Force range useful range
  zlim_los[0] = zlim_los[0] > zlim[0]
  zlim_los[1] = zlim_los[1] < zlim[1]
  zlim_los0 = zlim_los          ; for gaps and emission lines

  done = 0
  if keyword_set(debug) then print,''

  ;; Fractional contribution but account for sky lines
  if (zlim_los[1] lt zsky_lim[0,0] $                                      ; all below
      or zlim_los[0] gt zsky_lim[1,1]) $                                  ; all above
     or (zlim_los[0] gt zsky_lim[0,1] and zlim_los[1] lt zsky_lim[1,0]) $ ; all inbtw
  then begin
     ;; ;;;;;;;
     ;; Quick return righ there because no funny business necessary
     done = 1
     goto, skip_to_end

  endif else if keyword_set(debug) then $
     print,'sdss_calcdzlos() debug: have to modify zlim_los ',$
           zlim_los[0],zlim_los[1]

  ;; ;;;;;;;
  ;; There are 7 flavors of overlap to check and modify zlim_los
  ;; 7) LOS completely inside one sky line (2 setups) 
  if (zlim_los[0] ge zsky_lim[0,0] and zlim_los[1] le zsky_lim[0,1]) or $
     (zlim_los[0] ge zsky_lim[1,0] and zlim_los[1] le zsky_lim[1,1]) then begin
     nzpix = 0 
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: all inside one sky line.'
     return, -1                 ; EXIT
  endif 

  ;; 1) LOS encloses both sky lines completely (1 setup) 
  if zlim_los[0] lt zsky_lim[0,0] and zlim_los[1] gt zsky_lim[1,1] then begin
     zlim_los = [[zlim_los[0],zsky_lim[0,0]],$
                 [zsky_lim[0,1],zsky_lim[1,0]],$
                 [zsky_lim[1,1],zlim_los[1]]]
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: completely enclose both sky lines.'
     goto, skip_to_end          ; I know this is ugly
  endif 

  ;; 2) LOS begins in first sky line and ends in second (1) 
  if zlim_los[0] ge zsky_lim[0,0] and zlim_los[0] le zsky_lim[0,1] and $
     zlim_los[1] ge zsky_lim[1,0] and zlim_los[1] le zsky_lim[1,1] then begin
     zlim_los = [zsky_lim[0,1], zsky_lim[1,0]]
     done = 1                   ; EASY
     goto, skip_to_end
  endif 

  ;; 3) LOS begins or ends in one sky line and encloses the other (2) 
  ;; 3a) Begins in first and encloses second
  if (zlim_los[0] ge zsky_lim[0,0] and zlim_los[0] le zsky_lim[0,1] and $
      zlim_los[1] gt zsky_lim[1,1]) then begin
     zlim_los = [[zsky_lim[0,1],zsky_lim[1,0]],$
                 [zsky_lim[1,1],zlim_los[1]]]
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: begins in first and encloses second'
     goto, skip_to_end
  endif 
  ;; 3b)
  if (zlim_los[0] lt zsky_lim[0,0] and $
      zlim_los[1] ge zsky_lim[1,0] and zlim_los[1] le zsky_lim[1,1]) then begin
     zlim_los = [[zlim_los[0],zsky_lim[0,0]],$
                 [zsky_lim[0,1],zsky_lim[1,0]]]
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: encloses first and ends in second'
     goto, skip_to_end
  endif 

  ;; 4) LOS encloses one sky line completely (2) 
  ;; 4a) First one enclosed and end before second
  if zlim_los[0] lt zsky_lim[0,0] and zlim_los[1] gt zsky_lim[0,1] and $
     zlim_los[1] lt zsky_lim[1,0] then begin
     zlim_los = [[zlim_los[0], zsky_lim[0,0]],$
                 [zsky_lim[0,1],zlim_los[1]]]
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: completely enclose first sky line only.'
     goto, skip_to_end
  endif 
  ;; 4b) Second one enclosed and begins after first
  if zlim_los[0] gt zsky_lim[0,1] and $
     zlim_los[0] lt zsky_lim[1,0] and zlim_los[1] gt zsky_lim[1,1] then begin
     zlim_los = [[zlim_los[0], zsky_lim[1,0]], $
                 [zsky_lim[1,1],zlim_los[1]]]
     if keyword_set(debug) then $
        print,'sdss_calcdzlos() debug: completely enclose second sky line only.'
     goto, skip_to_end
  endif 

  ;; 5) LOS begins in one sky line and does not enclose other (2) 
  ;; 5a) Begins in first and ends before second; already checked
  ;; not completely enclosed
  if zlim_los[0] ge zsky_lim[0,0] and zlim_los[0] le zsky_lim[0,1] and $
     zlim_los[1] lt zsky_lim[1,0] then begin
     zlim_los[0] = zsky_lim[0,1]
     done = 1
     goto, skip_to_end
  endif 
  ;; 5b) Begins in second 
  if zlim_los[0] ge zsky_lim[1,0] and zlim_los[0] le zsky_lim[1,1] then begin
     zlim_los[0] = zsky_lim[1,1]
     done = 1
     goto, skip_to_end
  endif 

  ;; 6) LOS ends in one sky line and does not enclose other (2) 
  ;; 6a) Ends in first; Already checked not completely enclosed
  if zlim_los[1] ge zsky_lim[0,0] and zlim_los[1] le zsky_lim[0,1] then begin
     zlim_los[1] = zsky_lim[0,0]
     done = 1
     goto, skip_to_end
  endif 
  ;; 6b) Ends in second but begins before first
  if zlim_los[0] gt zsky_lim[0,1] and $
     zlim_los[1] ge zsky_lim[1,0] and zlim_los[1] le zsky_lim[1,1] then begin
     zlim_los[1] = zsky_lim[1,0]
     done = 1
     goto, skip_to_end
  endif 


  ;; Now process results
  skip_to_end:        

  ;; Adjust sdss_calcdzlos() results based on actual spectrum due to: 
  ;; i) gaps or 
  ;; ii) emission lines or other fixed regions
  if keyword_set(dvem) then begin
     ;; Emission-line rest wavelength; dv range
     ;; either [wrest, dvlo, dvhi] or 
     ;; transpose([[wrest1,dvlo1,dvhi1],[wrest2,dvlo2,dvhi2]])
     n_em = size(dvem,/n_dim)

     ;; Have to insert this in to zlim_los
     if size(zlim_los,/n_dim) gt 1 then $  ; could have gap
        nbin = n_elements(zlim_los[0,*]) $ ; not transposed yet
     else nbin = 1
     if keyword_set(debug) then print,'sdss_calcdzlos() debug: dvem=; nbin = ',nbin

     tmp_zlim_los = dblarr(2,nbin+2*n_em,/nozero)   ; not transposed yet
     tmp_zlim_los[*,0:nbin-1] = zlim_los            ; zlim_los n_dim won't matter
     if n_em eq 1 then dvem_use = transpose(dvem) $ ; [1,3]
     else dvem_use = dvem                           ; [n_em,3]
     for ii=0,n_em-1 do begin
        zem_dblt = dvem_use[ii,0]/dblt.wvI*(1+z_qso)
        tmp_zlim_los[0,2*ii+nbin] = zlim_los0[0]
        tmp_zlim_los[1,2*ii+nbin] = dvem_use[ii,1]*cinv * (1 + zem_dblt) + zem_dblt   ; em. line begin
        tmp_zlim_los[0,2*ii+1+nbin] = dvem_use[ii,2]*cinv * (1 + zem_dblt) + zem_dblt ; em. line end
        tmp_zlim_los[1,2*ii+1+nbin] = zlim_los0[1]
     endfor                                                                           ; loop ii=n_em
     
     ;; Stay within bounds
     gd = where(tmp_zlim_los[0,nbin:*] lt tmp_zlim_los[1,nbin:*],n_em) ; overwrite n_em!
     if n_em gt 0 then begin
        tmp = dblarr(2,nbin+n_em,/nozero)
        tmp[*,0:nbin-1] = tmp_zlim_los[*,0:nbin-1] ;  orig
        tmp[*,nbin:*] = tmp_zlim_los[*,nbin+gd]    ; emission lines
        tmp_zlim_los = tmp

        ;; Replicate code from keyword_set(spec_fil) or keyword_set(wave) below
        srt = sort(tmp_zlim_los[0,*]) ; going to look for gaps again
        tmp_zlim_los[0,*] = tmp_zlim_los[0,srt] 
        tmp_zlim_los[1,*] = tmp_zlim_los[1,srt] 
        srt = sort(tmp_zlim_los[1,*]) ; makes a difference!

        ;; Brute force (but hopefully short) loop
        new_zlim_los = replicate(-1d,2,nbin+n_em)
        count = 0 
        new_zlim_los[*,count] = tmp_zlim_los[*,0] ; very beginning
        while srt[0] ne -1 do begin
           ;; take next highest upper limit and place it where count is
           gd = where(tmp_zlim_los[1,srt] gt new_zlim_los[0,count],ngd,$
                      complement=bd,ncomplement=nbd)
;           if ngd ne 0 then begin
              new_zlim_los[1,count] = tmp_zlim_los[1,srt[gd[0]]]
              ;; Next lowest lower limit above last upper limit at count++
              gd = where(tmp_zlim_los[0,*] gt new_zlim_los[1,count],ngd)
              if ngd ne 0 and count ne nbin+n_em-1 then begin
                 count++
                 new_zlim_los[0,count] = tmp_zlim_los[0,gd[0]] 
              endif else srt = -1 ; nothing higher than last upper
;           endif else srt = -1    ; nothing greater than last lower
        endwhile                  ; srt[0] eq -1 

        if keyword_set(debug) then begin
           print,'sdss_calcdzlos() debug: dvem=; Old and new LOS z limits, excluding emission line(s)'
           printcol,tmp_zlim_los[0,*],tmp_zlim_los[1,*],$
                    new_zlim_los[0,*],new_zlim_los[1,*],$
                    format='(2(f8.5,1x),1x,"|",1x,2(f8.5,1x))'
        endif 

        ;; Clean up
        if count gt 0 then $
           zlim_los = new_zlim_los[*,0:count]
        if size(zlim_los,/n_dim) gt 1 then done = 0 $ ; must transpose
        else done = 1 
     endif                      ; n_em > 0
  endif                         ; dvem=

  if keyword_set(spec_fil) or keyword_set(wave) or keyword_set(dvem) then begin
     if keyword_set(spec_fil) then $
        parse_sdss,spec_fil,flux,wave,sig=error ; overrides input wave=, flux=, error=
     redshift = wave/dblt.wvI - 1.
     
     ;; Could also add in search for observed wavelength within certain
     ;; given limtis (e.g., CIV emission line) to this search
     badpix = where(error le 0. and redshift ge min(zlim_los,max=mx) and $
                    redshift le mx,nbadpix,complement=goodpix,ncomplement=ngoodpix) ; hope SDSS and sdss_calcnormerr() did this right
     if nbadpix ne 0 then begin
        if keyword_set(debug) then zlim_los0_safe = zlim_los0
        zlim_los0[0] = zlim_los0[0] > redshift[goodpix[0]] ; won't need zlim_lim0
        zlim_los0[1] = zlim_los0[1] < redshift[goodpix[ngoodpix-1]]
        if keyword_set(debug) then $
           print,'sdss_calcdzlos() debug; spec_fil or wave=; Old and new zlim_los0 = ',$
                 string(zlim_los0_safe,zlim_los0,format='(2(2(f9.4,1x),1x))')
        
        istart = badpix[where(badpix ne shift(badpix,1)+1,ngaps)]
        istop = badpix[where(badpix ne shift(badpix,-1)-1,ntest)]
        
        if ngaps ne ntest then $
           stop,'sdss_calcdzlos() stop: # istart != # istop (should not be here)'

        ;; Have to insert this in to zlim_los
        if size(zlim_los,/n_dim) gt 1 then $  ; could have gap
           nbin = n_elements(zlim_los[0,*]) $ ; not transposed yet
        else nbin = 1
        if keyword_set(debug) then print,'sdss_calcdzlos() debug: spec_fil or wave=; nbin = ',nbin

        tmp_zlim_los = dblarr(2,nbin+ngaps+1,/nozero)
        tmp_zlim_los[*,0:nbin-1] = zlim_los ; zlim_los n_dim won't matter
        tmp_zlim_los[*,nbin] = [zlim_los0[0],redshift[istart[0]]]
        for ii=1,ngaps-1 do $
           tmp_zlim_los[*,ii+nbin] = [redshift[istop[ii-1]],redshift[istart[ii]]]
        tmp_zlim_los[*,nbin+ngaps] = [redshift[istop[ngaps-1]],zlim_los0[1]]

        ;; Stay within bounds
        gd = where(tmp_zlim_los[0,nbin:*] lt tmp_zlim_los[1,nbin:*],ngaps) ; overwrite ngaps!
        if ngaps gt 0 then begin
           tmp = dblarr(2,nbin+ngaps,/nozero)
           tmp[*,0:nbin-1] = tmp_zlim_los[*,0:nbin-1] ;  orig
           tmp[*,nbin:*] = tmp_zlim_los[*,nbin+gd]    ; emission lines
           tmp_zlim_los = tmp

           ;; Replicate code from keyword_set(dvem) above
           srt = sort(tmp_zlim_los[0,*]) ; going to look for gaps again
           tmp_zlim_los[0,*] = tmp_zlim_los[0,srt] 
           tmp_zlim_los[1,*] = tmp_zlim_los[1,srt] 
           srt = sort(tmp_zlim_los[1,*]) ; makes a difference!

           ;; Brute force (but hopefully short) loop
           new_zlim_los = replicate(-1d,2,nbin+ngaps)
           count = 0 
           new_zlim_los[*,count] = tmp_zlim_los[*,0] ; very beginning
           while srt[0] ne -1 do begin
              ;; take next highest upper limit and place it where count is
              gd = where(tmp_zlim_los[1,srt] gt new_zlim_los[0,count],ngd,$
                         complement=bd,ncomplement=nbd)
;              if ngd ne 0 then begin
                 new_zlim_los[1,count] = tmp_zlim_los[1,srt[gd[0]]]
                 ;; Next lowest lower limit above last upper limit at count++
                 gd = where(tmp_zlim_los[0,*] gt new_zlim_los[1,count],ngd)
                 if ngd ne 0 and count ne nbin+ngaps-1 then begin
                    count++
                    new_zlim_los[0,count] = tmp_zlim_los[0,gd[0]] 
                 endif else srt = -1 ; nothing higher than last upper
;              endif else srt = -1    ; nothing greater than last lower
           endwhile                  ; srt[0] eq -1 

           if keyword_set(debug) then begin
              print,'sdss_calcdzlos() debug: spec_fil or wave=; Old and new LOS z limits, excluding gaps'
              printcol,tmp_zlim_los[0,*],tmp_zlim_los[1,*],$
                       new_zlim_los[0,*],new_zlim_los[1,*],$
                       format='(2(f8.5,1x),1x,"|",1x,2(f8.5,1x))'
           endif 

           ;; Clean up
           if count eq 0 then zlim_los = new_zlim_los[*,0:count] 
           if size(zlim_los,/n_dim) gt 1 then done = 0 $ ; must transpose
           else done = 1 
        endif                   ; ngaps > 0
     endif                      ; nbadpix != 0
  endif                         ; spec_fil or wave= set

  
  ;; Calculate pathlength
  if min(zlim_los,max=mx) lt zlim[0] or mx gt zlim[1] then $
     stop,'sdss_calcdzlos() stop: zlim_los spills outside input zlim (should not get here)'
  if done eq 1 then begin
     ;; no need to transpose zlim_los
     xlim_los = cosm_xz(zlim_los,/silent,/exact,/noinit)
     dz_los = sdss_calcfracpix(zlim_los, z_global, binsize=zbinsize, $
                               npix=nzpix, index_lim=iz_los, array2=dx_los, $
                               limits2=xlim_los, global_arr2=x_global,$
                               dglobal_arr2=dx_global)
     rz_los = replicate(1,nzpix)
  endif else begin
     bd = where(zlim_los[*,0] gt zlim_los[*,1])
     if bd[0] ne -1 then $
        stop, 'sdss_calcdzlos() stop: >= 1 zlim_los[*,0] > zlim_los[*,1] (should not get here)'

     ;; 1) (3), 3) , 4) 
     ;; must have zlim_los transposed
     zlim_los = transpose(zlim_los) 
     nzlim = (size(zlim_los,/dim))[0]  
     xlim_los = zlim_los
     for zz=0,nzlim-1 do begin
        xlim_los[zz,*] = cosm_xz(zlim_los[zz,*],/silent,/exact,/noinit)
        if zz gt 0 then $
           if zlim_los[zz,0] lt zlim_los[zz-1,1] then $
              stop,'sdss_calcdzlos() stop: zlim_los[zz,0] < zlim_los[zz-1,1] (should not get here)'
        dz_tmp = sdss_calcfracpix(zlim_los[zz,*], z_global, binsize=zbinsize, $
                                  npix=nztmp, index_lim=iz_tmp, array2=dx_tmp, $
                                  limits2=xlim_los[zz,*], global_arr2=x_global,$
                                  dglobal_arr2=dx_global)
        bd = where(dz_tmp lt 0 or dx_tmp lt 0)
        if bd[0] ne -1 then $
           stop,'sdss_calcdzlos() stop: sdss_calcfracpix() returned dz < 0 (should not get here)'
        rz_tmp = replicate(1,nztmp)
        if zz eq 0 then begin
           ;; Instantiate
           dz_los = dz_tmp 
           nzpix = nztmp
           iz_los = iz_tmp
           dx_los = dx_tmp
           rz_los = rz_tmp
        endif else begin
           ;; Save but take notice of the "seams"
           if iz_los[1] eq iz_tmp[0] then begin
              ;; The sky line gap falls in same bin so sum
              ;; *assumes* zlim_los[] is in numerical order
              dz_los[nzpix-1] = dz_los[nzpix-1] + dz_tmp[0]
              dx_los[nzpix-1] = dx_los[nzpix-1] + dx_tmp[0]
              if nztmp gt 1 then begin
                 dz_los = [dz_los, dz_tmp[1:*]]
                 dx_los = [dx_los, dx_tmp[1:*]]
                 rz_los = [rz_los, rz_tmp[1:*]]
              endif 
              nzpix = nzpix + nztmp - 1
           endif else begin
              ;; Have some gaps to fill in
              nmiss = iz_tmp[0] - iz_los[1] - 1
              if nmiss gt 0 then begin
                 fill = fltarr(nmiss) ; 0
;                 stop,'there is a gap'
                 dz_los = [dz_los, fill, dz_tmp]
                 dx_los = [dx_los, fill, dx_tmp]
                 rz_los = [rz_los, fill, rz_tmp]
              endif else begin
;                 stop,'there is not a gap but a continuation?'
                 dz_los = [dz_los, dz_tmp]
                 dx_los = [dx_los, dx_tmp]
                 rz_los = [rz_los, rz_tmp]                 
              endelse 
              nzpix = nzpix + nmiss + nztmp
           endelse 
           iz_los[1] = iz_tmp[1]

        endelse                 ; zz > 0
     endfor                     ; loop zz=nzlim

     if zlim_los[nzlim-1,1] lt zlim_los[0,0] then stop
  endelse                 
  
  return, dz_los

end                             ; sdss_calcdzlos()




;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_fndqsopairs, sdsssum, nclose=nclose, _extra=extra
  ;; Find the sightlines close to each other
  ;; Maybe eventually evolve to look at correlated absorbers
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fndqsopairs( sdsssum, [nclose=, _extra=])'
     return,-1
  endif 
  if not keyword_set(nclose) then nclose = 5
  cosmology = sdss_setcosmology(_extra=extra) ; _extra= cosmology=

  if size(sdsssum,/type) eq 7 then $
     sdsstab = xmrdfits(sdsssum,1,/silent) $
  else sdsstab = sdsssum
  nqso = (size(sdsstab,/dim))[0] > 1 ; foil singularity
  
  ;; Get names
  spec_fil = sdss_getname(sdsstab,root=qso_name)

  strct = {qso_name1:qso_name, $
           z_qso1:sdsstab.z, $
           qso_name2:strarr(nqso,nclose), $ ; closest neighbors
           z_qso2:dblarr(nqso,nclose,/nozero), $
           das:fltarr(nqso,nclose,/nozero), $
           dMpc:fltarr(nqso,nclose,/nozero), $ ; which ever is foreground
           dvqso:fltarr(nqso,nclose,/nozero) $ ; (z2-z1)/(1+z1)*c
          }

  deg_hr = 24/360.
  as_am = 1./60                 ; arcsec per arcmin
  c = 299792.458                ; km/s
  index = lindgen(nqso)
  for qq=0L,nqso-1 do begin
     case qq of 
        0: rng = index[qq+1:nqso-1]
        nqso-1: rng = index[0:nqso-2]
        else: rng = [index[0:qq-1],index[qq+1:*]]
     endcase

     ;; Calculate distance in arcsec
     gcirc,1,sdsstab[qq].ra*deg_hr,sdsstab[qq].dec,$
           sdsstab[rng].ra*deg_hr,sdsstab[rng].dec,dist
     srt = (sort(dist))[0:nclose-1]
     
     ;; Load up structure
     strct.qso_name2[qq,*] = qso_name[rng[srt]]
     strct.das[qq,*] = dist[srt]
     strct.z_qso2[qq,*] = sdsstab[rng[srt]].z

     ;; Distance in kpc for whichever is foreground
     fg1 = where(strct.z_qso1[qq] lt strct.z_qso2[qq,*],complement=fg2)
     if fg1[0] ne -1 then begin
        ;; Primary object is in foreground
        ;; amin is Mpc per arcmin
        lstar = x_calclstar(strct.z_qso1[qq],amin=amin,om=cosmology[1],$
                            h0=cosmology[0],ov=cosmology[2],/silent)
        strct.dMpc[qq,fg1] = strct.das[qq,fg1] * as_am * amin
;        strct.dvqso[qq,fg1] = c*(strct.z_qso2[qq,fg1]-strct.z_qso1[qq])/$
;                              (1+strct.z_qso1[qq])
     endif 
     if fg2[0] ne -1 then begin
        ;; Secondary object is in foreground
        lstar = x_calclstar(strct.z_qso1[qq],amin=amin,om=cosmology[1],$
                            h0=cosmology[0],ov=cosmology[2],/silent)
        strct.dMpc[qq,fg2] = strct.das[qq,fg2] * as_am * amin
;        strct.dvqso[qq,fg2] = c*(strct.z_qso1[qq]-strct.z_qso2[qq,fg2])/$
;                              (1+strct.z_qso2[qq,fg2])
     endif
     strct.dvqso[qq,*] = c*(strct.z_qso2[qq,*]-strct.z_qso1[qq])/$
                         (1+strct.z_qso1[qq])

  endfor                        ; loop qq=nqso
  
  return,strct
end                             ; sdss_fndqsopairs()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_mkzarr, zlim, zbinsize, header=header, nzbin=nzbin
  ;; consistent way of (re)constructing the global redshift array (binned)
  if n_params() ne 2 and not keyword_set(header) then begin
     print,'Syntax - sdss_mkzarr( zlim, zbinsize, [header=, nzbin=])'
     return, -1
  endif 

  if keyword_set(header) then begin
     ;; Recreate necessary values
     zmin = sxpar(header,'ZMIN')
     zmax = sxpar(header,'ZMAX')
     zlim = [zmin,zmax]
     
     zbinsize = sxpar(header,'ZBINSZ')

     nzbin = sxpar(header,'NAXIS1')
  endif else nzbin = ceil((zlim[1] - zlim[0])/zbinsize) ; + 1L ; why +1?
  zglobal = zlim[0] + findgen(nzbin)*zbinsize
  
  return, zglobal

end                             ; sdss_mkzarr()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_mkewarr, ewlim, ewbinsize, header=header, newbin=newbin
  ;; consistent way of (re)constructing the global EWlim array (binned)
  if n_params() ne 2 and not keyword_set(header) then begin
     print,'Syntax - sdss_mkewarr( ewlim, ewbinsize, [header=, newbin=])'
     return, -1
  endif 

  if keyword_set(header) then begin
     ;; Recreate necessary values
     ewmin = sxpar(header,'EWMIN')
     ewmax = sxpar(header,'EWMAX')
     ewlim = [ewmin,ewmax]
     
     ewbinsize = sxpar(header,'EWBINSZ')

     newbin0 = sxpar(header,'NAXIS2')
  endif 

  ewglobal = sdss_mkzarr(ewlim, ewbinsize, nzbin=newbin)
  if keyword_set(newbin0) then $
     if newbin0 ne newbin then $
        stop,'sdss_mkewarr() stop: header and mkzarr() return different newbin',$
  newbin0, newbin
  
  return, ewglobal

end                             ; sdss_mkewarr()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcsigpoiss, number, sigma=sigma, cl=cl, verbose=verbose
  ;; Re-hash x_poisscl() to handle arrays more efficiently and to
  ;; actually take the difference
  if n_params() ne 1 then begin
     print,'Syntax - sdss_calcsigpoiss(number, [cl=, sigma=, /verbose])'
     return, -1
  endif 

  if keyword_set(sigma) then begin
     ;; DRM Change: Should be double-sided CL passed to x_poisscl()
     cl = 2.*gauss_pdf(sigma) - 1 ; double-sided CL corresponding to gaussian sigma
     if keyword_set(verbose) then $
        print,'sdss_calcsigpoiss(): Given sigma = '+strtrim(sigma,2)+$
              '; using CL = ',cl
  endif 
  if not keyword_set(cl) then begin
     ;; DRM Change: Should be double-sided CL passed to x_poisscl()
     cl = 2.*gauss_pdf(1.)-1.   ; = .6827  (ie, 1-sigma)
  endif
  
  num = (size(number,/dim))[0] > 1 ; foil singularity
  sig = fltarr(num,2,/nozero)      ; temporary
  
  ;; x_poisscl() can't handle arrays and will make a
  ;; Gaussian approx if the observed number > 120 so just do the
  ;; array operation here and loop over the rest
  ipoiss = where(number le 120,complement=igauss)
  if igauss[0] ne -1 then begin
     ;; DRM Change: Transform to single-sided CL
     icl = cl + (1.-cl)/2.      ; single-sided CL
     S = abs(gauss_cvf(icl))
     high = (number[igauss]+1.) * $
            ( 1. - 1./(9*(number[igauss]+1)) + $
              S/(3*sqrt(number[igauss]+1.)) )^3
     low = number[igauss] - S*sqrt(number[igauss]) + $
           (S^2 - 1)/3.
     sig[igauss,0] = number[igauss] - low
     sig[igauss,1] = high - number[igauss] 
  endif                         ; igauss[0] ne -1
  
  
  if ipoiss[0] ne -1 then begin
     ;; Minimize looping
     ipoiss = ipoiss[sort(number[ipoiss])]
     unq = uniq(number[ipoiss])
     nunq = (size(unq,/dim))[0] > 1
     for jj=0L,nunq-1 do begin
        ;; _extra includes nspl= and maxx= for number of spline
        ;; points and max point for making numbers; may have to
        ;; play with this
        p = x_poisscl(number[ipoiss[unq[jj]]], cl, /silent)
        if jj eq 0 then rng = lindgen(unq[jj]+1) $
        else rng = lindgen(unq[jj]-unq[jj-1]) + unq[jj-1] + 1
        
        sig[ipoiss[rng],0] = number[ipoiss[unq[jj]]] - p[1]
        sig[ipoiss[rng],1] = p[0] - number[ipoiss[unq[jj]]]
     endfor                     ; loop jj=nuq
     
  endif                         ; ipoiss[0] ne -1
        
  sig = reform(sig)             ; remove dimen 1
  return, sig
end                             ; sdss_calcsigpoiss()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcsigbinom, ninput_2darr, nrec_2darr, sigma=sigma, $
                            cl=cl, verbose=verbose
  ;; Binomial error esimate of completeness
  ;; C(z,W) proportional to nrec_2darr / ninput_2darr
  ;; Wilson score interval:
  ;; Lower: (p + 1/(2n) * zscore^2 - zscore * sqrt( p*(1-p)/n +
  ;;                                 zscore^2/(4n^2))) / (1 +
  ;;                                 1/n*zscore^2)
  ;; Upper: (p + 1/(2n) * zscore^2 + zscore * sqrt( p*(1-p)/n -
  ;;                                 zscore^2/(4n^2))) / (1 +
  ;;                                 1/n*zscore^2)
  ;; where p = Nrec / Ninput (estimate of probability)
  ;; and n = number of trials (Ninput)
  ;; and zscore = abs(gauss_cvf(1-0.5*(1-cl)))
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcsigbinom(ninput_2darr, nrec_2darr, [cl=, sigma=, /verbose])'
     return, -1
  endif 

  if keyword_set(sigma) then begin
     ;; DRM Change: Should be double-sided CL to start, according to
     ;; the definition of ZSCORE used below.
     cl = 2.*gauss_pdf(sigma) - 1. ; double-sided CL
     if keyword_set(verbose) then $
        print,'sdss_calcsigbinom(): Given sigma = '+strtrim(sigma,2)+$
              '; using CL = ',cl
  endif 
  ;; DRM Change: Should be double-sided CL to start, according to
  ;; the definition of ZSCORE used below.
  if not keyword_set(cl) then cl = 2.*gauss_pdf(1.)-1.  ; .682689  

  zscore = abs(gauss_cvf(1-0.5*(1-cl)))

  num = (size(ninput_2darr,/dim))[0] > 1 ; foil singularity
  sigczw_2darr = fltarr(num,2,/nozero)   ; low, high

  gd = where(ninput_2darr ne 0.,ngd,complement=bd)

  if gd[0] ne -1 then begin
     n_inv = 1. / ninput_2darr[gd] ; faster to calculate this once
     prob = nrec_2darr[gd] * n_inv

     numer1 = prob + 0.5 * n_inv * zscore^2
     numer2 = zscore * sqrt(prob * (1-prob) * n_inv + 0.25*zscore^2 * n_inv^2)
     denom_inv = 1. / (1 + n_inv * zscore^2)

     plo = (numer1 - numer2) * denom_inv
     phi = (numer1 + numer2) * denom_inv
     
     ;; Take absolute value because can get rounding errors that make
     ;; verys mall negative numbers
     sigczw_2darr[gd,0] = abs(prob - plo)
     sigczw_2darr[gd,1] = abs(phi - prob)
  endif 
  if bd[0] ne -1 then $
     sigczw_2darr[bd,*] = 0.    ; must instantiate

  return, sigczw_2darr

end                             ; sdss_calcsigbinom()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_medianw, x, w, even=even, double=double, percent=percent,$
                       silent=silent
  ;; Calculate the weighted median as described here:
  ;; http://stackoverflow.com/questions/9794558/weighted-median-computation
  ;; which quotes the R documentation which describes the weighted
  ;; median as:
  ;; "For the n elements x = c(x[1], x[2], ..., x[n]) with positive
  ;; weights w = c(w[1], w[2], ..., w[n]) such that sum(w) = S, the
  ;; weighted median is defined as the element x[k] for which initial
  ;; the total weight of all elements x[i] < x[k] is less or equal to
  ;; S/2 and for which the total weight of all elements x[i] > x[k] is
  ;; less or equal to S/2." 
  ;; So trying to follow the Stack Overflow example
  if n_params() ne 2 then begin
     print,'Syntax - sdss_medianw(x, w, [/even, percent=])'
     return,-1
  endif
  if not keyword_set(percent) then percent = 0.5

  nx = (size(x,/dim))[0] > 1

  if nx eq 1 then begin
     if not keyword_set(silent) then $
        print,'sdss_medianw(): cannot median scalar; return value'
     return,x
  endif

  if nx eq 2 then begin
     if not keyword_set(silent) then $
        print,'sdss_medianw(): cannot median two elements; return weighted mean'
     return,total(w*x)/total(w,double=double)
  endif

  unq = uniq(w,sort(w))         ; handles case of only 1 weight
  if n_elements(unq) eq 1 then begin
     ;; Actually equal weighting is like no weighting
     ;; IDL's should be faster
     if percent eq 0.5 then $
        return, median(x, even=even, double=double) 

     ;; Else do it by hand
     srt = sort(x)
     cumdist = total(x[srt],/cum,double=double,/nan)
     cumdist = cumdist/cumdist[nx-1]

     med = interpol(x[srt],cumdist,percent)
     return, (med > x[srt[0]]) < x[srt[nx-1]] ; don't interpolate beyond
  endif


  ;; Now for the real weighted median
  twgt = total(w,double=double,/nan)       ; total weight
  srt = sort(x)
  x_srt = x[srt]
  w_srt = w[srt]
  
  ;; *Do* need sum_hi limit 
  sum_lo = total(w_srt,/cum,double=double,/nan)  ; sum(x, 0, k)
  sum_lo = [0,(shift(sum_lo,1))[1:nx-1]]
  sum_hi = twgt - total(w_srt,/cum,double=double,/nan) ; sum(x, ngd-1, k)
  gd = where(sum_lo le percent*twgt and sum_hi le (1.-percent)*twgt,ngd)

  if ngd eq 0 then stop, 'sdss_medianw() stop: no matching percent.'

  if ngd eq 1 or not keyword_set(even) then begin
     med = x_srt[gd[0]]     ; last number
;     print,'sdss_medianw(): IDL quickie would return ',x_srt[gd[0]]
  endif else begin
     ;; Take the weighted mean
     med = total(x_srt[gd]*w_srt[gd],double=double,/nan)/$
           total(w_srt[gd],double=double,/nan)
;     print,'sdss_medianw(): multiple possible answers, '+$
;           string(percent,format='(f4.2)')+'*twgt = ',percent*twgt
;     printcol,gd,x_srt[gd],w_srt[gd],sum_lo[gd],sum_hi[gd]
;     print,'   weighted mean would be', med           
  endelse


;  ii = 0L
;  ;; sum is the total weight of all x[i] > x[k]
;  sum = twgt - w_srt[ii]
;  print,ii,sum
;  while sum gt (1.-percent)*twgt do begin ; x[i] > x[k]
;     ii++
;     sum = sum - w_srt[ii]
;     print,ii,sum
;  endwhile                      ; 
;
;  print,'sdss_medianw(): loop answer',x_srt[ii]
;;  stop
  return, med

end                             ; sdss_medianw()


;+
; NAME:
;  medianw
; PURPOSE: (one line)
;  minimize the weighted average deviation.
; DESCRIPTION:
;  This routine estimates the average data value, xmed, such that
;  the average deviation, total(weight * abs(x-xmed) ), is minimized.
;  For equally-weighted points, this is the median.
;
;  The statistics are robust in that the result is insensitive
;  to outliers.  
; CATEGORY:
;  Statistics
; CALLING SEQUENCE:
;  mean = medianw(x, w)
; INPUTS:
;  x   - Input data to be analyzed.
; OPTIONAL INPUT PARAMETERS:
;  w - weights.  Assumed equally weighted if not passed.
; INPUT KEYWORD PARAMETERS:
; OUTPUT KEYWORD PARAMETERS:
; OUTPUTS:
; COMMON BLOCKS:
;  None.
; SIDE EFFECTS:
;  None.
; RESTRICTIONS:
;  None.
; PROCEDURE:
;
;
; Minimizing f(xmed) = total(weight * abs(x-xmed)) is equivalent to finding
; df(xmed)/d xmed = 0 = -total(weight * sgn(x-xmed)).  
; If x is sorted, and xmed=x[imed], this can be written as
; 0 = total(weight[0:imed-1] * sgn(x[0:imed-1]-xmed)) 
;   + total(weight[imed+1:n-1] * sgn(x[imed+1:n-1]-xmed))
; or
; 0 = - total(weight[0:imed-1]) + total(weight[imed+1:n-1])
;
; or
; total(weight[0:imed-1]) = total(weight[imed+1:n-1])
;
;
; MODIFICATION HISTORY:
;  Written by Leslie A. Young, Soutwest Research Institute, 2002 Jul 23.
;  Modified LAY Oct 2002.  Better medianwTEST; 
;                          Avoid eq for comparing floats in test for one median value
;                          (equivalent to an odd number of equally weighted values)
;                          Average two neighbors if xmed is not one of the listed x's
;                          (equivalent to an even number of equally weighted values)
;  Modified LAY Oct 2002.  Changed name to medianw; 
;-
function sdss_minwad, x, w
  ;; Incoporated into sdss_functions.pro for sdss_stackciv.pro by KLC
  ;; on 27 Jul 2014
  if n_params() lt 1 then begin
     print,'Syntax - sdss_minwad(x,[w])'
     return,-1
  endif
  n = n_elements(x)
  noweight = (n_params() lt 2)
  if noweight then  w = replicate(1.,n)

  if noweight then begin
    med = median(x, /EVEN)
  end else begin
    s = sort(x)
    xx = float(x[s])
    ww = w[s]/total(w)
    i = 0
    totw = 0.
    while totw lt (1-(totw+ww[i])) and i lt n do begin
      totw = totw + ww[i]
      i = i + 1
    end
    if (abs( totw - (1-(totw+ww[i])) ) lt ww[i]*1e-6) then begin
      med = xx[i]
    end else begin
      med = ( 0.5*xx[i] + 0.5*xx[i-1] )
    end
  end
  
  return, med
  
end


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcncmplt, cxw_civ, sigcxw_civ, signcmplt=signcmplt, _extra=extra
  ;; Calculate weighted quantity (either completeness to return
  ;; completeness-corrected number or dX(W) to return dN/dz)
  if n_params() lt 1 then begin
     print,'Syntax - sdss_calcncmplt( cxw_civ, [sigcxw_civ, signcmplt=, _extra=])'
     return,-1
  endif 

  ;; Up-weight the numbers to reflect completeness
  ncmplt = total(1./cxw_civ)
  ;; Errors just in Number = sum(1 / C(W)) and technically should
  ;; account for Poisson errors here
  nobs = (size(cxw_civ,/dim))[0] > 1 ; foil singularity
  ;; _extra includes sigma= or cl=, /verbose
  signcmplt = sdss_calcsigpoiss(float(nobs),_extra=extra)
  if keyword_set(sigcxw_civ) then begin
     signcmplt[0] = sqrt(signcmplt[0]^2 + total((sigcxw_civ[*,0]/cxw_civ^2)^2))
     signcmplt[1] = sqrt(signcmplt[1]^2 + total((sigcxw_civ[*,1]/cxw_civ^2)^2))
  endif

  return, ncmplt
end                             ; sdss_calcncmplt()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcncolm, ncolm_civ, cxw_civ, signcolm_civ, sigcxw_civ, $
                         signcolm=signcolm, _extra=extra
  ;; Calculate weighted quantity (either completeness to return
  ;; completeness-corrected number or dX(W) to return dN/dz)
  if n_params() lt 3 then begin
     print,'Syntax - sdss_calcncolm( ncolm_civ, cxw_civ, signcolm_civ, '
     print,'                        [sigcxw_civ, signcolm=, _extra=])'
     return,-1
  endif 
  nobs = (size(ncolm_civ,/dim))[0] > 1 ; foil singularity
  if size(signcolm_civ,/n_dim) eq 2 then $
     sign_civ = signcolm_civ $               ; already correct dimensionality
  else sign_civ = rebin(signcolm_civ,nobs,2) ; expand

  if median(ncolm_civ) lt 20. then begin
     n_civ = 10.^ncolm_civ     ; must be in linear space
     sign_civ = sign_civ / (alog(10)*n_civ)
  endif else n_civ = ncolm_civ

  ;; Up-weight the numbers to reflect completeness
  ncolm = total(n_civ/cxw_civ)
  ;; _extra includes sigma= or cl=, /verbose
  signcolm = sdss_calcsigpoiss(float(nobs),_extra=extra)
  if keyword_set(sigcxw_civ) then begin
     ;; Numerically nicer to scale both
     signcolm[0] = signcolm[0]^2 + $
                   total((n_civ/cxw_civ)^2 * ((sign_civ[*,0]/n_civ)^2 + $
                         (sigcxw_civ[*,0]/cxw_civ)^2),/double) 
     signcolm[1] = signcolm[1]^2 + $
                   total((n_civ/cxw_civ)^2 * ((sign_civ[*,1]/n_civ)^2 + $
                         (sigcxw_civ[*,1]/cxw_civ)^2),/double) 
  endif else begin
     ;; Otherwise just fold in the column error
     signcolm[0] = signcolm[0]^2 + total((sign_civ[*,0]/cxw_civ)^2,/double)
     signcolm[1] = signcolm[1]^2 + total((sign_civ[*,1]/cxw_civ)^2,/double)
  endelse 
  signcolm = sqrt(signcolm) 

  return, ncolm
end                             ; sdss_calcncolm()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_setncolmflg, strct_fil, ncolm=ncolm, signcolm=signcolm, final=final
  ;; Do a more reasonable meta-analysis of the column density limits
  if n_params() ne 1 then begin
     print,'Syntax - sdss_setncolmflg(strct_fil, ncolm=, signcolm=, /final)'
     return,-1
  endif 

  ;; Open structure
  if size(strct_fil,/type) eq 7 then $
     civstr = xmrdfits(strct_fil, 1, /silent) $
  else civstr = strct_fil
  nstrct = (size(civstr,/dim))[0]
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endif else begin
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endelse 

  ewrto = civstr.(ewtag)[0]/civstr.(ewtag)[1] 
  sigewrto = abs(ewrto)*sqrt((civstr.(sigewtag)[0]/civstr.(ewtag)[0])^2+$
                             (civstr.(sigewtag)[1]/civstr.(ewtag)[1])^2)

  ;; Very conservative but probably realistic to think that EW should
  ;; be 
  unsat = where(ewrto gt 2.-sigewrto,complement=sat)

  ncolmflg = replicate(sdss_getlimflg(),nstrct)
  if sat[0] ne -1 then ncolmflg[sat] = sdss_setlimflg(ncolmflg[sat],/lower)

  if keyword_set(ncolm) then begin
     if median(ncolm) lt 20. then begin
        n = 10.^ncolm
        sn = signcolm * alog(10.) * n
     endif else begin
        n = ncolm
        sn = signcolm
     endelse 

     upper = where(n/sn lt 3.)
     if upper[0] ne -1 then $
        ncolmflg[upper] = sdss_setlimflg(ncolmflg[upper],/upper) ; may have contradiction
  endif

  return,ncolmflg
end                             ; sdss_setncolmflg()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcncolmdblt,strct_fil,dblt_name,signcolm=signcolm,log=log,$
                            subscript=subscript,final=final,single=single,$
                            eq_weight=eq_weight, ncolmflg=ncolmflg,$
                            estflg=estflg,silent=silent,debug=debug
  ;; Return error-weighted average of two doublet values
  ;; Adoped from civ_calcewn_ndblt().
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcncolmdblt(strct_fil,dblt_name,[signcolm=,/log,'
     print,'                            subscript=,/final,/eq_weight, ncolmflg=,'
     print,'                            /estflg,/silent,/debug])'
     return,-1
  endif 

  ;; Open structure
  if size(strct_fil,/type) eq 7 then $
     strct = xmrdfits(strct_fil, 1, /silent) $
  else strct = strct_fil
  nstrct = (size(strct,/dim))[0]
  tags = tag_names(strct[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ncolmtag = (where(tags eq 'NCOLM_FINAL'))[0]
     signcolmtag = (where(tags eq 'SIGNCOLM_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ncolmtag = (where(tags eq 'NCOLM_ORIG'))[0]
     signcolmtag = (where(tags eq 'SIGNCOLM_ORIG'))[0]
  endelse 
  flgtag = (where(tags eq 'NCOLMFLG'))[0]        

  if keyword_set(debug) then silent = debug
  
  fgood = sdss_getlimflg()           ; 1 
  flow = sdss_getlimflg(/lower)      ; 2
  fup = sdss_getlimflg(/upper)       ; 4
  fboth = sdss_setlimflg(fup,/lower) ; 6

  ;; Doublet
  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  ;; Figure out the location of the doublet
  subscript = replicate(-1,nstrct,2)
  ii = 0
  gd = where(strct.wrest[ii] eq dblt.wvI)
  if gd[0] ne -1 then subscript[gd,0] = ii
  ii++
  if not keyword_set(single) then begin
     gd = where(strct.wrest[ii] eq dblt.wvII)
     if gd[0] ne -1 then subscript[gd,1] = ii
     ii++
  endif 
  bd = where(subscript[gd,0] eq -1,nbd)
  while nbd ne 0 do begin
     gd = where(strct.wrest[ii] eq dblt.wvI)
     if gd[0] ne -1 then subscript[gd,0] = ii
     ii++
     if not keyword_set(single) then begin
        gd = where(strct.wrest[ii] eq dblt.wvII)
        if gd[0] ne -1 then subscript[gd,1] = ii
        ii++
     endif 
     bd = where(subscript[gd,0] eq -1,nbd)
  endwhile                      ; nbd ne 0

  if keyword_set(single) then begin
     ncolm = 10.^strct.(ncolmtag)[subscript[*,0]]
     signcolm = alog(10.)*strct.(signcolmtag)[subscript[*,0]]*ncolm
     ncolmflg = strct.(flgtag)[subscript[*,0]]
     bd = where(ncolm/signcolm lt 3.)
     if bd[0] ne -1 then ncolmflg[bd] = sdss_setlimflg(ncolmflg[bd],/upper)
     goto,skip_dblt
  endif else begin
     ;; Error-weighted average of AODM column density
     nI = fltarr(nstrct,/nozero)
     nII = nI
     signI = nI
     signII = nI
     ncolm = fltarr(nstrct,/nozero) ; want to make sure these get instantiated
     signcolm = fltarr(nstrct,/nozero)
     ncolmflg = intarr(nstrct)
  endelse 

  ;; Set up arrays
  ;; WVI
  unq = uniq(subscript[*,0],sort(subscript[*,0]))
  nindxI = n_elements(unq)      ; may be 1
  if nindxI eq 1 then begin
     flgI = strct.(flgtag)[subscript[0,0]] 
     gd = where(strct.(ncolmtag)[subscript[0,0]] gt 0.,complement=bd)
     if gd[0] ne -1 then begin
        nI[gd] = 10.^strct[gd].(ncolmtag)[subscript[0,0]]
        signI[gd] = strct[gd].(signcolmtag)[subscript[0,0]]*alog(10.)*nI[gd]
     endif
     if bd[0] ne -1 then begin
        nI[bd] = strct[bd].(signcolmtag)[subscript[0,0]] ; already linear
        signI[bd] = nI[bd]
        flgI[bd] = flgI[bd] - (flgI[bd] and fgood)
     endif 
     ;; Error-weighted average
     if keyword_set(eq_weight) then signI[*] = 1.d 
  endif else $
     stop,'sdss_calcncolmdblt() stop: not ready to handle crazy subcripts (I)'

  ;; WVII
  unq = uniq(subscript[*,1],sort(subscript[*,1]))
  nindxII = n_elements(unq)     ; may be 1
  if nindxII eq 1 then begin
     flgII = strct.(flgtag)[subscript[0,1]] 
     gd = where(strct.(ncolmtag)[subscript[0,1]] gt 0.,complement=bd)
     if gd[0] ne -1 then begin
        nII[gd] = 10.^strct[gd].(ncolmtag)[subscript[0,1]]
        signII[gd] = strct[gd].(signcolmtag)[subscript[0,1]]*alog(10.)*nII[gd]
     endif
     if bd[0] ne -1 then begin
        nII[bd] = strct[bd].(signcolmtag)[subscript[0,1]] ; already linear
        signII[bd] = nII[bd]
        flgII[bd] = flgII[bd] - (flgII[bd] and fgood)
     endif 
     ;; Error-weighted average
     if keyword_set(eq_weight) then signII[*] = 1.d
  endif else $
     stop,'sdss_calcncolmdblt() stop: not ready to handle crazy subcripts (II)'


  ;; ;;;;;;;
  ;; Check flags (1 = analyze; 2 = lower limit; 4 = upper limit);
  ;; see sdss_getlimflg()
  gd = where(flgI eq fgood and flgII eq fgood,complement=bd)
  if gd[0] ne -1 then begin
     ;; Take error-weighted value
     wgt = (1./signI[gd]^2 + 1./signII[gd]^2) 
     ncolm[gd] = (nI[gd]/signI[gd]^2 + nII[gd]/signII[gd]^2)/wgt
     signcolm[gd] = sqrt(1./wgt)
     ncolmflg[gd] = fgood
  endif                         ; both good measurements

  ;; ;;;;;;;
  if bd[0] ne -1 then begin
     ;; Look for one good measurement
     gd = where(flgI[bd] eq fgood,complement=bd2)
     if gd[0] ne -1 then begin
        ;; Take wvI values
        gd = bd[gd]
        ncolm[gd] = nI[gd]
        signcolm[gd] = signI[gd]
        ncolmflg[gd] = fgood

        ;; Checks
        chk = where((flgII[gd] and flow) eq flow and nI[gd] lt nII[gd],nchk)
        if chk[0] ne -1 and not keyword_set(silent) then begin
           chk = gd[chk]
           print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvII,2),$
                 ' lower limit > ',strtrim(dblt.wvI,2),' measurement',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nII[chk]-nI[chk])/signI[chk]
        endif                   ; nI < nII
        chk = where((flgII[gd] and fup) eq fup and nI[gd] gt nII[gd],nchk)
        if chk[0] ne -1 and not keyword_set(silent) then begin
           chk = gd[chk]
           print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvII,2),$
                 ' upper limit < ',strtrim(dblt.wvI,2),' measurement',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nI[chk]-nII[chk])/signI[chk]
        endif                            ; nI > nII
     endif                               ; using nI
     if bd2[0] ne -1 then bd = bd[bd2] $ ; shrink
     else goto,skip_dblt                 ; done
     gd = where(flgII[bd] eq fgood,complement=bd2)
     if gd[0] ne -1 then begin
        ;; Take wvII values
        gd = bd[gd]
        ncolm[gd] = nII[gd]
        signcolm[gd] = signII[gd]
        ncolmflg[gd] = fgood

        ;; Checks
        chk = where((flgI[gd] and flow) eq flow and nII[gd] lt nI[gd],nchk)
        if chk[0] ne -1 and not keyword_set(silent) then begin
           chk = gd[chk]
           print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvI,2),$
                 ' lower limit > ',strtrim(dblt.wvII,2),' measurement',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nI[chk]-nII[chk])/signII[chk]
        endif                   ; nII < nI
        chk = where((flgI[gd] and fup[gd]) eq fup[gd] and nII[gd] gt nI[gd],nchk)
        if chk[0] ne -1 and not keyword_set(silent) then begin
           chk = gd[chk]
           print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvI,2),$
                 ' upper limit <  ',strtrim(dblt.wvII,2),' measurement',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nII[chk]-nI[chk])/signII[chk]
        endif                            ; nII > nI
     endif                               ; using nII
     if bd2[0] ne -1 then bd = bd[bd2] $ ; shrink
     else goto,skip_dblt                 ; done
  endif                                  ; looking for one good measurement
  if bd[0] eq -1 then goto,skip_dblt
  

  ;; ;;;;;;;
  ;; No precise measurement so look for agreement in limits
  ;; Check upper limit first b/c may have values that are both upper
  ;; and lower limits but trust the upper limit more b/c based on
  ;; value being < 3sigma
  gd = where((flgI[bd] and fup) eq fup and (flgII[bd] and fup) eq fup,$
             complement=bd2)
  if gd[0] ne -1 then begin
     gd = bd[gd]
     wgt = (1./signI[gd]^2 + 1./signII[gd]^2) 
     ncolm[gd] = (nI[gd]/signI[gd]^2 + nII[gd]/signII[gd]^2)/wgt
     signcolm[gd] = sqrt(1./wgt)
     ncolmflg[gd] = flgI[gd] or flgII[gd] ; intersection

;     ;; Take the lower upper limit
;     useI = where(nI[gd] lt nII[gd],complement=useII)
;     if useI[0] ne -1 then begin
;        useI = gd[useI]
;        ncolm[useI] = nI[useI]
;        signcolm[useI] = signI[useI]
;        ncolmflg[useI] = flgI[useI]
;     endif
;     if useII[0] ne -1 then begin
;        useII = gd[useII]
;        ncolm[useII] = nII[useII]
;        signcolm[useII] = signII[useII]
;        ncolmflg[useII] = flgII[useII]
;     endif 
  endif                               ; both upper limits
  if bd2[0] ne -1 then bd = bd[bd2] $ ; shrink
  else goto,skip_dblt

  gd = where((flgI[bd] and flow) eq flow and (flgII[bd] and flow) eq flow,$
             complement=bd2)
  if gd[0] ne -1 then begin
     gd = bd[gd]
     wgt = (1./signI[gd]^2 + 1./signII[gd]^2) 
     ncolm[gd] = (nI[gd]/signI[gd]^2 + nII[gd]/signII[gd]^2)/wgt
     signcolm[gd] = sqrt(1./wgt)
     ncolmflg[gd] = flgI[gd] or flgII[gd] ; intersection

;     ;; Take the higher lower limit
;     useI = where(nI[gd] gt nII[gd],complement=useII)
;     if useI[0] ne -1 then begin
;        useI = gd[useI]
;        ncolm[useI] = nI[useI]
;        signcolm[useI] = signI[useI]
;        ncolmflg[useI] = flgI[useI]
;     endif
;     if useII[0] ne -1 then begin
;        useII = gd[useII]
;        ncolm[useII] = nII[useII]
;        signcolm[useII] = signII[useII]
;        ncolmflg[useII] = flgII[useII]
;     endif 
  endif                               ; taking lower upper
  if bd2[0] ne -1 then bd = bd[bd2] $ ; shrink
  else goto,skip_dblt

  ;; ;;;;;;;
  ;; Check for ambiguity (simultaneously lower and upper limit)
  gd = where((flgI[bd] and flow) eq flow and (flgII[bd] and fup) eq fup,$
             complement=bd2)
  if gd[0] ne -1 then begin
     gd = bd[gd]
     ave = where(nI[gd] lt nII[gd],complement=chk)
     if ave[0] ne -1 then begin
        ave = gd[ave]
        ncolm[ave] = 0.5*(nI[ave] + nII[ave])
        signcolm[ave] = 0.5*(nII[ave] - nI[ave])
        ncolmflg[ave] = fgood or 8 ; just the flag for average
     endif 
     
     if chk[0] ne -1 then begin
        chk = gd[chk]
        sigcomb = sqrt(signI[chk]^2+signII[chk]^2)
        ave = where(nI[chk] lt nII[chk]+sigcomb,complement=chk2,ncomplement=nchk)
        ;; Use larger error but still take average
        if ave[0] ne -1 then begin
           ave = chk[ave]
           ncolm[ave] = 0.5*(nI[ave] + nII[ave])
           signcolm[ave] = sigcomb[ave]
           ncolmflg[ave] = fgood or 8 ; just the flag for average
        endif 

        if chk2[0] ne -1 then begin
           chk = chk[chk2]
           if not keyword_set(silent) then $
              print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvI,2),$
                    ' lower limit > ',strtrim(dblt.wvII,2),' upper limit + 1sig',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nI[chk]-nII[chk])/sigcomb[chk]
        endif 
        ;; Any in chk2 is just going to puff off now b/c won't
        ;; be in bd2
     endif                            ; check ambiguity
  endif                               ; take average, nI < nII
  if bd2[0] ne -1 then bd = bd[bd2] $ ; shrink
  else goto,skip_dblt


  gd = where((flgI[bd] and fup) eq fup and (flgII[bd] and flow) eq flow,$
             complement=bd2)
  if gd[0] ne -1 then begin
     gd = bd[gd]
     ave = where(nI[gd] gt nII[gd],complement=chk)
     if ave[0] ne -1 then begin
        ave = gd[ave]
        ncolm[ave] = 0.5*(nI[ave] + nII[ave])
        signcolm[ave] = 0.5*(nI[ave] - nII[ave])
        ncolmflg[ave] = fgood or 8 ; just the flag for average
     endif 
     
     if chk[0] ne -1 then begin
        chk = gd[chk]
        sigcomb = sqrt(signI[chk]^2+signII[chk]^2)
        ave = where(nI[chk]+sigcomb gt nII[chk],complement=chk2,ncomplement=nchk)
        ;; Use larger error but still take average
        if ave[0] ne -1 then begin
           ave = chk[ave]
           ncolm[ave] = 0.5*(nI[ave] + nII[ave])
           signcolm[ave] = sigcomb[ave]
           ncolmflg[ave] = fgood or 8 ; just the flag for average
        endif 

        if chk2[0] ne -1 then begin
           chk = chk[chk2]
           if not keyword_set(silent) then $
              print,'sdss_calcncolmdblt(): ',strtrim(dblt.wvI,2),$
                    ' upper limit + 1sig < ',strtrim(dblt.wvII,2),' lower limit',nchk
           if keyword_set(debug) then $
              printcol,strct[chk].qso_name,strct[chk].(ztag)[0],$
                       flgI[chk],flgII[chk],(nI[chk]-nII[chk])/sigcomb[chk]
        endif 
        ;; Any in chk2 is just going to puff off now b/c won't
        ;; be in bd2
     endif 
  endif                               ; take average, nI > nII

  skip_dblt:

  ;; ;;;;;;;
  bd = where(ncolmflg eq 0,nbd)
  if bd[0] ne -1 then begin
     ;; Clean up of messy cases; take the lower error
     gd = where(nI[bd]/signI[bd] gt nII[bd]/signII[bd],complement=bd2)
     if gd[0] ne -1 then begin
        ;; nI is the winner
        gd = bd[gd]
        ncolm[gd] = nI[gd]
        signcolm[gd] = signI[gd]
        ncolmflg[gd] = flgI[gd] or 16 ; flag just to point out problem
     endif                             ; take nI
     if bd2[0] ne -1 then begin
        ;; nII is the winner
        bd2 = bd[bd2]
        ncolm[bd2] = nII[bd2]
        signcolm[bd2] = signII[bd2]
        ncolmflg[bd2] = flgII[bd2] or 16 ; flag just to point out problem
     endif                                ; take nII
     if not keyword_set(silent) then $
        print,'sdss_calcncolmdlbt(): Leftover number not handled above; take highest S/N measurement/limit',nbd
     if keyword_set(debug) then $
        printcol,strct[bd].qso_name,strct[bd].(ztag)[0],$
                 flgI[bd],flgII[bd],nI[bd],nI[bd]/signI[bd],nII[bd],nII[bd]/signII[bd]
  endif                                   ; clean up

  if keyword_set(estflg) then begin
     ncolmflg0 = ncolmflg
     ncolmflg = sdss_setncolmflg(strct,ncolm=ncolm,signcolm=signcolm,final=final)
     ncolmflg = ncolmflg + (ncolmflg0 and 8)  ; flag for average
     ncolmflg = ncolmflg + (ncolmflg0 and 16) ; flag for general problem     
  endif

  ;; Log 
  if keyword_set(log) then begin
     gd = where(ncolm gt 0.,ngd)
     if ngd ne 0 then begin
        signcolm[gd] = signcolm[gd]/(alog(10.)*ncolm[gd])
        ncolm[gd] = alog10(ncolm[gd])
     endif 
  endif 

  return,ncolm
end                             ; sdss_calcncolmdblt()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getewclim, cmplt_fil, clim=clim, grid=grid, nlosmax=nlosmax, $
                         ncolm=ncolm, dblt_name=dblt_name, log=log, _extra=extra
  ;; Compute completeness limit
  if n_params() ne 1 then begin
     print,'Syntax - sdss_getewclim(cmplt_fil, [clim=, nlosmax=, /ncolm, '
     print,'                        dblt_name=, /log, /grid, _extra=])'
     return,-1
  endif 
  if not keyword_set(clim) then clim = 0.5 ; 50% complete

  if size(cmplt_fil,/type) eq 8 then cmpltstr = cmplt_fil $
  else cmpltstr = xmrdfits(cmplt_fil,1,/silent)

  ewclim = fltarr(cmpltstr.nzbin,/nozero)
  nlosmax = lonarr(cmpltstr.nzbin,/nozero)

  ew_global = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  index = lindgen(cmpltstr.newbin) * cmpltstr.nzbin

  for zz=0,cmpltstr.nzbin-1 do begin
     tmp = cmpltstr.czw_2darr[zz+index] - clim
     gd = where(tmp ge 0.)
     if keyword_set(grid) then begin
        ;; take the closest grid point with small buffer (for both)
        mn = min(tmp[0:gd[0]+1],imn,/abs) 
        ewclim[zz] = ew_global[imn]  
     endif else begin
        ;; May have dip at beginning 
        bd = where(gd ne shift(gd,1)+1,nbd) ; not consecutive with prev (always first one)
        if nbd gt 1 then begin
           iend = gd[bd[nbd-1]]
           if iend+1 eq cmpltstr.newbin then begin
              ;; Funny dip below clim at high ew_global
              print,'sdss_getewclim(): WARNING!!! dip below clim at high EW'
              iend = gd[0]      ; asusme OK
           endif
        endif else iend = gd[0]
        ;; Interpolate
        ;; _extra includes /lsquadratic, /quadratic, /spline
        ewclim[zz] = interpol(ew_global[0:iend+1], tmp[0:iend+1], $
                              0., _extra=extra)
     endelse 
;     if ewclim[zz] lt 0. or ewclim[zz] gt 1. then stop ; Angstrom
     if ewclim[zz] lt 0. then ewclim[zz] = ew_global[0] ; force sensible

     nlosmax[zz] = max(cmpltstr.rz_2darr[zz+index])

  endfor ; loop zz=cmpltstr.nzbin


  if keyword_set(ncolm) then begin
     if size(dblt_name,/type) eq 7 then $
        dblt = dblt_retrieve(dblt_name) $
     else dblt = dblt_name
     ncolmclim = ew_to_colm(replicate(dblt.wvI,cmpltstr.nzbin),$
                            ewclim*1.e3,/silent)
     ewclim = ncolmclim
  endif 

  if keyword_set(log) then ewclim = alog10(ewclim)

  return, ewclim
end                             ; sdss_getewclim()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_catcmpltstr, cmpltlo_fil, cmplthi_fil, ewseam, _extra=extra
  ;; Concatenate two completeness files
  if n_params() ne 3 then begin
     print,'Syntax - sdss_catcmpltstr(cmpltlo_fil, cmplthi_fil, ewseam, [_extra=])'
     return,-1
  endif
  
  ;; Read files
  if size(cmpltlo_fil,/type) eq 7 then $
     cmpltlo = xmrdfits(cmpltlo_fil,1,/silent) $
  else cmpltlo = cmpltlo_fil
  if size(cmplthi_fil,/type) eq 7 then $
     cmplthi = xmrdfits(cmplthi_fil,1,/silent) $
  else cmplthi = cmplthi_fil

  ;; Checks
  if cmpltlo.czn ne cmplthi.czn then $
     stop,'sdss_catcmpltstr() stop: czn flags do not match in the completeness structures'
  ewlim = [cmpltlo.ewlim[0], cmplthi.ewlim[1]] ; by design
  if cmpltlo.ewbinsize ne cmplthi.ewbinsize then begin
     print,'sdss_catcmpltstr(): NOTE!!! ewbinsize not equal in completeness files',$
           cmpltlo.ewbinsize, cmplthi.ewbinsize
     ewbinsize = cmpltlo.ewbinsize < cmplthi.ewbinsize 
  endif else ewbinsize = cmpltlo.ewbinsize
  ew_global = sdss_mkewarr(ewlim, ewbinsize, newbin=newbin)

  if cmpltlo.ewlim[1]+cmpltlo.ewbinsize lt ewseam then $
     stop,'sdss_catcmpltstr() stop: lower completeness structure below ewseam'
  if cmplthi.ewlim[0] gt ewseam then $
     stop,'sdss_catcmpltstr() stop: higher completeness structure excludes ewseam'
  bd = where(cmpltlo.zlim ne cmplthi.zlim)
  if bd[0] ne -1 then $
     print,'sdss_catcmpltstr(): redshift binning does not match in two completeness structures'
  bd = where(cmpltlo.cosmology ne cmplthi.cosmology)
  if bd[0] ne -1 then $
     print,'sdss_catcmpltstr() stop: cosmology does not match in two completeness structures'

  ;; Make arrays
  ew_global_lo = sdss_mkewarr(cmpltlo.ewlim, cmpltlo.ewbinsize, newbin=newbin_lo)
  ew_global_hi = sdss_mkewarr(cmplthi.ewlim, cmplthi.ewbinsize, newbin=newbin_hi)
  
  ;; Find best seam possible (bins defined on right-hand side)
  mn_lo = min(ew_global_lo+cmpltlo.ewbinsize - ewseam, imn_lo, /abs)
  mn_hi = min(ew_global_hi - ewseam, imn_hi, /abs)
  mn = min(ew_global - ewseam, imn, /abs)

  ;; ;;;;;;;
  ;; Create output (matches sdss_completeness_czw)
  f2darr = fltarr(newbin,/nozero)
  cmpltstr = { $
             ;; General parameters, so can re-create this
             list_fil:[cmpltlo.list_fil,cmplthi.list_fil], $
             dblt_name:cmpltlo.dblt_name, $
             rec_param:cmpltlo.rec_param or cmplthi.rec_param, $
             civobs_corr:cmpltlo.civobs_corr or cmplthi.civobs_corr, $
             userbias:cmpltlo.userbias or cmplthi.userbias, $
             czn:cmpltlo.czn, $
             ewseam:ewseam, $
             nzbin:cmpltlo.nzbin, $
             zlim:cmpltlo.zlim, $
             zbinsize:cmpltlo.zbinsize, $
             newbin:newbin, $
             ewlim:ewlim, $
             ewbinsize:ewbinsize, $
             cosmology:cmpltlo.cosmology, $
             ;; Now for actual results in this bin
             ninput_2darr:f2darr, $
             nrec_2darr:f2darr, $
             rz_2darr:f2darr, $ ; can reform(2darr, nzbin, newbin)
             dz_2darr:f2darr, $
             dx_2darr:f2darr, $
             czw_2darr:f2darr, $
             sigczw_2darr:fltarr(newbin,2,/nozero) $ ; [*,0] = lower, [*,1] = upper
             }

  ;; ;;;;;;;
  ;; lower completeness file
  if imn_lo-1 eq imn then begin
     ;; Since begins the same, can just paste in
     cmpltstr.ninput_2darr[0:imn-1] = cmpltlo.ninput_2darr[0:imn_lo-1]
     cmpltstr.nrec_2darr[0:imn-1] = cmpltlo.nrec_2darr[0:imn_lo-1]
     cmpltstr.rz_2darr[0:imn-1] = cmpltlo.rz_2darr[0:imn_lo-1]
     cmpltstr.dz_2darr[0:imn-1] = cmpltlo.dz_2darr[0:imn_lo-1]
     cmpltstr.dx_2darr[0:imn-1] = cmpltlo.dx_2darr[0:imn_lo-1]
     cmpltstr.czw_2darr[0:imn-1] = cmpltlo.czw_2darr[0:imn_lo-1]
     cmpltstr.sigczw_2darr[0:imn-1,0] = cmpltlo.sigczw_2darr[0:imn_lo-1,0]
     cmpltstr.sigczw_2darr[0:imn-1,1] = cmpltlo.sigczw_2darr[0:imn_lo-1,1]
  endif else begin
     ;; Interpolate
     ;; _extra includes /lsquadratic, /quadratic, /spline
     print,'sdss_catcmpltstr(): interpolating lower completeness structure'
     cmpltstr.ninput_2darr[0:imn-1] = interpol(cmpltlo.ninput_2darr[0:imn_lo-1], $
                                               ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.nrec_2darr[0:imn-1] = interpol(cmpltlo.nrec_2darr[0:imn_lo-1], $
                                             ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.rz_2darr[0:imn-1] = interpol(cmpltlo.rz_2darr[0:imn_lo-1], $
                                           ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.dz_2darr[0:imn-1] = interpol(cmpltlo.dz_2darr[0:imn_lo-1], $
                                           ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.dx_2darr[0:imn-1] = interpol(cmpltlo.dx_2darr[0:imn_lo-1], $
                                           ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.czw_2darr[0:imn-1] = interpol(cmpltlo.czw_2darr[0:imn_lo-1], $
                                            ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.sigczw_2darr[0:imn-1,0] = interpol(cmpltlo.sigczw_2darr[0:imn_lo-1,0], $
                                                 ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
     cmpltstr.sigczw_2darr[0:imn-1,1] = interpol(cmpltlo.sigczw_2darr[0:imn_lo-1,1], $
                                                 ew_global_lo[0:imn_lo-1], ew_global[0:imn-1], _extra=extra)
  endelse 

  ;; ;;;;;;;
  ;; higher completeness file
  if cmplthi.newbin-imn_hi eq newbin-imn then begin
     ;; Since begins the same, can just paste in
     cmpltstr.ninput_2darr[imn:*] = cmplthi.ninput_2darr[imn_hi:*]
     cmpltstr.nrec_2darr[imn:*] = cmplthi.nrec_2darr[imn_hi:*]
     cmpltstr.rz_2darr[imn:*] = cmplthi.rz_2darr[imn_hi:*]
     cmpltstr.dz_2darr[imn:*] = cmplthi.dz_2darr[imn_hi:*]
     cmpltstr.dx_2darr[imn:*] = cmplthi.dx_2darr[imn_hi:*]
     cmpltstr.czw_2darr[imn:*] = cmplthi.czw_2darr[imn_hi:*]
     cmpltstr.sigczw_2darr[imn:*,0] = cmplthi.sigczw_2darr[imn_hi:*,0]
     cmpltstr.sigczw_2darr[imn:*,1] = cmplthi.sigczw_2darr[imn_hi:*,1]
  endif else begin
     ;; Interpolate
     ;; _extra includes /lsquadratic, /quadratic, /spline
     print,'sdss_catcmpltstr(): interpolating higher completeness structure'
     cmpltstr.ninput_2darr[imn:*] = interpol(cmplthi.ninput_2darr[imn_hi:*], $
                                               ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.nrec_2darr[imn:*] = interpol(cmplthi.nrec_2darr[imn_hi:*], $
                                             ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.rz_2darr[imn:*] = interpol(cmplthi.rz_2darr[imn_hi:*], $
                                           ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.dz_2darr[imn:*] = interpol(cmplthi.dz_2darr[imn_hi:*], $
                                           ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.dx_2darr[imn:*] = interpol(cmplthi.dx_2darr[imn_hi:*], $
                                           ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.czw_2darr[imn:*] = interpol(cmplthi.czw_2darr[imn_hi:*], $
                                            ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.sigczw_2darr[imn:*,0] = interpol(cmplthi.sigczw_2darr[imn_hi:*,0], $
                                                 ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
     cmpltstr.sigczw_2darr[imn:*,1] = interpol(cmplthi.sigczw_2darr[imn_hi:*,1], $
                                                 ew_global_hi[imn_hi:*], ew_global[imn:*], _extra=extra)
  endelse 

;  stop
  return, cmpltstr
end                             ; sdss_catcmpltstr()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getdxw, cmplt_fil, zabs, ewabs, sigewabs=sigewabs, dz=dz, $
                      sigma=sigma, grid=grid, ewmax=ewmax, czw=czw, $
                      dxwmax=dxwmax, silent=silent, civobs_corr=civobs_corr, $
                      final=final, dxciv=dxciv, _extra=extra
  ;; For each zabs and ewabs, find the dX(W) in the correct redshift bin
  if n_params() ne 3 then begin
     print,'Syntax - sdss_getdxw( cmplt_fil, zabs, ewabs, [sigewabs=, /dz,'
     print,'                      dxwmax=, sigma=, /grid, /ewmax, dxciv=, '
     print,'                      /czw, /silent, civobs_corr=, /final, _extra=])'
     return,-1
  endif 

  if size(cmplt_fil,/type) eq 8 then cmpltstr = cmplt_fil $
  else cmpltstr = xmrdfits(cmplt_fil,1,/silent)
  tags = tag_names(cmpltstr)
  if keyword_set(dz) then $
     dxtag = (where(tags eq 'DZ_2DARR'))[0] $
  else dxtag = (where(tags eq 'DX_2DARR'))[0] 

  if keyword_set(czw) then begin
     ;; Return completess (0 to 1) 
     dxw_2darr = cmpltstr.czw_2darr 
     sigdxw_2darr = cmpltstr.sigczw_2darr
  endif else begin
     ;; Return completeness-corrected dX(W) or dz(W)
     dxw_2darr = cmpltstr.(dxtag)*cmpltstr.czw_2darr
     tmp = rebin(cmpltstr.(dxtag),cmpltstr.nzbin*cmpltstr.newbin,2) ; duplicate dimen
     sigdxw_2darr = tmp*cmpltstr.sigczw_2darr                       ; [*,0] = lower, [*,1] = upper
  endelse 
  dxwmax = cmpltstr.(dxtag)[0:cmpltstr.nzbin-1] ; 100% complete

  if keyword_set(civobs_corr) then begin
     ;; want to account for detected abosrbers blocking out pathlength
     if keyword_set(cmpltstr.czn) then $
        stop,'sdss_getdxw() stop: WARNING! civobs_corr being applied to column density completeness'

     if size(civobs_corr,/type) eq 8 then civstr = civobs_corr $
     else civstr = xmrdfits(civobs_corr,1,/silent)
     nciv = (size(civstr,/dim))[0] > 1 ; foil singularity

     tags = tag_names(civstr[0])
     if keyword_set(final) then begin
        wvlimtag = (where(tags eq 'WVLIM_FINAL'))[0] 
        zabstag = (where(tags eq 'ZABS_FINAL'))[0] 
        ewtag = (where(tags eq 'EW_FINAL'))[0] 
     endif else begin
        wvlimtag = (where(tags eq 'WVLIM_ORIG'))[0] 
        zabstag = (where(tags eq 'ZABS_ORIG'))[0] 
        ewtag = (where(tags eq 'EW_ORIG'))[0] 
     endelse 
     civstr = civstr[sort(civstr.(ewtag)[0])] ; make later things easier
     dblt = dblt_retrieve(civstr[0].wrest[0])
     zcivlo = civstr.(wvlimtag)[0,0] / dblt.wvI - 1.
     zcivhi = civstr.(wvlimtag)[1,1] / dblt.wvI - 1. ; 1550 upper limit

     cosmology = sdss_setcosmology(cosmology=cmpltstr.cosmology)
     dxciv = cosm_xz(zcivhi,zmin=zcivlo,/exact,/noinit,/silent)

     ;; Absorbers may overlap and don't want to double count
     ;; but not a quick search to work across sightlines
  endif

  if (size(cmpltstr.zlim,/n_dim))[0] gt 1 then begin
     z_global = cmpltstr.zlim[*,0] 
     ;; ignore cmpltstr.zbinsize b/c has z fudge factor
     dz_global = cmpltstr.zlim[*,1] - cmpltstr.zlim[*,0] 
  endif else begin
     z_global = sdss_mkzarr(cmpltstr.zlim, cmpltstr.zbinsize)
     dz_global = replicate(cmpltstr.zbinsize, cmpltstr.nzbin)
  endelse 
  ew_global = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  index = lindgen(cmpltstr.newbin) * cmpltstr.nzbin

  ;; Set up output arrays
  nabs = (size(zabs,/dim))[0] > 1
  ntest = (size(ewabs,/dim))[0] > 1
  if nabs ne ntest then $
     stop,'sdss_getdxw() stop: WARNING! redshift and EW array dimensions not equal'
  dxw_abs = fltarr(nabs,/nozero)
  sigma = fltarr(nabs,2,/nozero)

  ;; Test overall fit of objects
  bd = where(zabs lt z_global[0] or $
             zabs gt z_global[cmpltstr.nzbin-1]+dz_global[cmpltstr.nzbin-1],nbd)
  if nbd ne 0 then begin
     if not keyword_set(silent) then $
        print,'sdss_getdxw(): number out of z bounds',$
              z_global[0], z_global[cmpltstr.nzbin-1]+dz_global[cmpltstr.nzbin-1],$
           nbd
     dxw_abs[bd] = !values.f_infinity ; make everything zero
     sigma[bd,*] = 0.
  endif 

  ;; Loop over given redshift bins and interpolate values
  ;; Assumes nzbin the faster way to go
  for zz=0,cmpltstr.nzbin-1 do begin
     sub = where(zabs ge z_global[zz] and $
                 zabs lt z_global[zz]+dz_global[zz],nsub)

     if sub[0] eq -1 then continue ; nothing observed in this z bin

     if keyword_set(grid) then begin
        ;; Take whole box values
        iew_abs = floor((ewabs[sub] - cmpltstr.ewlim[0])/cmpltstr.ewbinsize)
        indx_abs = zz + iew_abs * cmpltstr.nzbin
        dxw_abs[sub] = dxw_2darr[indx_abs]
        sigma[sub,0] = sigdxw_2darr[indx_abs,0]
        sigma[sub,1] = sigdxw_2darr[indx_abs,1]

     endif else begin
        ;; Interpolate values
        ;; _extra includes /lsquadratic, /quadratic, /spline
        dxw_abs[sub] = $
           interpol(dxw_2darr[zz+index], ew_global, ewabs[sub], _extra=extra)
        sigma[sub,0] = $
           interpol(sigdxw_2darr[zz+index,0], ew_global, ewabs[sub], _extra=extra)
        sigma[sub,1] = $
           interpol(sigdxw_2darr[zz+index,1], ew_global, ewabs[sub], _extra=extra)
     endelse 

     ;; Fix overflow but don't return flat zeros
     gd = where(dxw_2darr[zz+index] gt 0.)
     lodefault = min(dxw_2darr[zz+index[gd]],ilodefault,$
                     max=hidefault,subscript_max=ihidefault)
     ilodefault = gd[ilodefault]
     ihidefault = gd[ihidefault]

     bdlo = where(ewabs[sub] lt cmpltstr.ewlim[0],nbdlo)
     if nbdlo ne 0 then begin
        if not keyword_set(silent) then $
           print,'sdss_getdxw(): number below lower EW limit',$
                 cmpltstr.ewlim[0],nbdlo
        dxw_abs[sub[bdlo]] = lodefault
        sigma[sub[bdlo],0] = sigdxw_2darr[zz+index[ilodefault],0]
        sigma[sub[bdlo],1] = sigdxw_2darr[zz+index[ilodefault],1]
     endif
     
     ;; Since can have stochastic fall-off (currently: 13 Dec 2011)
     ;; for high EW if not tested, the force the upper dX(W) to be on
     ;; par with the max ones; mostly works
     if keyword_set(ewmax) then $
        bdhi = where(ewabs[sub] gt ew_global[ihidefault],nbdhi) $
     else bdhi = where(ewabs[sub] gt cmpltstr.ewlim[1],nbdhi)
     if nbdhi ne 0 then begin
        if not keyword_Set(silent) then begin
           if keyword_set(ewmax) then $
              print,'sdss_getdxw(): number above forced EW max',$
                    ew_global[ihidefault],nbdhi $
           else $
              print,'sdss_getdxw(): number above EW limit',$
                    cmpltstr.ewlim[1],nbdhi
        endif                   ; /silent
        dxw_abs[sub[bdhi]] = hidefault
        sigma[sub[bdhi],0] = sigdxw_2darr[zz+index[ihidefault],0]
        sigma[sub[bdhi],1] = sigdxw_2darr[zz+index[ihidefault],1]
     endif

     if keyword_set(civobs_corr) then begin
        ;; Correct for absorbers in bin, only block out everything
        ;; that has a higher EW (assuming it's wider)
        iciv = where(civstr.(zabstag)[0] ge z_global[zz] and $
                     civstr.(zabstag)[0] lt z_global[zz]+dz_global[zz] and $
                     finite(civstr.(ewtag)[0]) eq 1,niciv) ; update 16 Aug 2013 b/c new sdss_ewciv()
        if iciv[0] ne -1 then begin
           if niciv gt 1 then begin
              ;; civstrct must be sorted by EW; want to subtract the dX
              ;; blocked by every *other* absorber with greater EW
              dx_blocked = total(dxciv[iciv]) - $
                           [0.,(total(dxciv[iciv],/cum))[0:niciv-2]]
              ;; interpolate for the requested values
              dxw_rm = interpol(dx_blocked, civstr[iciv].(ewtag)[0], ewabs)
           endif else begin
              ;; Step function; let the checks below handle
              dx_blocked = total(dxciv[iciv])
              dxw_rm = ewabs * 0. ; zero out array
           endelse 

           bd = where(ewabs lt civstr[iciv[0]].(ewtag)[0])
           if bd[0] ne -1 then dxw_rm[bd] = dx_blocked[0] ; max off
           bd = where(ewabs gt civstr[iciv[niciv-1]].(ewtag)[0])
           if bd[0] ne -1 then dxw_rm[bd] = dx_blocked[niciv-1] ; full path

           ;; Sanity check
           test = where(dxw_rm lt 0.)
           if test[0] ne -1 then $
              stop,'sdss_getdxw() stop: WARNING!!! Negative dxw_rm'

           ;; Remember: dxw_abs is a completeness-corrected value so
           ;; technically, the subtracted part can be bigger than the
           ;; *completeness-corrected* path available to the absorber.
           ;; Take the fraction obscured, which makes adding back in
           ;; an absorber's dX illogical.
           dxw_abs[sub] = dxw_abs[sub] * (1 - dxw_rm / dxwmax[zz])
           sigma[sub,0] = sigma[sub,0] * (1 - dxw_rm / dxwmax[zz])
           sigma[sub,1] = sigma[sub,1] * (1 - dxw_rm / dxwmax[zz])
        endif                                   ; iciv[0] ne -1
     endif                                      ; civobs_corr= 

  endfor                        ; loop zz=cmpltstr.nzbin


  if keyword_set(sigewabs) then begin
     ;; Recursively call this function for the EW +/- sigma on the
     ;; given values but technically there's some error in
     ;; which absorbers have the small pathelength modification due to
     ;; dX_CIV_block if civobs_corr set...
     sigdxwlo = dxw_abs - sdss_getdxw( cmpltstr, zabs, ewabs-sigewabs, dz=dz, $
                                       grid=grid, ewmax=ewmax, czw=czw, $
                                       civobs_corr=civobs_corr, dxciv=0, $
                                       final=final, _extra=extra)
     sigdxwhi = sdss_getdxw( cmpltstr, zabs, ewabs+sigewabs, dz=dz, $
                             grid=grid, ewmax=ewmax, czw=czw, $
                             civobs_corr=civobs_corr, dxciv=0, final=final, $
                             _extra=extra) $
                - dxw_abs
     sigma[*,0] = sqrt(sigma[*,0]^2 + sigdxwlo^2)
     sigma[*,1] = sqrt(sigma[*,1]^2 + sigdxwhi^2)
  endif 

  return,dxw_abs
end                             ; sdss_getdxw()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntdndx, dndx_fil, skip_null=skip_null, out_fil=out_fil
  if n_params() ne 1 then begin
     print,'Syntax - sdss_prntdndx, dndx_fil, [/skip_null,out_fil=]'
     return
  endif 

  if size(dndx_fil,/type) eq 7 then $
     dndx = xmrdfits(dndx_fil,1,/silent) $
  else dndx = dndx_fil
  nstrct = (size(dndx,/dim))[0] > 1

  if keyword_set(out_fil) then $
     openw,1,out_fil

  for ss=0,nstrct-1 do begin
     if keyword_set(dndx[ss].czn) then $
        lin0 = strtrim(string(dndx[ss].ewlim[0],format='(f7.2)'),2)+' <= log N < '+$
               strtrim(string(dndx[ss].ewlim[1],format='(f7.2)'),2) $
     else $
        lin0 = strtrim(string(dndx[ss].ewlim[0],format='(f7.2)'),2)+' Ang <= EW < '+$
               strtrim(string(dndx[ss].ewlim[1],format='(f7.2)'),2)+' Ang'
     if dndx[ss].dz eq 1 then colhead = ['dz','dN/dz'] $
     else colhead = ['dX','dN/dX']
     lin = string('zlo','zhi','Num0','Errlo','Errhi','Num','Errlo','Errhi',$
                  colhead[0],'Errlo','Errhi',colhead[1],'Errlo','Errhi',$
                  format='(2(a6,1x),3(a5,1x),3(a8,1x),3(a10,1x),3(a8,1x))')
     if keyword_set(out_fil) then begin
        printf,1,''
        printf,1,lin0
        printf,1,lin
     endif else begin
        print,''
        print,lin0
        print,lin
     endelse

     nzbin = (size(dndx[ss].zcenter,/dim))[0] > 1
     for zz=0,nzbin-1 do begin
        if keyword_set(skip_null) and dndx[ss].dndx[zz] eq 0. then $
           continue 
        lin = $
           string(dndx[ss].zcenter[zz]-dndx[ss].sigzcenter[zz,0],$
                  dndx[ss].zcenter[zz]+dndx[ss].sigzcenter[zz,1],$
                  dndx[ss].numtot0[zz],dndx[ss].signumtot0[zz,0],$
                  dndx[ss].signumtot0[zz,1],dndx[ss].numtot[zz],$
                  dndx[ss].signumtot[zz,0],dndx[ss].signumtot[zz,1],$
                  dndx[ss].dxtot[zz],$
                  dndx[ss].sigdxtot[zz,0],dndx[ss].sigdxtot[zz,1],$
                  dndx[ss].dndx[zz],dndx[ss].sigdndx[zz,0],$
                  dndx[ss].sigdndx[zz,1],$
                  format='(2(f6.2,1x),i5,1x,2(f5.1,1x),f8.1,1x,2(f8.1,1x),3(f10.2,1x),3(f8.4,1x))')
        if keyword_set(out_fil) then printf,1,lin $
        else print,lin
     endfor                     ; loop zz=nzbin
  endfor                        ; loop ss=nstrct

  if keyword_set(out_fil) then begin
     close,1
     print,'sdss_prntdndx: created ',out_fil
  endif

end                             ; sdss_prntdndx


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltdndx, dndxstrct_fil, dndx2_fil=dndx2_fil, xrng=xrng,$
                  yrng=yrng, title=title, psize=psize, $
                  csize=csize, lthick=lthick, linsty=linsty, psfil=psfil,$
                  _extra=extra
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltdndx, dndxstrct_fil, [dndx2_fil=, xrng=,'
     print,'                    yrng=, title=, psize=, '
     print,'                    csize=, lthick=, linsty=, psfil=, _extra=]'
     return
  endif
  angstrom = STRING("305B)   
  gesign = ' !9'+string("263B)+'!X '
  lesign = ' !9'+string("243B)+'!X '
 
  if not keyword_set(csize) then csize = 2.
  if not keyword_set(psize) then psize = 1.5
  if not keyword_set(lthick) then lthick = 1.5
  if not keyword_set(scale) then scale = 1.
 
  if keyword_set(psfil) then begin
     ;; _extra= includes /encaps
     x_psopen,psfil,/maxs,/portrait,_extra=extra
     !p.multi = [1,1,1]
     !x.margin = [7.8,1.5]      ; left and right border
     !y.margin = [3.2,0.5]      ; bottom and top border
  endif 

  if size(dndxstrct_fil,/type) eq 7 then $
     dndx = xmrdfits(dndxstrct_fil,1,/silent) $
  else dndx = dndxstrct_fil
  clr = getcolor(/load)
  if not keyword_set(title) then $
     title = ['!8z!X!Dabs!N','d!8N!X/d!8X!X (!8W!X!Dr!N'+gesign+$
              string(dndx.ewlim[0],format='(f3.1)')+' '+angstrom+')']

  ;; Plot
  if not keyword_set(xrng) then $
     xrng = [min(dndx.zcenter-dndx.sigzcenter[*,0]),max(dndx.zcenter+dndx.sigzcenter[*,1])]
  if not keyword_set(yrng) then $
     yrng = [min(dndx.dndx-dndx.sigdndx[*,0]),max(dndx.dndx+dndx.sigdndx[*,1])]

  plot,[0],[0],/nodata,/xsty,/ysty,color=clr.black,$
       background=clr.white,thick=lthick,charsize=csize,$
       xtitle=title[0],ytitle=title[1],xrange=xrng,yrange=yrng

  oploterror,dndx.zcenter,dndx.dndx,dndx.sigzcenter[*,0],dndx.sigdndx[*,0],$
             /lobar, errcolor=clr.black,errthick=lthick,$
             psym=3,color=clr.black,thick=lthick,/nohat
  oploterror,dndx.zcenter,dndx.dndx,dndx.sigzcenter[*,1],dndx.sigdndx[*,1],$
             /hibar, errcolor=clr.black,errthick=lthick,$
             psym=3,color=clr.black,thick=lthick,/nohat

  if keyword_set(dndx2_fil) then begin
     ;; Plot the extra data
     if size(dndx2_fil,/type) eq 7 then $
        dndx2 = xmrdfits(dndx2_fil,1,/silent) $
     else dndx2 = dndx2_fil
     
     oploterror,dndx2.zcenter,dndx2.dndx,dndx2.sigzcenter[*,0],dndx2.sigdndx[*,0],$
                /lobar, errcolor=clr.red,errthick=lthick,$
                psym=3,color=clr.red,thick=lthick,/nohat
     oploterror,dndx2.zcenter,dndx2.dndx,dndx2.sigzcenter[*,1],dndx2.sigdndx[*,1],$
                /hibar, errcolor=clr.red,errthick=lthick,$
                psym=3,color=clr.red,thick=lthick,/nohat
  endif 

  if keyword_set(psfil) then begin
     x_psclose
     print,'sdss_pltdndx: created ',psfil
  endif

end                             ; sdss_pltdndx



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getuserbiasindex, biasuser_fil, zlim, tp_coeff=tp_coeff, $
                                fp_coeff=fp_coeff
  ;; For given redshift bin, figure out best-matching userbias bin
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getuserbiasindex(biasuser_fil, zlim)'
     return,-1
  endif 
  if size(zlim,/n_dim) ne 1 then begin
     print,'sdss_getuserbiasindex(): zlim must be [min, max] array'
     return,-1
  endif
  dz = (zlim[1]-zlim[0])
  zmid = 0.5*total(zlim)
  ztol = 1.e-3                  ; for matching

  if size(biasuser_fil,/type) eq 7 then $
     userbiasstr = xmrdfits(biasuser_fil,1,/silent) $
  else userbiasstr = biasuser_fil

  nzdim_user = size(userbiasstr.zlim,/n_dim)
  if nzdim_user eq 1 then nzbin_user = 1 $
  else nzbin_user = (size(userbiasstr.zlim[*,0],/dim))[0] 

  if nzbin_user eq 1 then begin
     ;; Only one option
     print,'sdss_getuserbiasindex(): zbin, user_indx',zlim,0,$
           format='(a,1x,2(f6.4,1x),i2)'
     tp_coeff = userbiasstr.coeff_tp ; [2,3]
     fp_coeff = userbiasstr.coeff_afp ; [2,3]
  endif else begin
     ;; Multiple options, find best overlap 
     dzrto = dz/(userbiasstr.zlim[*,1]-userbiasstr.zlim[*,0]) 
     invdzrto = (userbiasstr.zlim[*,1]-userbiasstr.zlim[*,0])/dz
     
     ;; If completeness zbin of interest is outside the bounds
     ;; of the userbias zbins, zero it out
     bd = where(zmid lt userbiasstr.zlim[*,0]+ztol or $
                zmid ge userbiasstr.zlim[*,1]-ztol,complement=gd)
     if bd[0] ne -1 then dzrto[bd] = 0. ; no overlap!
     if gd[0] eq -1 then stop,'sds_getuserbiasindex() stop: no non-zero dz ratio!'
     
     ;; if the completeness zbin is bigger than userbias zbin, then
     ;; look at the inverse
     bd = where(dzrto gt 1.)
     if bd[0] ne -1 then dzrto[bd] = invdzrto[bd]

     ;; Now find the biggest overlap
     mx = max(dzrto,user_indx)
     print,'sdss_getuserbiasindex(): zbin, user_indx',zlim,user_indx,$
           format='(a,1x,2(f6.4,1x),i2)'
     tp_coeff = reform(userbiasstr.coeff_tp[user_indx,*,*]) ; [2,3]
     fp_coeff = reform(userbiasstr.coeff_afp[user_indx,*,*]) ; [2,3]
  endelse

  return, user_indx
end                             ; sdss_getuserbiasindex()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcdndx, civstrct_fil, cmplt_fil, ewlim, final=final, dz=dz, $
                        silent=silent, afp_flg=afp_flg, czn=czn, _extra=extra
  ;; Calculate loads of stuff associated with dN/dX (or dN/dz)
  if n_params() ne 3 then begin
     print,'Syntax - sdss_calcdndx(civstrct_fil, cmplt_fil, ewlim, [/final, /dz, '
     print,'                        /silent, /afp_flg, /czn, _extra=])'
     return,-1
  endif 

  ;; Read in data
  if size(cmplt_fil,/type) eq 7 then cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil
  tags = tag_names(cmpltstr)
  if keyword_set(dz) then begin
     xstr = 'z' 
     dxtag = (where(tags eq 'DZ_2DARR'))[0]
  endif else begin
     xstr = 'X'
     dxtag = (where(tags eq 'DX_2DARR'))[0]
  endelse 
  index = lindgen(cmpltstr.newbin)


  if size(civstrct_fil,/type) eq 8 then civstr = civstrct_fil $
  else civstr = xmrdfits(civstrct_fil,1,/silent)
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endelse 

  ;; May need column densities so just compile
  dblt = dblt_retrieve(civstr[0].wrest[0])
  ;; _extra incldues /estflg
  ncolm = sdss_calcncolmdblt(civstr,dblt,signcolm=signcolm,$
                            ncolmflg=ncolmflg,final=final,/silent,$
                            /log,_extra=extra)
  ncolmclim = sdss_getewclim(cmpltstr,/ncolm,dblt_name=dblt,/log) ; trim below

  ;; Set up the data
  if keyword_set(czn) then begin
     data = ncolm
     sigdata = signcolm
     dataflg = ncolmflg
  endif else begin
     data = civstr.(ewtag)[0]
     sigdata = civstr.(sigewtag)[0]
     dataflg = civstr.ewflg[0]
  endelse 

  dndx_strct = {$
               ewlim:ewlim, $
               dz:keyword_set(dz), $
               afp_flg:keyword_set(afp_flg), $
               czn:keyword_set(czn), $
               czn_cmplt:keyword_set(cmpltstr.czn), $
               zcenter:fltarr(cmpltstr.nzbin,/nozero), $
               sigzcenter:fltarr(cmpltstr.nzbin,2,/nozero), $
               zmed:fltarr(cmpltstr.nzbin,/nozero), $ ; from discrete data (cmplt weighted)
               sigzmed:fltarr(cmpltstr.nzbin,2,/nozero), $
               numtot0:fltarr(cmpltstr.nzbin,/nozero), $
               signumtot0:fltarr(cmpltstr.nzbin,2,/nozero), $
               numtot:fltarr(cmpltstr.nzbin,/nozero), $
               signumtot:fltarr(cmpltstr.nzbin,2,/nozero), $
               dxtot:fltarr(cmpltstr.nzbin,/nozero), $
               sigdxtot:fltarr(cmpltstr.nzbin,2,/nozero), $
               dndx:fltarr(cmpltstr.nzbin,/nozero), $
               sigdndx:fltarr(cmpltstr.nzbin,2,/nozero) $
               }
  if cmpltstr.nzbin gt 1 then begin
     dndx_strct.zcenter = (0.5*(cmpltstr.zlim[*,0] + cmpltstr.zlim[*,1])) 
     dndx_strct.sigzcenter[*,0] = dndx_strct.zcenter - cmpltstr.zlim[*,0]
     dndx_strct.sigzcenter[*,1] = cmpltstr.zlim[*,1] - dndx_strct.zcenter
  endif else begin
     dndx_strct.zcenter = mean(cmpltstr.zlim) 
     dndx_strct.sigzcenter[*,0] = dndx_strct.zcenter - cmpltstr.zlim[0]
     dndx_strct.sigzcenter[*,1] = cmpltstr.zlim[1] - dndx_strct.zcenter
  endelse 


  ;; Find what's in bins
  sub = where(data ge ewlim[0] and $
              data lt ewlim[1],nciv)
  
  if sub[0] eq -1 then begin
     print,'sdss_calcdndx(): no doublets with desired EW'
     return,-1
  endif 
  civstr = civstr[sub]          ; just truncate
  data = data[sub]
  sigdata = sigdata[sub]
  dataflg = dataflg[sub]
  ncolm = ncolm[sub]
  signcolm = signcolm[sub]
  ncolmflg = ncolmflg[sub]
  

  ;; Get data necessary
  ;; _extra includes /grid, /ewmax, /lsquadratic,
  ;; /quadratic, /spline 
  ;; FORCE civobs_corr=0 because pathlength should already be adjusted
  ;; for part blocked by other absorbers in z and EW bin
  if keyword_set(czn) ne keyword_set(cmpltstr.czn) then begin
     ;; Then there could be a problem
     print,'sdss_calcdndx(): WARNING! Mixing EW and column density between dN/dX and completeness'
     if keyword_set(cmpltstr.czn) then $ ; query completeness in column density space
        cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                              ncolm, sigewabs=signcolm, $
                              sigma=sigcxw_civ, dz=dz, /czw, $
                              civobs_corr=0, final=final, $
                              dxwmax=dxwmax, _extra=extra) $
     else $                     ; then query completeness in EW space
        cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                              civstr.(ewtag)[0], sigewabs=civstr.(sigewtag)[0], $
                              sigma=sigcxw_civ, dz=dz, /czw, $
                              civobs_corr=0, final=final, $
                              dxwmax=dxwmax, _extra=extra) 
  endif else $                  ; then completeness and data are consistent
     cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                           data, sigewabs=sigdata, $
                           sigma=sigcxw_civ, dz=dz, /czw, $
                           civobs_corr=0, final=final, $
                           dxwmax=dxwmax, _extra=extra) 


  ;; Loop through redshift bins
  for zz=0,cmpltstr.nzbin-1 do begin
     zlim = dndx_strct.zcenter[zz]+$
            [-dndx_strct.sigzcenter[zz,0],dndx_strct.sigzcenter[zz,1]]

     gd = where(civstr.(ztag)[0] ge zlim[0] and civstr.(ztag)[0] lt zlim[1],ngd)

     ;; Save all values generally
     dndx_strct.numtot0[zz] = ngd
     ;; _extra includes sigma=, cl=,
     dndx_strct.signumtot0[zz,*] = $
        sdss_calcsigpoiss(float(ngd),verbose=(not keyword_set(silent)),$
                          _extra=extra) ; included in sdss_calcncmplt() now
     if ngd eq 0 then begin
        ;; must instantiate
        dndx_strct.zmed[zz] = dndx_strct.zcenter[zz] ; default 
        dndx_strct.sigzmed[zz,*] = dndx_strct.sigzmed[zz,*] 
        dndx_strct.numtot[zz] = dndx_strct.numtot0[zz]
        dndx_strct.signumtot[zz,*] = dndx_strct.signumtot0[zz,*]
        dndx_strct.dxtot[zz] = dxwmax[zz] 
        dndx_strct.sigdxtot[zz,*] = 0.
        dndx_strct.dndx[zz] = 0. 
        dndx_strct.sigdndx[zz,0] = 0. 
        ;; depending on the upper limit on detection divided by the
        ;; pathlength 
        dndx_strct.sigdndx[zz,1] = dndx_strct.signumtot[zz,1] / $
                                   dndx_strct.dxtot[zz] 
     endif else begin
        ;; Weighted redshift
        dndx_strct.zmed[zz] = sdss_medianw(civstr[gd].(ztag)[0],1/cxw_civ[gd],/even)
        dndx_strct.sigzmed[zz,0] = dndx_strct.zmed[zz] - $
                                   (dndx_strct.zcenter[zz] - $
                                    dndx_strct.sigzcenter[zz,0]) ; handles best
        dndx_strct.sigzmed[zz,1] = (dndx_strct.zcenter[zz] + $
                                    dndx_strct.sigzcenter[zz,1]) - $
                                   dndx_strct.zmed[zz] 
                                   

        ;; Using completeness in numerator
        ;; Completeness-corrected number which includes lump Poisson error
        dndx_strct.numtot[zz] = sdss_calcncmplt(cxw_civ[gd],sigcxw_civ[gd,*],$
                                                signcmplt=signcmplt)

        ;; Errors just in N = sum(1 / C(W))
        dndx_strct.signumtot[zz,*] = signcmplt

        ;; Total dX(W) is the max possible in bin (not necessarily the
        ;; 100% complete limit because may not exist).
        ;; Add back in the contribution of each absorber in the bin.
        dndx_strct.sigdxtot[zz,*] = 0.
        dndx_strct.dxtot[zz] = dxwmax[zz] ;+ dndx_strct.sigdxtot[zz,0] OBSOLETE
        
        dndx_strct.dndx[zz] = dndx_strct.numtot[zz] / dndx_strct.dxtot[zz]
        
        ;; Error combining counting error for what's in bin (observed error)
        dndx_strct.sigdndx[zz,*] = dndx_strct.signumtot[zz,*] / dndx_strct.dxtot[zz]

     endelse 

     ;; Accepted false-positive adjustment
     if keyword_set(afp_flg) then begin
        if keyword_set(czn) or keyword_set(cmpltstr.czn) then $
           stop,'sdss_calcdndx() stop: accepted false-positive adjustment does not work for columns'

       ;; Correct a given dN/dX structure for accepted false positive
        ;; dN/dX_afp 
        ;; dN/dX = dN/dX_0 - dN/dX_afp
        ;; sig(dN/dX)^2 = sqrt( sig(dN/dX_0)^2 - sig(dN/dX_afp)^2 )
        ;; dN/dX structure is organized so that there are multiple
        ;; redshift bins which must be handled indepedently

        ;; For now, going to read in the file every loop
        ;; afndndx will only be for one bin
        ;; _extra includes stuff for sdss_calcdxzw()
        afpdndx = sdss_getafpdndx(zlim, dndx_strct.ewlim, dz=dndx_strct.dz, czn=czn, _extra=extra)
        
        nwdndx_strct = dndx_strct
        nwdndx_strct.numtot0[zz] = dndx_strct.numtot0[zz] - afpdndx.numtot0
        nwdndx_strct.signumtot0[zz,*] = sqrt(dndx_strct.signumtot0[zz,*]^2 + $
                                             afpdndx.signumtot0[0,*]^2)
        nwdndx_strct.numtot[zz] = dndx_strct.numtot[zz] - afpdndx.numtot
        nwdndx_strct.signumtot[zz,*] = sqrt(dndx_strct.signumtot[zz,*]^2 + $
                                            afpdndx.signumtot[0,*]^2)
        ;; dxtot, sigdxtot do *not* change
        nwdndx_strct.sigdndx[zz,*] = sqrt(dndx_strct.sigdndx[zz,*]^2 + $
                                          afpdndx.sigdndx[0,*]^2)
        if dndx_strct.dndx[zz] eq 0. then $
           ;; Upper limit (just inflating the error)
           nwdndx_strct.sigdndx[zz,0] = 0. $
        else nwdndx_strct.dndx[zz] = dndx_strct.dndx[zz] - afpdndx.dndx

        dndx_strct = nwdndx_strct
     endif                      ; /afp_flg
     
  endfor                        ; loop zz=cmpltstr.nzbin

  if not keyword_set(silent) then $
     sdss_prntdndx, dndx_strct, _extra=extra ; includes /skip_null, out_fil=

  return, dndx_strct
end                             ; sdss_calcdndx()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntfxw, fxw_fil, skip_null=skip_null, out_fil=out_fil, $
                  multiext=multiext
  if n_params() ne 1 then begin
     print,'Syntax - sdss_prntfxw, fxw_fil, [/skip_null]'
     return
  endif 

  if size(fxw_fil,/type) eq 7 then $
     fxw = xmrdfits(fxw_fil,1,/silent) $
  else fxw = fxw_fil
  if keyword_set(multiext) then begin
     ss = 1L
     if size(fxw_fil,/type) ne 7 then $
        stop,'sdss_prntfxw stop: input must be file to access multiple extensions.'
  endif else begin
     nstrct = (size(fxw,/dim))[0] > 1
     ss = 0L
  endelse 

  if keyword_set(out_fil) then $
     openw,1,out_fil
  
  done = 0
  while not done do begin
     if keyword_set(multiext) then begin
        subfxw = fxw
     endif else begin
        subfxw = fxw[ss]
     endelse 

     if keyword_set(subfxw.czn) then begin
        varstr = 'N' 
        ;; put in log space
        gd = where(subfxw.fxw ne 0.,complement=bd)
        if gd[0] ne -1 then begin
           subfxw.sigfxw[gd,0] = abs(subfxw.sigfxw[gd,0]/$
                                      (alog(10.)*subfxw.fxw[gd]))
           subfxw.sigfxw[gd,1] = abs(subfxw.sigfxw[gd,1]/$
                                      (alog(10.)*subfxw.fxw[gd]))
           subfxw.fxw[gd] = alog10(subfxw.fxw[gd])
        endif 
        if bd[0] ne -1 then begin
           subfxw.sigfxw[bd,1] = alog10(subfxw.sigfxw[bd,1])
           ;; Rest remain zero
        endif 
     endif else begin
        varstr = 'W'
     endelse 

     lin0 = strtrim(string(subfxw.zlim[0],format='(f7.4)'),2)+' <= z < '+$
            strtrim(string(subfxw.zlim[1],format='(f7.4)'),2)
     if subfxw.dz eq 1 then colhead = ['dz','f(z,'+varstr+')'] $
     else colhead = ['dX','f(X,'+varstr+')']
     lin = string(varstr+'lo',varstr+'hi','Num0','Errlo','Errhi',$
                  'Num','Errlo','Errhi',$
                  colhead[0],'Errlo','Errhi',colhead[1],'Errlo','Errhi',$
                  format='(2(a6,1x),3(a5,1x),3(a8,1x),3(a10,1x),3(a8,1x))')
     if keyword_set(out_fil) then begin
        printf,1,''
        printf,1,lin0
        printf,1,lin
     endif else begin
        print,''
        print,lin0
        print,lin
     endelse

     newbin = (size(subfxw.ewcenter,/dim))[0] > 1
     for ee=0,newbin-1 do begin
        if keyword_set(skip_null) and subfxw.fxw[ee] eq 0. then $
           continue 
        lin = string(subfxw.ewcenter[ee]-subfxw.sigewcenter[ee,0],$
                     subfxw.ewcenter[ee]+subfxw.sigewcenter[ee,1],$
                     subfxw.numtot0[ee],subfxw.signumtot0[ee,0],$
                     subfxw.signumtot0[ee,1],subfxw.numtot[ee],$
                     subfxw.signumtot[ee,0],subfxw.signumtot[ee,1],$
                     subfxw.dxtot[ee],$
                     subfxw.sigdxtot[ee,0],subfxw.sigdxtot[ee,1],$
                     subfxw.fxw[ee],subfxw.sigfxw[ee,0],$
                     subfxw.sigfxw[ee,1],$
                     format='(2(f6.2,1x),i5,1x,2(f5.1,1x),f8.1,1x,2(f8.1,1x),3(f10.2,1x),3(f8.4,1x))')
        if keyword_set(out_fil) then printf,1,lin $
        else print,lin
     endfor                     ; loop ee=newbin

     ss++
     if keyword_set(multiext) then begin
        fxw = xmrdfits(fxw_fil,ss,/silent)
        if size(fxw,/type) ne 8 then done = 1
     endif else begin
        if ss gt nstrct-1 then done = 1 
     endelse
  endwhile                      ; not done

  if keyword_set(out_fil) then begin
     close,1
     print,'sdss_prntfxw: created ',out_fil
  endif

end                             ; sdss_prntfxw



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcfxwfit, data, option=option, fitstrct_fil=fitstrct_fil, $
                          coeff=coeff, alpha=alpha, datanorm=datanorm, ndim=ndim
  if n_params() ne 1 then begin
     print,'Syntax - sdss_calcfxwfit(data, [option=, fitstrct_fil=, coeff=, alpha=, datanorm=, ndim=])'
     return,-1
  endif 

   if keyword_set(fitstrct_fil) then begin
      if size(fitstrct_fil,/type) eq 8 then $
         fitstrct = fitstrct_fil $
      else fitstrct = xmrdfits(fitstrct_fil,1,/silent)
      option = fitstrct.fittype
      coeff = fitstrct.coeff
      alpha = fitstrct.alpha
      datanorm = fitstrct.datanorm
   endif else if not keyword_set(option) then $
      stop,'sdss_calcfxwfit() stop: option keyword not set'
   
   ndim = 2
   case strlowcase(option) of
      'pow': rslt = (data/datanorm)^alpha 
      'exp': rslt = exp(data/datanorm * alpha)  
      'sch': begin
         rslt = (data/datanorm)^alpha *exp(-data/datanorm) 
         ndim = 3
      end
   endcase
   rslt = coeff * rslt
   
   return, rslt
end                             ; sdss_calcfxwfit()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcncolmfit, data, cmplt_fil, flgdata=flgdata, sigdata=sigdata, $
                            seed=seed, oseed=oseed, intlim=intlim, debug=debug, _extra=extra

  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcncolmfit(data, cmplt_fil, [flgdata=, sigdata=, intlim=, seed=, oseed=, /debug, _extra=])'
     return,-1
  endif 

  ndata = (size(data,/dim))[0] > 1

  ;; completeness better be column density
  if size(cmplt_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil
  
  if not keyword_set(intlim) then intlim = [min(data),10.^18]
  if keyword_set(flgdata) then begin
     if not keyword_set(sigdata) then $
        stop,'sdss_calcncolmfit() stop: must set sigdata to use flgdata' 
  endif else $
     flgdata = replicate(sdss_getlimflg(/lower),ndata) ; assume all lower limits
  
  ;; Define range
  data_range = [min(data,max=mx),mx]

  ;; Set up column density array for numerical integration
  tmp = sdss_mkewarr(alog10([data_range[1],intlim[1]]),cmpltstr.ewbinsize,$
                     newbin=ntmp)
  ncolm_global = [data,10.^tmp[1:*]] ; first should be input and last one
  nncolm_global = ndata + ntmp - 1
  dncolm_global = shift(ncolm_global,-1) - ncolm_global
  dncolm_global[nncolm_global-1] = 10.^(tmp[ntmp-1]+cmpltstr.ewbinsize)-ncolm_global[nncolm_global-1]
  if keyword_set(sigdata) then $
     signcolm_global = [sigdata,cmpltstr.ewbinsize/(alog(10)*10.^tmp[1:*])] $
  else signcolm_global = 0
  z_global = replicate(mean(cmpltstr.zlim),nncolm_global)

  dxw_global = sdss_getdxw(cmpltstr, z_global, alog10(ncolm_global), $
                           sigewabs=signcolm_global/(alog(10)*ncolm_global),$
                           sigma=sigdxw_global,dxwmax=dxwmax)

  ;; assuming Schechter function with given fit_params
  ;; f(N) = k (N/N*)^alpha exp(-N/N*).
  ;; f(N) < f(N_lim) due to monotonic nature of Schechter function.

  ;; Cumulative f(N)*dN*dX give probability for drawing actual column
  ;; for a limit, under the assumption of the Schechter model. Going
  ;; to create the CDF empirically, once, and the scale as needed for
  ;; each limit.
  ;; Number(i) = Number(i-1) + f(N_i)*(N_i-N_i-1)*dX(N_i)
  ;; _extra= includes option=, fitstrct_fil=, coeff=, alpha=, datanorm=, ndim=
  fN_global = sdss_calcfxwfit(ncolm_global,_extra=extra)
  cdf_global = total(fN_global*dncolm_global*dxw_global,/cum)
  rnum = randomu(seed,ndata)
  new_ncolm_civ = dblarr(ndata,/nozero)
  for ii=0L,ndata-1 do begin
     if (flgdata[ii] and sdss_getlimflg(/lower)) eq sdss_getlimflg(/lower) then begin
        ;; Customized normalization
        cdf_tmp = (cdf_global[ii:*]-cdf_global[ii])/(cdf_global[nncolm_global-1]-cdf_global[ii])
        new_ncolm_civ[ii] = interpol(ncolm_global[ii:*], cdf_tmp, rnum[ii], _extra=extra)
     endif else $
        new_ncolm_civ[ii] = data[ii] + randomn(seed,1)*sigdata
  endfor                        ; loop ii=ndata

  if keyword_set(debug) then $
     x_splot,ncolm_civ,new_ncolm_civ,psym1=4,xtitle='Column Density Limit',$
             ytitle='Expected Column Density Drawn from Assumed f(N)',$
             xtwo=[0,1e20],ytwo=[0,1e20],/block

  oseed = seed[0]

  return, new_ncolm_civ
end                             ; sdss_calcncolmfit()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcdndxfit, ewlim, option=option, fitstrct_fil=fitstrct_fil, $
                           coeff=coeff, alpha=alpha, datanorm=datanorm, ndim=ndim
  if n_params() ne 1 then begin
     print,'Syntax - sdss_calcdndxfit(data, option,[fitstrct_fil=, coeff=, alpha=, datanorm=])'
     return,-1
  endif 

   if keyword_set(fitstrct_fil) then begin
      if size(fitstrct_fil,/type) eq 8 then $
         fitstrct = fitstrct_fil $
      else fitstrct = xmrdfits(fitstrct_fil,1,/silent)
      option = fitstrct.fittype
      coeff = fitstrct.coeff
      alpha = fitstrct.alpha
      datanorm = fitstrct.datanorm
   endif else if not keyword_set(option) then $
      stop,'sdss_calcdndxfit() stop: option keyword not set'
   if n_elements(ewlim) eq 1 then ewrng = [ewlim,!values.f_infinity] $
   else ewrng = ewlim           ; assume
   if size(ewrng,/n_dim) eq 1 then begin
      ewrng = transpose(ewrng)  ; access [*,0:1]
      flg_ewrng = 1
   endif 


   ndim = 2
   case strlowcase(option) of
      'pow': rslt = 1./(1+alpha)*(ewrng[*,1]^(1+alpha)-ewrng[*,0]^(1+alpha))/$
                    datanorm^alpha
      'exp': rslt = datanorm/alpha*(exp(ewrng[*,1]/datanorm * alpha) $
                                    - exp(ewrng[*,0]/datanorm * alpha))
      'sch': begin
         ;; incomplete gamma function igamma(2+alpha, data/datanorm)
         ;; dN/dX = -coeff*datanorm*(igamma(2+alpha,ewrng[*,1]/datanorm)
         ;; - igamma(2+alpha,ewrng[*,0]/datanorm)
         if finite(ewrng[*,1]) eq 0 then $
            rslt = datanorm*igamma(2+alpha, ewrng[*,0]/datanorm,/double) $
         else $
            rslt = -1.*datanorm*($
                   igamma(2+alpha, ewrng[*,1]/datanorm,/double) - $
                   igamma(2+alpha, ewrng[*,0]/datanorm,/double)) 
         ndim = 3
      end
   endcase

   if keyword_set(flg_ewrng) then rslt = rslt[0] ; expected dimensionality

   rslt = coeff * rslt
   
   return, rslt
end                             ; sdss_calcdndxfit()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_expgrowthfunc, var, coeff, func, pder
  ;; C(W) = C0*( 1-exp(C1 * (W - C2)) )
  ncoeff = (size(coeff,/dim))[0] > 1

  if ncoeff eq 3 then nwvar = var - coeff[2] 

  func = coeff[0] * (1. - exp(coeff[1]*nwvar))
  
  if n_params() eq 4 then begin
     if ncoeff eq 3 then $
        pder = [$
               [1-exp(coeff[1]*nwvar)], $
               [-1.*coeff[0]*nwvar*exp(coeff[1]*nwvar)], $
               [coeff[0]*coeff[1]*exp(coeff[1]*nwvar)] $
               ] $
     else $
        pder = [ $
               [1-exp(coeff[1]*var)], $
               [-1.*coeff[0]*nwvar*exp(coeff[1]*var)] $
               ]
  endif 
  
end                             ; sdss_expgrowthfunc


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcexpgrowthfunc, xfit, coeff, sigma=sigma
  ;; Actually compute the fit
  ;; C(W) = C0*(1-exp(C1*(W-C2)))
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcexpgrowthfunc(xfit, coeff, [sigma=])'
     return,-1
  endif 
  ncoeff = (size(coeff[*,0],/dim))[0] > 1
  sdss_expgrowthfunc, xfit, coeff[*,0], yfit, pder
  nfit = (size(xfit,/dim))[0] > 1

  ;; Error: sigma^2 = sigC0^2*(1-exp(C1*(W-C2)))^2 +
  ;; sigC1^2*(C0*(W-C2)*-exp(C1*(W-C2)) +
  ;; sigC2^2*(C0*C1*exp(C1*(W-C2)))^2
  ;; 
  sigma = fltarr(nfit,ncoeff)  
  sigma[*,0] = coeff[0,1]^2 * pder[*,0]^2 + coeff[1,1]^2 * pder[*,1]^2 + $
               2*coeff[0,1]*coeff[1,1]*pder[*,0]*pder[*,1]
  sigma[*,1] = coeff[0,2]^2 * pder[*,0]^2 + coeff[1,2]^2 * pder[*,1]^2 + $
               2*coeff[0,2]*coeff[1,2]*pder[*,0]*pder[*,1]           
  if ncoeff eq 3 then begin
     sigma[*,0] = sigma[*,0] + coeff[2,1]^2 * pder[*,2]^2 + $
                  2*coeff[2,1]*coeff[1,1]*pder[*,2]*pder[*,1] + $
                  2*coeff[2,1]*coeff[0,1]*pder[*,2]*pder[*,0]
     sigma[*,1] = sigma[*,1] + coeff[2,2]^2 * pder[*,2]^2 + $
                  2*coeff[2,2]*coeff[1,2]*pder[*,2]*pder[*,1] + $
                  2*coeff[2,2]*coeff[0,2]*pder[*,2]*pder[*,0]
  endif 
  sigma = sqrt(sigma)

  return, yfit
end                             ; sdss_calcexpgrowthfunc()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getcivstrct, strct_fil, ext, dblt_name=dblt_name, dvem=dvem, $
                           rating=rating, zlim=zlim, ewlim=ewlim, $
                           nlim=nlim, dvgal=dvgal, dvqso=dvqso, noBAL=noBAL, $
                           unblend=unblend, final=final, default=default, $
                           dropbox=dropbox, civstrct_fil=civstrct_fil, $
                           count=count, help=help, _extra=extra
  if keyword_set(help) then begin
     print,'Syntax - sdss_getcivstrct([strct_fil, ext, dblt_name=, rating=, zlim=, ewlim=, '
     print,'                           nlim=, dvem=, dvgal=, dvqso=, /noBAL,/unblend, /final, /dropbox,'
     print,'                           /default, civstrct_fil=, count=, /help, _extra=])'
     return,-1
  endif 

  sdssdir = sdss_getsdssdir()
  if keyword_set(default) then begin
     ext = 1 ; fits file extension
     ;; Force rating = 2, 3, none in visual BALs
     if keyword_set(rating) then begin
        print,'WARNING!!! '
        print,'sdss_getcivstrct(): Asking for /default but overriding rating with',rating
        print,'WARNING!!! '
     endif else rating = [2,3]
     noBAL = 1
     if keyword_set(dvqso) then begin
        print,'WARNING!!! '
        print,'sdss_getcivstrct(): Asking for /default but overriding dvqso with',dvqso 
        print,'WARNING!!! '
     endif else dvqso = -5000.  ; km/s
  endif 
 
  if not keyword_set(ext) then ext = 1 $ ; default 
  else begin
     if ext lt 1 then ext = 0  ; only way to handle wanting zeroth extension
  endelse 

  c = 299792.458                ; km s^-1
  if keyword_set(strct_fil) then civstrct_fil = strct_fil $
  else begin
     if keyword_set(dblt_name) then begin
        if size(dblt_name,/type) eq 8 then dblt = dblt_name $
        else dblt = dblt_retrieve(dblt_name)
     endif else dblt = dblt_retrieve('')
     case dblt.ion of 
        'MgII': $
           if keyword_set(dropbox) then $
              civstrct_fil = getenv('HOME')+$
                             '/Dropbox/CIV/MIT/PaperMgII/sdss_mgiirate_hyb_all_rtgge2.fit' $
           else civstrct_fil = sdssdir+'candidates/rated/sdss_mgiirate_hyb_all_rtgge2.fit'
        else: $                 ; CIV always default
           if keyword_set(dropbox) then $
              civstrct_fil = getenv('HOME')+$
                             '/Dropbox/CIV/MIT/Paper1/sdss_civrate_hyb_all.fit' $
           else civstrct_fil = sdssdir+'candidates/rated/sdss_civrate_hyb_all.fit'
     endcase
  endelse 
  if size(civstrct_fil,/type) eq 7 then $
     substrct = xmrdfits(civstrct_fil,ext,/silent) $
  else substrct = civstrct_fil
  count = (size(substrct,/dim))[0] > 1
  
  tags = tag_names(substrct)
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     ncolmtag = (where(tags eq 'NCOLM_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     ncolmtag = (where(tags eq 'NCOLM_ORIG'))[0]
  endelse 

  ;; ;;;;;;;
  if keyword_set(rating) then begin
     nrtg = n_elements(rating)
     substrct0 = substrct
     for rr=0,nrtg-1 do begin
        gd = where(rating[rr] eq substrct0.rating[0],ngd)
        if ngd eq 0 then begin
           print,'sdss_getcivstrct(): no doublets with rating = '+$
                 strtrim(rating[rr],2)
           done = 0
        endif else begin
           if not keyword_set(done) then begin
              substrct = substrct0[gd]
              count = ngd
              done = 1          ; don't come here again
           endif else begin
              substrct = [substrct,substrct0[gd]]
              count = count + ngd
           endelse 
        endelse 
     endfor                     ; rr = nrtg
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets with rating(s) = ',rating
        return, -1              ; EXIT
     endif 
  endif                         ; rating=

  ;; ;;;;;;;
  if keyword_set(zlim) then begin
     gd = where(substrct.(ztag)[0] ge zlim[0] and substrct.(ztag)[0] lt zlim[1],$
                count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets within redshift limits ',zlim
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse 
  endif                         ; zlim=

  ;; ;;;;;;;
  if keyword_set(ewlim) then begin
     gd = where(substrct.(ewtag)[0] ge ewlim[0] and $
                substrct.(ewtag)[0] lt ewlim[1],count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets within EW limits ',ewlim
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse 
  endif                         ; ewlim=
        
  ;; ;;;;;;;
  if keyword_set(nlim) then begin
     dblt = dblt_retrieve(substrct[0].wrest[0])
     if dblt.ion eq '' then begin
        print,'sdss_getcivstrct(): NOTE!!! applying nlim= to 0th ion'
        gd = where(substrct.(ncolmtag)[0] ge nlim[0] and $
                   substrct.(ncolmtag)[0] lt nlim[1],count) 
     endif else begin 
        ;; _extra= includes /eq_weight, /silent, /debug
        ;; just going to assume if nlim[0] is large (> 100),
        ;; it's linear and doing that for the lower limit
        ;; because the upper one likely infinity.
        ncolm = sdss_calcncolmdblt(substrct,dblt,final=final,$
                                   log=(nlim[0] lt 100.),_extra=extra)
        gd = where(ncolm ge nlim[0] and ncolm lt nlim[1],count)
     endelse
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets within column limits ',nlim
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse 
  endif                         ; nlim=

  ;; ;;;;;;;
  if keyword_set(noBAL) then begin
;     gd = where((substrct.balflg and sdss_getbalflg(/visual)) eq 0,count)
     gd = where(substrct.balflg eq 0,count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets not in BAL sightlines '
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse      
  endif                         ; /noBAL

 ;; ;;;;;;;
  if keyword_set(unblend) then begin
     gd = where(substrct.rating[9] eq 0,count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no un-blended doublets'
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse      
  endif                         ; /unblend

  ;; ;;;;;;;
  if keyword_set(dvgal) or keyword_set(dvqso) then begin
     ;; Assume Lya-forest or OI/SiII forest limits already handled
     if not keyword_set(dvgal) then dvgal = 0.
     if not keyword_set(dvqso) then dvqso = 0.
     dvqso_civ = (substrct.(ztag)[0] - substrct.z_qso)/(1. + substrct.z_qso) * c
     dvgal_civ = substrct.(ztag)[0] * c
     gd = where(dvgal_civ gt dvgal and dvqso_civ lt dvqso,count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets > dvgal, < dvqso ',dvgal, dvqso
        return, -1
     endif else begin
        substrct = substrct[gd]
     endelse 
  endif                         ; dvgal= or dvqso=

  ;; ;;;;;;;
  if keyword_set(dvem) then begin
     ;; Emission-line rest wavelength; dv range
     ;; either [wrest, dvlo, dvhi] or 
     ;; transpose([[wrest1,dvlo1,dvhi1],[wrest2,dvlo2,dvhi2]])
     n_em = size(dvem,/n_dim)
     if n_em eq 1 then dvem_use = transpose(dvem) $ ; [1,3]
     else dvem_use = dvem                           ; [n_em,3]
     mask = intarr(count)
     for ii=0,n_em-1 do begin
        zem_civ = dvem_use[ii,0]/substrct.wrest[0]*(1+substrct.z_qso) - 1.
        dvem_civ = c*(substrct.zabs_orig[0]-zem_civ)/(1+zem_civ)
        bd = where(dvem_civ ge dvem_use[ii,1] and dvem_civ lt dvem_use[ii,2],nbd)
        if nbd ne 0 then mask[bd]++
     endfor                     ; loop ii=n_em
     gd = where(mask eq 0,count)
     if count eq 0 then begin
        print,'sdss_getcivstrct(): no doublets outside dvem='
        printcol,dvem_use[*,0],dvem_use[*,1],dvem_use[*,2]
        return,-1
     endif else begin
        substrct = substrct[gd]
     endelse 
  endif                         ; dvem=

  ;; Finite values only
  gd = where(finite(substrct.ew_orig[0]) eq 1,count)
  substrct = sdss_srtcivstrct(substrct[gd])
  return, substrct
end                             ; sdss_getcivstrct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_estczn, cmplt_fil, civstrct_fil, nbinsize=nbinsize, $
                      ewlim=ewlim, nlim=nlim, final=final, _extra=extra
  ;; Estimate C(N) from C(W) by drawing the C(W) directly with the
  ;; given data
  if n_params() ne 2 then begin
     print,'Syntax - sdss_estczn(cmplt_fil, civstrct_fil, [nbinsize=, ewlim=, nlim=, /final, _extra=])'
     return,-1
  endif
  
  ;; _extra= includes rating=, zlim=, dvgal=, dvqso=, /noBAL,
  ;; /unblend, /final, /default, /dropbox, civstrct_fil=, /help,
  ;; _extra=
  civstrct = sdss_getcivstrct(civstrct_fil, ewlim=ewlim, nlim=nlim, $
                              count=nciv, final=final, _extra=extra)
  tags = tag_names(civstrct)
  if keyword_set(final) then begin
     ;; Use final zabs and ew
     ztag = where(tags eq 'ZABS_FINAL')
     ewtag = where(tags eq 'EW_FINAL')
     sigewtag = where(tags eq 'SIGEW_FINAL')
  endif else begin
     ;; Use original zabs and ew
     ztag = where(tags eq 'ZABS_ORIG')
     ewtag = where(tags eq 'EW_ORIG')
     sigewtag = where(tags eq 'SIGEW_ORIG')
  endelse 
  
  if size(cmplt_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil
  if cmpltstr.nzbin ne 1 then begin
     print,'sdss_estczn(): WARNING! cannot work on mulit-z completeness curves'
     return,-1
  endif 
  ew_global = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  
  if not keyword_set(nbinsize) then nbinsize = 0.1 
  if keyword_set(ewlim) then ewseam = ewlim[0] $
  else ewseam = min(civstrct.(ewtag)[0])

  ;; _extra includes /estflg,/debug
  dblt = dblt_retrieve(civstrct[0].wrest[0])
  ncolm = sdss_calcncolmdblt(civstrct,dblt,signcolm=signcolm,/log,$
                             ncolmflg=ncolmflg,/silent,final=final,$
                             _extra=extra)

  ;; Need the completeness 
  ;; _extra includes /dz, /grid, ewmax=, civobs_corr=, dxwmax=, dxciv= 
  czw = sdss_getdxw(cmpltstr, civstrct.(ztag)[0], civstrct.(ewtag)[0], $
                    sigewabs=civstrct.(sigewtag)[0], sigma=sigczw, /czw, $
                    /silent, _extra=extra)

  ;; Round the binning down just to make it even
  if keyword_set(nlim) then ncolmlim = nlim $ ; may be changed
  else ncolmlim = [min(ncolm,max=mx),mx] 
  ncolmlim[0] = floor(ncolmlim[0]/nbinsize)*nbinsize
  ncolm_global = sdss_mkewarr(ncolmlim, nbinsize, newbin=nnbin)
  
  ;; ;;;;;;;
  ;; Create output (matches sdss_completeness_czw)
  f2darr = fltarr(nnbin,/nozero)
  cznstr = { $
           ;; General parameters, so can re-create this
           list_fil:cmpltstr.list_fil, $
           dblt_name:cmpltstr.dblt_name, $
           rec_param:cmpltstr.rec_param, $
           civobs_corr:cmpltstr.civobs_corr, $
           userbias:cmpltstr.userbias, $
           czn:1, $             ; by construction
           ewseam:ewseam, $
           nzbin:cmpltstr.nzbin, $
           zlim:cmpltstr.zlim, $
           zbinsize:cmpltstr.zbinsize, $
           newbin:nnbin, $
           ewlim:ncolmlim, $
           ewbinsize:nbinsize, $
           cosmology:cmpltstr.cosmology, $
           ;; Now for actual results in this bin
           ninput_2darr:f2darr, $
           nrec_2darr:f2darr, $
           rz_2darr:f2darr, $   ; can reform(2darr, nzbin, newbin)
           dz_2darr:f2darr, $
           dx_2darr:f2darr, $
           czw_2darr:f2darr, $
           sigczw_2darr:fltarr(nnbin,2,/nozero) $ ; [*,0] = lower, [*,1] = upper
           }

  ;; Loop!
  for nn=0,cznstr.newbin-1 do begin
     gd = where(ncolm ge ncolm_global[nn] and $
                ncolm lt ncolm_global[nn]+cmpltstr.ewbinsize,ngd)
     cznstr.ninput_2darr[nn] = ngd
     case ngd of 
        0: begin
           ;; Must instantiate
           cznstr.czw_2darr[nn] = !values.f_nan ; interpolate later
        end
        1: begin
           cznstr.czw_2darr[nn] = czw[gd[0]]
           cznstr.sigczw_2darr[nn,*] = sigczw[gd[0],*]
        end
        else: begin
           ;; Fundamentally, the best f(N) to be made from C(W) is
           ;; f(N) = sum(1/C(W))/dN/dXtot 
           ;; since sum(1/C(W))/dW/dXtot ~ Nobs/dW/<dXtot*C(W)>
           ;; I can prove that the best approximation for C(N) is
           ;; <C(N)> ~ Nobs/sum(1/C(W))... and the square is just a
           ;; stronger weighting
           cznstr.czw_2darr[nn] = sqrt(ngd/total(1/czw[gd]^2))
           ;; Take weighted sum and something representing 1 sigma
;           cznstr.czw_2darr[nn] = $
;              total(czw[gd]/(0.5*(sigczw[gd,0]+sigczw[gd,1]))^2)/$
;              total(1./(0.5*(sigczw[gd,0]+sigczw[gd,1]))^2)
;              total(czw[gd]*(1./sigczw[gd,0]^2+1./sigczw[gd,1]^2))/$
;              total(1./sigczw[gd,0]^2+1./sigczw[gd,1]^2)
           ;; Take median and something representing 1 sigma
;           cznstr.czw_2darr[nn] = gmean(czw[gd])
;           cznstr.czw_2darr[nn] = median(czw[gd],/even)
;           cznstr.czw_2darr[nn] = sdss_minwad([czw[gd],czw[gd]],$
;                                          [1./sigczw[gd,0]^2,1./sigczw[gd,1]^2])
;           print,ncolm_global[nn],cznstr.czw_2darr[nn],median(czw[gd],/even),ngd/total(1/czw[gd]) ;,gmean(czw[gd]) ;sdss_minwad([czw[gd],czw[gd]],[1./sigczw[gd,0]^2,1./sigczw[gd,1]^2])
           ;; range
           srt = gd[sort(czw[gd])]
           frac = (findgen(ngd)+1)/ngd
           sig = interpol(czw[srt], frac, [1-gauss_pdf(1),gauss_pdf(1)]) ; +/-1sig
           cznstr.sigczw_2darr[nn,0] = cznstr.czw_2darr[nn] - sig[0]
           cznstr.sigczw_2darr[nn,1] = sig[1] - cznstr.czw_2darr[nn]
         end 
     endcase                    ; of ngd

     ;; R(z), dz, dX
     if ngd gt 0 then begin
        indx = floor((civstrct[gd].(ewtag)[0]-cmpltstr.ewlim[0])/cmpltstr.ewbinsize)
        gd2 = where(indx ge 0,ngd2)
        if gd2[0] eq -1 then $
           stop,'sdss_estczn() stop: why no EW matches to est. R(z), etc?'
        if ngd2 eq 1 then begin
           cznstr.rz_2darr[nn] = cmpltstr.rz_2darr[indx[gd2[0]]]
           cznstr.dz_2darr[nn] = cmpltstr.dz_2darr[indx[gd2[0]]]
           cznstr.dx_2darr[nn] = cmpltstr.dx_2darr[indx[gd2[0]]]
        endif else begin
           ;; Median? Min? Max? Unclear
           ;; Also might all be the same... They do seem to be the same
           cznstr.rz_2darr[nn] = median(cmpltstr.rz_2darr[indx[gd2]],/even)
           cznstr.dz_2darr[nn] = median(cmpltstr.dz_2darr[indx[gd2]],/even)
           cznstr.dx_2darr[nn] = median(cmpltstr.dx_2darr[indx[gd2]],/even)
        endelse 
     endif 

  endfor                        ; loop nn=cmpltstr.newbin


  ;; Handle extreme ends and Clean up by interpolating
  bd = where(finite(cznstr.czw_2darr) eq 0,complement=gd)
  if gd[0] ne 0 then begin
;     stop,'sdss_estczn() stop: why not flush to lower column?'
     ;; Floor this out
     cznstr.czw_2darr[0:gd[0]-1] = min(czw,imn)
     cznstr.sigczw_2darr[0:gd[0]-1,0] = sigczw[imn,0]
     cznstr.sigczw_2darr[0:gd[0]-1,1] = sigczw[imn,1]
     cznstr.rz_2darr[0:gd[0]-1] = max(cmpltstr.rz_2darr,imn)
     cznstr.dz_2darr[0:gd[0]-1] = cmpltstr.dz_2darr[imn]
     cznstr.dx_2darr[0:gd[0]-1] = cmpltstr.dx_2darr[imn]
     bd = where(finite(cznstr.czw_2darr) eq 0,complement=gd)     
  endif

  if bd[0] ne -1 then begin
     test = where(bd eq shift(bd,1)+1 and bd eq shift(bd,-1)-1)
     if test[0] ne -1 then begin
        ;; Have to max out from this to end b/c all should be largest value
        istrt = bd[test[0]]-1     
        cznstr.czw_2darr[istrt:*] = max(czw,imx)
        cznstr.sigczw_2darr[istrt:*,0] = sigczw[imx,0]
        cznstr.sigczw_2darr[istrt:*,1] = sigczw[imx,1]
        cznstr.rz_2darr[istrt:*] = max(cmpltstr.rz_2darr,imx)
        cznstr.dz_2darr[istrt:*] = cmpltstr.dz_2darr[imx]
        cznstr.dx_2darr[istrt:*] = cmpltstr.dx_2darr[imx]
        bd = where(finite(cznstr.czw_2darr) eq 0,complement=gd)
     endif

     if bd[0] ne -1 then begin
        cznstr.czw_2darr[bd] = interpol(cznstr.czw_2darr[gd], $
                                        ncolm_global[gd],$
                                        ncolm_global[bd],_extra=extra)
        cznstr.sigczw_2darr[bd,0] = interpol(cznstr.sigczw_2darr[gd,0], $
                                             ncolm_global[gd],$
                                             ncolm_global[bd],_extra=extra)
        cznstr.sigczw_2darr[bd,1] = interpol(cznstr.sigczw_2darr[gd,1], $
                                             ncolm_global[gd],$
                                             ncolm_global[bd],_extra=extra)
        cznstr.rz_2darr[bd] = interpol(cznstr.rz_2darr[gd], $
                                       ncolm_global[gd],$
                                       ncolm_global[bd],_extra=extra)
        cznstr.dz_2darr[bd] = interpol(cznstr.dz_2darr[gd], $
                                       ncolm_global[gd],$
                                       ncolm_global[bd],_extra=extra)
        cznstr.dx_2darr[bd] = interpol(cznstr.dx_2darr[gd], $
                                       ncolm_global[gd],$
                                       ncolm_global[bd],_extra=extra)
     endif 
  endif 

  cznstr.nrec_2darr = cznstr.ninput_2darr*cznstr.czw_2darr

  return, cznstr
end                             ; sdss_estczn()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_pltfxw, fxwstrct_fil, fw2_fil=fw2_fil, fxwfit_fil=fxwfit_fil, xrng=xrng,$
                 yrng=yrng, xlog=xlog, ylog=ylog, title=title, psize=psize, $
                 csize=csize, lthick=lthick, linsty=linsty, psfil=psfil,$
                 scale=scale, pltxtr=pltxtr, pclr=pclr, _extra=extra
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltfxw, fxwstrct_fil, [fxwfit_fil=, fw2_fil=, xrng=,'
     print,'                    yrng=, /xlog, /ylog, title=, psize=, scale=,'
     print,'                    csize=, lthick=, linsty=, psfil=, /pltxtr, pclr=, _extra=]'
     return
  endif
  angstrom = STRING("305B)   
 
  if not keyword_set(csize) then csize = 2.
  if not keyword_set(psize) then psize = 1.5
  if not keyword_set(lthick) then lthick = 3
  if not keyword_set(scale) then scale = 1.
 
  if keyword_set(psfil) then begin
     ;; _extra= includes /encaps
     x_psopen,psfil,/maxs,/portrait,_extra=extra
     !p.multi = [1,1,1]
     !x.margin = [7.8,1.5]      ; left and right border
     !y.margin = [3.2,0.5]      ; bottom and top border
  endif 

  if size(fxwstrct_fil,/type) eq 7 then $
     fxw = xmrdfits(fxwstrct_fil,1,/silent) $
  else fxw = fxwstrct_fil
  clr = getcolor(/load)
  if not keyword_set(pclr) then pclr = clr.black
  if not keyword_set(title) then begin
     title = ['!8W!X!Dr!N ('+angstrom+')','!8f!X(!8W!X!Dr!N)'] 
     if keyword_set(xlog) then title[0] = 'log '+title[0] 
     if keyword_set(ylog) then title[1] = 'log '+title[1]
  endif 


  ;; ;;;;;;;
  ;; Plot data
  ydat = fxw.fxw
  sigydat = fxw.sigfxw
  if keyword_set(ylog) then begin
     gd = where(fxw.fxw ne 0.,ngd,complement=bdy)
     ydat[gd] = alog10(fxw.fxw[gd])
     sigydat[gd,0] = fxw.sigfxw[gd,0]/abs(alog(10.)*fxw.fxw[gd])
     sigydat[gd,1] = fxw.sigfxw[gd,1]/abs(alog(10.)*fxw.fxw[gd])
     if bdy[0] ne -1 then begin
        ydat[bdy] = !values.f_nan
        sigydat[bdy,*] = !values.f_nan
        sub = where(bdy lt gd[ngd-1] and bdy gt gd[0])
        if sub[0] ne -1 then begin ; don't draw too manx arrows
           bdy = bdy[sub]
           sigydat[bdy,1] = alog10(fxw.sigfxw[bdy,1]) ; upper limit
        endif 
     endif      
  endif else bdy = -1                        ; /ylog

  xdat = fxw.ewcenter
  sigxdat = fxw.sigewcenter
  if keyword_set(xlog) and not keyword_set(czn) then begin
     gd = where(fxw.ewcenter ne 0.,complement=bdx)
     xdat[gd] = alog10(fxw.ewcenter[gd])
     sigxdat[gd,0] = fxw.sigewcenter[gd,0]/abs(alog(10.)*fxw.ewcenter[gd])
     sigxdat[gd,1] = fxw.sigewcenter[gd,1]/abs(alog(10.)*fxw.ewcenter[gd])
     if bdx[0] ne -1 then begin
        xdat[bdx] = !values.f_nan
        sigxdat[bdx,*] = !values.f_nan
     endif                 
     ;; For small values normal approximation breaks down (EW <
     ;; 1)
     bd = where(fxw.ewcenter[gd] lt 1.)
     if bd[0] ne -1 then begin
        xlo = fxw.ewcenter[bd] - fxw.sigewcenter[bd,0]
        xhi = fxw.ewcenter[bd] + fxw.sigewcenter[bd,1]
        gd = where(xlo ne 0.)
        sigxdat[bd[gd],0] = xdat[bd[gd]] - alog10(xlo[gd])
        sigxdat[bd,1] = alog10(xhi) - xdat[bd]
     endif 
  endif else bdx = -1                         ; /xlog


  ;; Plot
  if not keyword_set(xrng) then $
     xrng = [min(xdat-sigxdat[*,0]),max(xdat+sigxdat[*,1])]
  if not keyword_set(yrng) then $
     yrng = [min(ydat-sigydat[*,0]),max(ydat+sigydat[*,1])]

  if not keyword_set(pltxtr) then $
     plot,[0],[0],/nodata,/xsty,/ysty,color=clr.black,$
          background=clr.white,thick=lthick,charsize=csize,$
          xtitle=title[0],ytitle=title[1],xrange=xrng,yrange=yrng

  oploterror,xdat,ydat,sigxdat[*,0],sigydat[*,0],$
             /lobar, errcolor=pclr,errthick=lthick,$
             psym=3,color=pclr,thick=lthick,/nohat
  oploterror,xdat,ydat,sigxdat[*,1],sigydat[*,1],$
             /hibar, errcolor=pclr,errthick=lthick,$
             psym=3,color=pclr,thick=lthick,/nohat

  ;; Upper limit
  if bdy[0] ne -1 then begin
     plotsym,1,color=pclr
     oplot,xdat[bdy],sigydat[bdy,1],symsize=psize,thick=lthick,$
           color=pclr,psym=8
  endif                         ; bdy[0] ne -1

  if keyword_set(fw2_fil) then $
     sdss_pltfxw, fw2_fil, xlog=xlog, ylog=ylog, /pltxtr, pclr=clr.red

  if keyword_set(fxwfit_fil) then begin
     ;; Plot the fit
     if keyword_set(xlog) then $
        xfit = 10.^xrng[0] + findgen(500)*(10.^xrng[1]-10.^xrng[0])/500. $
     else xfit = xrng[0] + findgen(500)*(xrng[1]-xrng[0])/500.
     
     yfit = sdss_calcfxwfit(xfit,fitstrct_fil=fxwfit_fil) * scale
     if keyword_set(ylog) then yfit = alog10(yfit)
     if keyword_set(xlog) then xfit = alog10(xfit)

     oplot,xfit,yfit,color=clr.red,thick=lthick,linestyle=linsty
  endif 

  if keyword_set(psfil) then begin
     x_psclose
     print,'sdss_pltfxw: created ',psfil
  endif

end                             ; sdss_pltfxw


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Author: Shigenobu Hirose at JAMSTEC
;; based on original paper
;; Shimazaki and Shinomoto, Neural Computation 19, 1503-1527, 2007
;;
function sshist, data, x=x, cost=cost, nbin=nbin

  COMPILE_OPT idl2

  nbin_min = 2
  nbin_max = 200

  ntrial = nbin_max - nbin_min + 1

  nbin  = INDGEN(ntrial) + nbin_min

  delta = FLTARR(ntrial)
  cost  = FLTARR(ntrial)

  for n = 0, ntrial-1  do begin
     delta[n] = (MAX(data) - MIN(data)) / (nbin[n] - 1)

     k = HISTOGRAM(data, nbins=nbin[n])

     kmean = MEAN(k)
     kvari = MEAN((k - kmean)^2)
     cost[n] = (2. * kmean - kvari) / delta[n]^2
  endfor

  ;; Change this up to return something useful
;  n = (WHERE(cost eq MIN(cost)))[0]
  mn = MIN(cost,n)
  nbin = nbin[n]
  k = HISTOGRAM(data, nbins=nbin, locations=x)

  return, k

end



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcfxw, civstrct_fil, cmplt_fil, zlim, ewbinsize=ewbinsize, $
                       n_per_bin=n_per_bin, final=final, dz=dz, $
                       ewlim=ewlim, ewdndx=ewdndx, silent=silent, $
                       bigewbin=bigewbin, c_ech=c_ech, czn=czn, $
                       afp_flg=afp_flg, use_flg=use_flg, optimal=optimal, _extra=extra
  ;; Calculate loads of stuff associated with f(W) or f(N) in redshift
  ;; or pathlength space
  if n_params() ne 3 then begin
     print,'Syntax - sdss_calcfxw(civstrct_fil, cmplt_fil, zlim, [/final, /dz=, '
     print,'                      ewbinsize=, n_per_bin=, ewlim=, /silent, bigewbin=, '
     print,'                      ewdndx=, /c_ech, /czn, /afp_flg, /use_flg, /optimal, _extra=])'
     return,-1
  endif 

  ;; Read in data
  ;; Extract info from completeness file
  if size(cmplt_fil,/type) eq 7 then cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil
  tags = tag_names(cmpltstr)
  if keyword_set(dz) then begin
     xstr = 'z' 
     dxtag = (where(tags eq 'DZ_2DARR'))[0]
  endif else begin
     xstr = 'X'
     dxtag = (where(tags eq 'DX_2DARR'))[0]
  endelse 

  ;; Default
  if keyword_set(use_flg) then begin
     if not keyword_set(c_ech) then begin
        print,'sdss_calcfxw(): NOTE!!! for /use_flg, must set /c_ech'
        c_ech = 1
     endif
     if keyword_set(bigewbin) then begin
        print,'sdss_calcfxw(): NOTE!!! for /use_flg, bigewbin= not allowed'
        bigewbin = 0
     endif
  endif

  if size(civstrct_fil,/type) eq 8 then civstr = civstrct_fil $
  else civstr = xmrdfits(civstrct_fil,1,/silent)
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endelse 

  ;; Find what's in bins
  sub = where(civstr.(ztag)[0] ge zlim[0] and $
              civstr.(ztag)[0] lt zlim[1],nciv)
  
  if sub[0] eq -1 then begin
     print,'sdss_calcfxw(): no doublets with desired z'
     return,-1
  endif 
  civstr = civstr[sub]          ; just truncate
  ;; May need column densities so just compile
  dblt = dblt_retrieve(civstr[0].wrest[0])
  ;; _extra includes /estflg
  ncolm = sdss_calcncolmdblt(civstr,dblt,signcolm=signcolm,$
                             ncolmflg=ncolmflg,final=final,/silent,$
                             /log,_extra=extra)
  ncolmclim = sdss_getewclim(cmpltstr,/ncolm,dblt_name=dblt,/log) ; trim below

  ;; Set up the data
  if keyword_set(czn) then begin
     data = ncolm
     sigdata = signcolm
     dataflg = ncolmflg
  endif else begin
     data = civstr.(ewtag)[0]
     sigdata = civstr.(sigewtag)[0]
     dataflg = civstr.ewflg[0]
  endelse 


  ;; Handle binning
  if not keyword_set(ewlim) then begin
     if keyword_set(czn) then begin
        ewlim = [12.6,max(data+alog10(1.+alog(10.)*sigdata))] ; + 1 sigma
     endif else $
        ewlim = [0.,max(data+sigdata)] ; Ang
  endif 

  ;; ew_global defined on the left-hand side
  if keyword_set(optimal) then begin
     hist = sshist(data, x=ew_global, nbin=newbin)
     ewbinsize = ew_global[1]-ew_global[0] ; override input ewbinsize= or n_per_bin=
     dew_global = replicate(ewbinsize,newbin)
  endif 

  if not keyword_set(ewbinsize) then begin
     if keyword_set(n_per_bin) then begin
        ;; Use unequal binning
        hist = sdss_histogram(data, n_per_bin,location=ew_global,szloc=dew_global,$
                              nbin=newbin)
     endif else begin
        ewbinsize = cmpltstr.ewbinsize
     endelse                    ; default ewbisize
  endif                         ; ewbinsize=0
  
  if not keyword_set(ew_global) then begin
     ew_global = sdss_mkewarr(ewlim,ewbinsize,newbin=newbin)
     dew_global = replicate(ewbinsize,newbin)
  endif 

  if keyword_set(bigewbin) then begin
     ;; Make one final big bin
     mn = min(bigewbin-ew_global,imn,/abs) ; find closest regularly gridded point
     if imn eq 0 then begin
        newbin = 1
        ew_global = [ewlim[0] + dew_global[0]]
        dew_global = [ewlim[1] - bigewbin]
     endif else begin
        dew_global = [dew_global[0:imn-1],$
                      ew_global[newbin-1]+dew_global[newbin-1]-ew_global[imn]]
        ew_global = ew_global[0:imn]
        newbin = imn + 1
     endelse 
  endif                         ; bigewbin=
  
  if keyword_set(use_flg) then begin
     ;; All bins go from a range of the left to the very limit
     dew_global = ew_global[newbin-1] + dew_global[newbin-1] - ew_global
  endif

  ;; Prepare dN/dX (civstr is already truncated) in redshift; default
  ;; to 0.6 Ang (50% limit) or log N = 14
  if not keyword_set(ewdndx) then begin
     ewdndx = [0.6,!values.f_infinity]
     if keyword_set(czn) then ewdndx[0] = 14.
  endif
  dndx = sdss_calcdndx(civstr, cmpltstr, ewdndx, final=final, dz=dz,$
                       /silent, czn=czn, _extra=extra)
  gddndx = (where(dndx.dndx gt 0.,ngddndx))[0]
  if ngddndx ne 1 then $
     stop,'sdss_calcfxw() stop: dN/dX should be one value',dndx.zcenter[gddndx]


  ;; Set up structure
  fxw_strct = {$
              zlim:zlim, $
              dz:keyword_set(dz), $
              afp_flg:keyword_set(afp_flg), $
              use_flg:keyword_set(use_flg), $
              c_ech:keyword_set(c_ech), $ ; method
              czn:keyword_set(czn), $
              czn_cmplt:keyword_set(cmpltstr.czn), $
              ewcenter:(ew_global+0.5*dew_global), $
              sigewcenter:rebin(0.5*dew_global,newbin,2), $ ; inflate
              ewmed:fltarr(newbin,/nozero), $               ; from discrete data (cmplt weighted)
              sigewmed:fltarr(newbin,2,/nozero), $
              numtot0:fltarr(newbin,/nozero), $ ; straight-up number in bin
              signumtot0:fltarr(newbin,2,/nozero), $
              numtot:fltarr(newbin,/nozero), $ ; completeness corrected #
              signumtot:fltarr(newbin,2,/nozero), $
              dxtot:fltarr(newbin,/nozero), $ ; = dX(ewcenter)
              sigdxtot:fltarr(newbin,2,/nozero), $
              ewlim_dndx:dndx.ewlim, $ ; ewdndx
              dndx:dndx.dndx[gddndx], $
              sigdndx:reform(dndx.sigdndx[gddndx,*]), $ ; collapse
              fxw:fltarr(newbin,/nozero), $ 
              sigfxw:fltarr(newbin,2,/nozero), $
              fxwflg:intarr(newbin) $ ; upper or lower limit
              }

  ;; Set up some useful stuff
  zcent = mean(zlim)
  if (size(cmpltstr.zlim,/n_dim))[0] gt 1 then begin
     z_global = cmpltstr.zlim[*,0] 
     ;; ignore cmpltstr.zbinsize b/c has z fudge factor
     dz_global = cmpltstr.zlim[*,1] - cmpltstr.zlim[*,0] 
  endif else begin
     z_global = sdss_mkzarr(cmpltstr.zlim, cmpltstr.zbinsize)
     dz_global = replicate(cmpltstr.zbinsize, cmpltstr.nzbin)
  endelse 
  gd = (where(zcent ge z_global and zcent lt z_global+dz_global,ngd))[0] ; use gd indices below
  if keyword_set(czn) then ncolmclim = ncolmclim[gd]
  
;  ;; This was used to test results
;  if keyword_set(czn) ne keyword_set(cmpltstr.czn) then begin
;     if keyword_set(czn) then begin
;        newcmpltstr = sdss_estczn(cmpltstr, civstr, final=final, _extra=extra)
;        cmpltstr = newcmpltstr 
;        print,'sdss_calcfxw(): NOTE! Estimated C(N) from C(W)'
;     endif
;  endif 

  ;; Get necessary completeness correction values (0 to 1)
  ;; _extra includes /grid, /lsquadratic,
  ;; /quadratic, /spline
  ;; FORCE civobs_corr=0 because pathlength should already be adjusted
  ;; for part blocked by other absorbers in z and EW bin
  if keyword_set(czn) ne keyword_set(cmpltstr.czn) then begin
     ;; Then there could be a problem
     print,'sdss_calcfxw(): WARNING! Mixing EW and column density between freq. distr. and completeness'
     if keyword_set(c_ech) then begin
        ;; Up-weighting the numerator counts so need the max (100%
        ;; completeness) and can always do this; just query the right way
        if keyword_set(cmpltstr.czn) then $
           cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                                 ncolm, sigewabs=signcolm, $
                                 sigma=sigcxw_civ, dz=dz, /czw, $
                                 civobs_corr=0, final=final, dxwmax=dxwmax, _extra=extra) $
        else $
           ;; Completeness structure is compiled in EW space
           cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                                 civstr.(ewtag)[0], sigewabs=civstr.(sigewtag)[0], $
                                 sigma=sigcxw_civ, dz=dz, /czw, $
                                 civobs_corr=0, final=final, dxwmax=dxwmax, _extra=extra)

     endif else begin
        ;; Will be taking median(dX(W)) in bins for f(N)
        ;; the below should be the same as taking dxwmax * cxw_civ
        ;; (array question)
        if keyword_set(cmpltstr.czn) then $
           stop,'sdss_calcfxw() stop: WARNING! Really want to estimate C(W) from C(N) to make this symmetric.' $
;           dxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
;                                 ncolm, sigewabs=signcolm, $
;                                 sigma=sigdxw_civ, dz=dz, czw=0, $ ; <--
;                                 civobs_corr=0, final=final, dxwmax=dxwmax, _extra=extra) $
        else begin
           ;; May not be ideal to use the full distribution but have no
           ;; right to cut
           ;; _extra= includes nlim=, nbinsize=, ewlim=, _extra=
           newcmpltstr = sdss_estczn(cmpltstr, civstr, final=final, _extra=extra)
           cmpltstr = newcmpltstr 
           print,'sdss_calcfxw(): NOTE! Estimated C(N) from C(W)'
;           dxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
;                                 civstr.(ewtag)[0], sigewabs=civstr.(sigewtag)[0], $
;                                 sigma=sigdxw_civ, dz=dz, czw=0, $ ; <--
;                                 civobs_corr=0, final=final, dxwmax=dxwmax, _extra=extra) 
        endelse                 ; czn = 0
     endelse                    ; c_ech = 0 
  endif                         ; czn != cmpltstr.czn

  ;; Since may have changed cmpltstr, check again
  if keyword_set(czn) eq keyword_set(cmpltstr.czn) then begin
     ;; There is no problem b/c all consistent
     if keyword_set(c_ech) then begin
        ;; Up-weighting the numerator counts so need the max (100%
        ;; completeness) 
        cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                              data, sigewabs=sigdata, $
                              sigma=sigcxw_civ, dz=dz, /czw, $
                              civobs_corr=0, final=final, dxwmax=dxwmax, _extra=extra) 
     endif else begin
        ;; Instead of up-weighting the numerator; take the average
        ;; pathelength available in the EW bin (for given redshift bin)
        ;; Find the dX(Wbin) but at slightly larger z to make sure to get
        ;; right grid
        dxw_bin = sdss_getdxw(cmpltstr, replicate(zcent,newbin), $
                              fxw_strct.ewcenter, sigewabs=0.5*dew_global,$
                              sigma=sigdxw_bin, dz=dz, czw=0, $ ; <--
                              civobs_corr=0, dxwmax=dxwmax, _extra=extra)
     endelse 
  endif                         ; czn == cmpltstr.czn

  if keyword_set(c_ech) then begin
     dxwmax = dxwmax[gd]        ; used below
     fxw_strct.dxtot = dxwmax   ; all the same values for z bin
     fxw_strct.sigdxtot = 0.    ; must instantiate
  endif else begin
     if keyword_set(dxw_bin) then begin
        fxw_strct.dxtot = dxw_bin
        fxw_strct.sigdxtot = sigdxw_bin
     endif                      ; else going to take median in bins
  endelse 


  ;; Bin it up the hard way
  for ee=0,newbin-1 do begin
     gddef = where(data ge ew_global[ee] and $ ; deinitely in the bin
                   data lt ew_global[ee]+dew_global[ee],ngddef) 
     if keyword_set(use_flg) then begin ; maximum contribution
        gdmax = where($
                ;; stored value in the bin
                (data ge ew_global[ee] and $ 
                 data lt ew_global[ee]+dew_global[ee]) or $
                ;; or stored value is below the bin and a lower limit
                ((dataflg and sdss_getlimflg(/lower)) eq sdss_getlimflg(/lower) and $
                 data lt ew_global[ee]) or $
                ;; or stored value is above the bin and an upper limit
                ;; but with new changes to /use_flg forcing
                ;; /c_ech shouldn't matter
                ((dataflg and sdss_getlimflg(/upper)) eq sdss_getlimflg(/upper) and $
                 data ge ew_global[ee]+dew_global[ee]), $
                ngdmax) 
        if ngdmax ne ngd then fxw_strct.fxwflg[ee] = sdss_getlimflg(/upper)
        gd = gdmax
        ngd = ngdmax
     endif else begin           ; /use_flg
        gd = gddef
        ngd = ngddef
     endelse

     ;; Save all values generally
     fxw_strct.numtot0[ee] = ngd
     ;; _extra= includes sigma=, cl=; already included in sdss_calcncmplt()
     fxw_strct.signumtot0[ee,*] = $
        sdss_calcsigpoiss(fxw_strct.numtot0[ee],verbose=(not keyword_set(silent)),$
                          _extra=extra)

     if ngd eq 0 then begin
        ;; must instantiate
        fxw_strct.ewmed[ee] = fxw_strct.ewcenter[ee] ; default
        fxw_strct.sigewmed[ee,*] = fxw_strct.ewcenter[ee,*]
        fxw_strct.numtot[ee] = fxw_strct.numtot0[ee] ; yes, it's just 0
        fxw_strct.signumtot[ee,*] = fxw_strct.signumtot0[ee,*]
        fxw_strct.fxw[ee] = 0. 
        fxw_strct.sigfxw[ee,0] = 0. ; there must be a limit
        if keyword_set(czn) then begin
           if not keyword_set(c_ech) then begin
              if fxw_strct.ewcenter[ee] lt ncolmclim and $
                 not keyword_set(cmpltstr.czn) then begin ; if cmpltstr in EW space
                 ;; Assume linear COG and find value b/c better upper
                 ;; limit
                 ewtmp = fltarr(3,/nozero)
                 ewtmp[0] = ew_to_colm(dblt.wvI,10.^fxw_strct.ewcenter[ee],$
                                       /rvrs,/silent) * 1.e-3 ; Ang
                 ewtmp[1] = ew_to_colm(dblt.wvI,10.^fxw_strct.ewcenter[ee]*$
                                       (1-alog(10)*fxw_strct.sigewcenter[ee,0]),$
                                       /rvrs,/silent) * 1.e-3 ; Ang
                 ewtmp[1] = ewtmp[0] - ewtmp[1]
                 ewtmp[2] = ew_to_colm(dblt.wvI,10.^fxw_strct.ewcenter[ee]*$
                                       (1+alog(10)*fxw_strct.sigewcenter[ee,1]),$
                                       /rvrs,/silent) * 1.e-3 ; Ang
                 ewtmp[2] = ewtmp[2] - ewtmp[0]
                 ;; Do not need to pass in civstrct_fil to account for
                 ;; detection masking b/c no detections
                 fxw_strct.dxtot[ee] = $
                    sdss_getdxw(cmpltstr,zcent,ewtmp[0],dz=dz,sigewabs=ewtmp[1],$
                                civobs_corr=0,sigma=sig, _extra=extra)
                 fxw_strct.sigdxtot[ee,0] = sig[0]
                 tmp = sdss_getdxw(cmpltstr,zcent,ewtmp[0],dz=dz,sigewabs=ewtmp[2],$
                                   civobs_corr=0,sigma=sig, _extra=extra)
                 fxw_strct.sigdxtot[ee,1] = sig[1]
              endif else begin
                 ;; Take it as it comes out if > ncolmclim or
                 ;; cmpltstr.czn = 1
                 fxw_strct.fxw[ee] = 0.
                 fxw_strct.dxtot[ee] = dxwmax[0] ; default --- is this right?!!!
                 fxw_strct.sigdxtot[ee,*] = 0.
              endelse
           endif                ; c_ech = 0
           ;; log bins for column
           fxw_strct.sigfxw[ee,1] = fxw_strct.signumtot[ee,1] / $ ; limit
                                    (10.^ew_global[ee]*(10.^dew_global[ee] - 1.) $
                                     * fxw_strct.dxtot[ee]) 
        endif else $
           ;; linear bins for EW
           fxw_strct.sigfxw[ee,1] = fxw_strct.signumtot[ee,1] / $ ; limit
           (dew_global[ee] * fxw_strct.dxtot[ee])

     endif else begin
        ;; Can actually calculate
        
        if keyword_set(c_ech) then begin
           ;; Up-weight the numbers to reflect completeness
           fxw_strct.ewmed[ee] = sdss_medianw(civstr[gd].(ewtag)[0], cxw_civ[gd])
           
           ;; includes lump Poisson error
           fxw_strct.numtot[ee] = sdss_calcncmplt(cxw_civ[gd],sigcxw_civ[gd,*],$
                                                  signcmplt=signumtot)
           ;; Errors just in Number = sum(1 / C(W))
           fxw_strct.signumtot[ee,*] = signumtot
        endif else begin
           ;; Bin dX in denominator
           if not keyword_set(dxw_bin) then begin 
              stop,'sdss_calcfxw() stop: should never here, where want to take median(dxw_civ).'
              ;; mean we have EW-column mixture and will be using
              ;; individual numbers to define bin quantity
;              if ngd eq 1 then begin
;                 fxw_strct.dxtot[ee] = dxw_civ[gd[0]]
;                 fxw_strct.sigdxtot[ee,*] = sigdxw_civ[gd[0],*] ; no other choice
;              endif else begin
;                 fxw_strct.dxtot[ee] = median(dxw_civ[gd],/even)
;                 fxw_strct.sigdxtot[ee,*] = stddev(dxw_civ[gd]) ; probably larger than sigdxw_civ
;              endelse 
           endif                ; /czn

           ;; Can't weight individually (whole bin assumed to have same)
           fxw_strct.ewmed[ee] = median(civstr[gd].(ewtag)[0],/even)

           ;; Since dX(W) reflects average completeness, keep the real
           ;; numbers in numerator
           fxw_strct.numtot[ee] = fxw_strct.numtot0[ee]
           fxw_strct.signumtot[ee,*] = fxw_strct.signumtot0[ee,*]
        endelse                 ; c_ech = 0
        ;; Finish median EW 
        fxw_strct.sigewmed[ee,0] = fxw_strct.ewmed[ee] - $
                                   (fxw_strct.ewcenter[ee] - $
                                    fxw_strct.sigewcenter[ee,0]) ; handles best
        fxw_strct.sigewmed[ee,1] = (fxw_strct.ewcenter[ee] + $
                                    fxw_strct.sigewcenter[ee,1]) - $
                                   fxw_strct.ewmed[ee] 
        
        
        ;; f(X,W) = Number / (dW * dX(W))
        if keyword_set(czn) then $                      ; log
           fxw_strct.fxw[ee] = fxw_strct.numtot[ee] / $ ; log space
                               (10.^ew_global[ee]*(10.^dew_global[ee] - 1.) $
                                * fxw_strct.dxtot[ee]) $
        else fxw_strct.fxw[ee] = fxw_strct.numtot[ee] / $ ; linear
                                 (dew_global[ee] * fxw_strct.dxtot[ee])

        ;; Error combining counting error for what's in bin
        ;; (numtot0) and the completeness correction and the dX(W)
        fxw_strct.sigfxw[ee,0] = $
           fxw_strct.fxw[ee]*$
           sqrt( (fxw_strct.sigdxtot[ee,0]/fxw_strct.dxtot[ee])^2 + $
                 (fxw_strct.signumtot[ee,0]/fxw_strct.numtot[ee])^2 ) 
        fxw_strct.sigfxw[ee,1] = $
           fxw_strct.fxw[ee]*$
           sqrt( (fxw_strct.sigdxtot[ee,1]/fxw_strct.dxtot[ee])^2 + $
                 (fxw_strct.signumtot[ee,1]/fxw_strct.numtot[ee])^2 )
        ;; Actually might be better to just take the stddev of the
        ;; values in the bin (if sufficient numbers); grabs the
        ;; character better... 
     endelse 
     
  endfor                        ; loop ee=newbin


  ;; Accepted false-positive adjustment
  if keyword_set(afp_flg) then begin
     if keyword_set(czn) or keyword_set(cmpltstr.czn) then $
        stop,'sdss_calcfxw() stop: accepted false-positive adjustment does not work for columns'

     ;; _extra= includes stuff for sdss_getdxzw
     afpdndx = sdss_getafpdndx(fxw_strct.zlim, fxw_strct.ewlim_dndx, $
                               dz=fxw_strct.dz, _extra=extra)
     nwfxw_strct = fxw_strct

     ;; Add information to structure
     nwfxw_strct.numtot0 = fxw_strct.numtot0 - afpdndx.numtot0
     nwfxw_strct.signumtot0 = sqrt(fxw_strct.signumtot0^2 + afpdndx.signumtot0^2)
     nwfxw_strct.numtot = fxw_strct.numtot - afpdndx.numtot
     nwfxw_strct.signumtot = sqrt(fxw_strct.signumtot^2 + afpdndx.signumtot^2)
     ;; dxtot, sigdxtot do *not* change
     ;; scale = (1 - (dN/dX_afp)/(dN/dX_0))
     ;; f(W) = scale * f(W)_0
     ;; var(f(W)) = var(dN/dX_afp)/(dN/dX_0)^2 * f(W)_0^2 + 
     ;;             var(dN/dX_0)*( (dN/dX_afp)/(dN/dX_0)^2 * f(W)_0 )^2
     ;;             + var(f(W)_0) * scale^2 +
     ;;             2*sig(dN/dX_0)*sig(f(W)_0)*((dN/dX_afp)/(dN/dX_0)^2*scale*f(W)_0)
     ;; Last term is co-variance, dN/dX_afp acts opposite of dN/dX_0
     ;; and f(W)_0 w.r.t. how f(W) goes up or down.
     scale = 1. - afpdndx.dndx[0]/fxw_strct.dndx
     nwfxw_strct.fxw = scale * fxw_strct.fxw
     bd = where(fxw_strct.sigfxw[*,0] eq 0.)
     nwfxw_strct.sigfxw[*,0] = $
        (afpdndx.sigdndx[0,1]/fxw_strct.dndx * fxw_strct.fxw)^2 + $
        (fxw_strct.sigdndx[0]*afpdndx.dndx[0]/fxw_strct.dndx^2 * fxw_strct.fxw)^2 + $
        (fxw_strct.sigfxw[*,0] * scale)^2 + $
        2*fxw_strct.sigdndx[0]*fxw_strct.sigfxw[*,0] * $
        (afpdndx.dndx[0]/fxw_strct.dndx^2 * scale * fxw_strct.fxw)
     nwfxw_strct.sigfxw[*,1] = $
        (afpdndx.sigdndx[0,0]/fxw_strct.dndx * fxw_strct.fxw)^2 + $
        (fxw_strct.sigdndx[1]*afpdndx.dndx[0]/fxw_strct.dndx^2 * fxw_strct.fxw)^2 + $
        (fxw_strct.sigfxw[*,1] * scale)^2 + $
        2*fxw_strct.sigdndx[1]*fxw_strct.sigfxw[*,1] * $
        (afpdndx.dndx[0]/fxw_strct.dndx^2 * scale * fxw_strct.fxw)
     nwfxw_strct.sigfxw = sqrt(nwfxw_strct.sigfxw)
     if bd[0] ne -1 then $
        nwfxw_strct.sigfxw[*,0] =  0. ; to indicate limit; upper limit inflated

     ;; Re-calculate dN/dX (self-consistent), _extra= includes /noBAL
     dndx = sdss_calcdndx(civstr, cmpltstr, fxw_strct.ewlim_dndx, final=final, dz=dz,$
                          /silent, afp_flg=afp_flg, _extra=extra)
     nwfxw_strct.dndx = dndx.dndx[0]
     nwfxw_strct.sigdndx = reform(dndx.sigdndx[0,*])

     fxw_strct = nwfxw_strct
  endif                         ; /afp_flg

  if not keyword_set(silent) then $
     sdss_prntfxw, fxw_strct, _extra=extra ; /skip_null, out_fil=

  return, fxw_strct
end                             ; sdss_calcfxw()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_mkfxwstacks, civstrct_fil, outroot, zlim, ewbinsize=ewbinsize, $
                      n_per_bin=n_per_bin, final=final, $
                      ewlim=ewlim, bigewbin=bigewbin, $
                      afp_flg=afp_flg, _extra=extra
  ;; Calculate loads of stuff associated with dN/dX (or dN/dz)
  if n_params() ne 3 then begin
     print,'Syntax - sdss_mkfxwstacks, civstrct_fil, outroot, zlim, '
     print,'                           [/final, ewbinsize=, n_per_bin=, ewlim=,'
     print,'                          bigewbin=, /afp_flg, _extra=]'
     return
  endif 

  ;; Read in data
  if size(civstrct_fil,/type) eq 8 then civstr = civstrct_fil $
  else civstr = xmrdfits(civstrct_fil,1,/silent)
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endelse 

  ;; Find what's in redshift range
  sub = where(civstr.(ztag)[0] ge zlim[0] and $
              civstr.(ztag)[0] lt zlim[1],nciv)
  
  if sub[0] eq -1 then begin
     print,'sdss_mkfxwstacks: no doublets with desired z'
     return
  endif 
  civstr = civstr[sub]          ; just truncate

  ;; Handle EW binning
  if not keyword_set(ewlim) then begin
        ewlim = [0.,max(civstr.(ewtag)[0]+civstr.(sigewtag)[0])] ; Ang
  endif 
  if not keyword_set(ewbinsize) then begin
     if keyword_set(n_per_bin) then begin
        ;; Use unequal binning
        hist = sdss_histogram(civstr.(ewtag)[0], n_per_bin,location=ew_global,szloc=dew_global,$
                              nbin=newbin)
     endif else begin
        ewbinsize = cmpltstr.ewbinsize
        
        ew_global = sdss_mkewarr(ewlim,ewbinsize,newbin=newbin)
        dew_global = replicate(ewbinsize,newbin)
     endelse                    ; default ewbisize
  endif                         ; ewbinsize=0

  if keyword_set(bigewbin) then begin
     ;; Make one final big bin
     mn = min(bigewbin-ew_global,imn,/abs) ; find closest regularly gridded point
     dew_global = [dew_global[0:imn-1],$
                   ew_global[newbin-1]+dew_global[newbin-1]-ew_global[imn]]
     ew_global = ew_global[0:imn]
     newbin = imn + 1
  endif                         ; bigewbin=
  
  ew_global[0] = ew_global[0] > 0. ; sanity

  ;; Bin it up the hard way
  for ee=0,newbin-1 do begin
     gd = where(civstr.(ewtag)[0] ge ew_global[ee] and $
                civstr.(ewtag)[0] lt ew_global[ee]+dew_global[ee],ngd) 
     if ngd eq 0 then continue 

     ewstr = string(ew_global[ee],ew_global[ee]+dew_global[ee],$
                    format='(f4.2,"w",f4.2)')
     if keyword_set(afp_flg) then begin
        ;; Randomly exclude the appropriate amount of spectra, I
        ;; suppose 
        stop,'sdss_mkfxwstacks stop: /afp_flg option not functioning'
     endif 

     ;; Create stack
     ;; _extra= includes /debug, /clobber, gwave=, wvmnx=, wvnrm=,
     ;; wvmsk=, cmplt_fil=, /nowgt, /civobs_corr, [everyn=, sset=,
     ;; maxrej=, lower=, upper=, nord=, /groupbadpix, /sticky,
     ;; bsplmask=, /silent]
     ofil = outroot + '_' + ewstr + '.fit'
     sdss_stackciv, civstr[gd], ofil, final=final, _extra=extra

  endfor                        ; loop ee=newbin

  print,'sdss_mkfxwstacks: All done!'
end                             ; sdss_mkfxwstacks



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcfxwchisq, civstrct_fil, cmplt_fil, zlim=zlim, $
                            final=final, _extra=extra
  ;; Calc chi^2 from fit and binned distribution
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcfxwchisq(civstrct_fil, cmplt_fil, [_extra=])'
     return,-1
  endif 

  if size(civstrct_fil,/type) eq 8 then $
     civstr = civstrct_fil $
  else civstr = xmrdfits(civstrct_fil,1,/silent) 
  if size(cmplt_fil,/type) eq 8 then $
     cmpltstr = cmplt_fil $
  else cmpltstr = xmrdfits(cmplt_fil,1,/silent) 

  ;; Chi^2 for just the data used in the fit
  ;; _extra= includes /dz, /final, /silent, n_per_bin=
  if not keyword_set(zlim) then zlim = reform(cmpltstr.zlim)
  fxw_strct = sdss_calcfxw(civstr, cmpltstr, zlim, _extra=extra)
  ;; _extra= includes option=, fitstrct_fil=, coeff=, alpha=,
  ;; datanorm=
  if keyword_set(fxw_strct.czn) then $
     yfit = sdss_calcfxwfit(10.^fxw_strct.ewcenter,ndim=ndim,_extra=extra) $
  else yfit = sdss_calcfxwfit(fxw_strct.ewcenter,ndim=ndim,_extra=extra)
  err = fxw_strct.sigfxw[*,0]         ; lower errors for points above fit
  gd = where(fxw_strct.fxw lt yfit)   ; upper errors for points below fit
  if gd[0] ne -1 then err[gd] = fxw_strct.sigfxw[gd,1] ; upper errors
  gd = where(fxw_strct.fxw ne 0.,ndof)                 ; only good bins
  chi_sqr = total((fxw_strct.fxw[gd]-yfit[gd])^2/err[gd]^2)
  chi_sqr = [chi_sqr, ndof - ndim, 0.] ; -1 extra in dof?
  chi_sqr[2] = chisqr_pdf(chi_sqr[0], chi_sqr[1])

  return, chi_sqr
end                             ; sdss_calcfxwchisq()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_prntomega, omega_fil, skip_null=skip_null, out_fil=out_fil
  if n_params() ne 1 then begin
     print,'Syntax - sdss_prntomega, omega_fil, [/skip_null,out_fil=]'
     return
  endif 
  scale = 1.e8                  ; scale up Omega's for printing
  linf = 'Omega scaled up by '+string(scale,format='(e8.1)')
  if size(omega_fil,/type) eq 7 then $
     omega = xmrdfits(omega_fil,1,/silent) $
  else omega = omega_fil
  nstrct = (size(omega,/dim))[0] > 1

  if keyword_set(out_fil) then $

     openw,1,out_fil

  for ss=0,nstrct-1 do begin
     lin0 = strtrim(string(omega[ss].ewlim[0],format='(f7.2)'),2)+' Ang <= EW < '+$
            strtrim(string(omega[ss].ewlim[1],format='(f7.2)'),2)+' Ang; ' + $
            strtrim(string(omega[ss].ncolmlim[0],format='(f7.2)'),2)+' <= log N < '+$
            strtrim(string(omega[ss].ncolmlim[1],format='(f7.2)'),2)
     lin = string('zlo','zhi','Num','Errlo','Errhi',$
                  'dX(W)','Errlo','Errhi',$
                  'log Nmin','log Ntot','Errlo','Errhi','Omega','Errlo','Errhi',$
                  format='(2(a6,1x),3(a5,1x),3(a8,1x),4(a8,1x),3(a8,1x))')
     if keyword_set(out_fil) then begin
        printf,1,''
        printf,1,lin0
        printf,1,lin
     endif else begin
        print,''
        print,lin0
        print,lin
     endelse

     nzbin = (size(omega[ss].zcenter,/dim))[0] > 1
     for zz=0,nzbin-1 do begin
        if keyword_set(skip_null) and omega[ss].omega[zz] eq 0. then $
           continue 
        lin = $
           string(omega[ss].zcenter[zz]-omega[ss].sigzcenter[zz,0],$
                  omega[ss].zcenter[zz]+omega[ss].sigzcenter[zz,1],$
                  omega[ss].numtot[zz],$
                  omega[ss].signumtot[zz,0],omega[ss].signumtot[zz,1],$
                  omega[ss].dxtot[zz],$
                  omega[ss].sigdxtot[zz,0],omega[ss].sigdxtot[zz,1],$
                  omega[ss].ncolmmin[zz],omega[ss].ncolmtot[zz],$
                  omega[ss].signcolmtot[zz,0],omega[ss].signcolmtot[zz,1],$
                  omega[ss].omega[zz]*scale,omega[ss].sigomega[zz,0]*scale,$
                  omega[ss].sigomega[zz,1]*scale,$
                  format='(2(f6.2,1x),i5,1x,2(f5.1,1x),3(f8.2,1x),4(f8.3,1x),3(f8.4,1x))')
        if keyword_set(out_fil) then printf,1,lin $
        else print,lin
     endfor                     ; loop zz=nzbin

     if keyword_set(out_fil) then printf,1,linf $
     else  print,linf
  endfor                        ; loop ss=nstrct


  if keyword_set(out_fil) then begin
     close,1
     print,'sdss_prntomega: created ',out_fil
  endif

end                             ; sdss_prntomega



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_calcomega, civstrct_fil, cmplt_fil, c_ech=c_ech, $
                         final=final, ewlim=ewlim, nlim=nlim, silent=silent,$
                         afp_flg=afp_flg, _extra=extra
  ;; 13 Nov 2012
  print,'sdss_calcomega(): NOTE! sdss_calcfxw() has changed and this function has not matched; likely need to call sdss_estczn() too.'

  ;; Calculate the mass density
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcomega( civstrct_fil, cmplt_fil, [/c_ech, '
     print,'                         /final, ewlim=, nlim=, /silent, /afp_flg, _extra=])'
     return, -1
  endif 

  if not keyword_set(ewlim) then ewlim = [0.,!values.f_infinity]
  if not keyword_set(nlim) then nlim = [10.^13,!values.f_infinity] ;must be linear
  if nlim[0] lt 100 or nlim[1] lt 100 then $
     stop,'sdss_calcomega() stop: column limits must be linear',nlim

  ;; Read in data
  if size(cmplt_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil

  if size(civstrct_fil,/type) eq 7 then $
     civstr = xmrdfits(civstrct_fil,1,/silent) $
  else civstr = civstrct_fil
  nciv = (size(civstr,/dim))[0]
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
     ncolmtag = (where(tags eq 'NCOLM_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
     ncolmtag = (where(tags eq 'NCOLM_ORIG'))[0]
  endelse 

  ;; figure out the doublet
  dblt = dblt_retrieve(civstr[0].wrest[0])
  ;; Use information from both lines to set doublet column density 
  ;; This is slow b/c loops over everyone individual system
  ;; _extra= includes /estflg
  ncolm = sdss_calcncolmdblt(civstr,dblt,signcolm=signcolm,final=final,$
                             ncolmflg=ncolmflg,/silent,log=0,$ ; must be linear
                             _extra=extra)


  ;; Constants
  cosmology = sdss_setcosmology(_extra=extra)
  rhocrit = 1.89e-29 * (cosmology[0]/100.)^2 ; g/cm^3
  h_invs = cosmology[0] / (1.e6 * 3.09e13)   ; s^-1
  sol = 2.99792458e10                        ; cm/s
  amu = 1.66053886e-24                       ; g
  case dblt.ion of
     'SiIV': mciv = 28.0855 * amu 
     'CIV': mciv = 12.0107 * amu ; g 
     'MgII': mciv = 24.305 * amu
     else: stop,'sdss_calcomega() stop: unknown doublet ',dblt.ion
  endcase
  om_coeff = h_invs * mciv / (sol * rhocrit)
  
  ;; ;;;;;;;
  ;; Trim
  gd = where(civstr.(ewtag)[0] ge ewlim[0] and $
             civstr.(ewtag)[0] lt ewlim[1] and $
             ncolm ge nlim[0] and ncolm lt nlim[1],nciv)
  if nciv eq 0 then $
     stop,'sdss_calcomega() stop: no doublets in EW and column limits'
  civstr = civstr[gd]
  ncolm = ncolm[gd]
  signcolm = signcolm[gd]
  ncolmflg = ncolmflg[gd]

  ;; Prepare dN/dX structure
  if keyword_set(cmpltstr.czn) then $
     dndx = sdss_calcdndx(civstr, cmpltstr, nlim, final=final, /silent) $
  else $
     dndx = sdss_calcdndx(civstr, cmpltstr, ewlim, final=final, /silent)

  ;; ;;;;;;;
  ;; Set output
  omega_strct = {$
                ewlim:ewlim, $
                afp_flg:keyword_set(afp_flg), $
                ncolmlim:alog10(nlim), $
                c_ech:keyword_set(c_ech), $ ; method
                czn_cmplt:keyword_set(cmpltstr.czn), $
                zcenter:fltarr(cmpltstr.nzbin,/nozero), $
                sigzcenter:fltarr(cmpltstr.nzbin,2,/nozero), $
                numtot:fltarr(cmpltstr.nzbin,/nozero), $
                signumtot:fltarr(cmpltstr.nzbin,2,/nozero), $
                ncolmmin:fltarr(cmpltstr.nzbin,/nozero), $
                ncolmtot:fltarr(cmpltstr.nzbin,/nozero), $
                signcolmtot:fltarr(cmpltstr.nzbin,2,/nozero), $
                dxtot:fltarr(cmpltstr.nzbin,/nozero), $
                sigdxtot:fltarr(cmpltstr.nzbin,2,/nozero), $
                dndx:dndx.dndx, $ ; ewlim of Omega
                sigdndx:dndx.sigdndx, $
                omega:fltarr(cmpltstr.nzbin,/nozero), $
                sigomega:fltarr(cmpltstr.nzbin,2,/nozero) $
                }
  if cmpltstr.nzbin gt 1 then begin
     omega_strct.zcenter = (0.5*(cmpltstr.zlim[*,0] + cmpltstr.zlim[*,1])) 
     omega_strct.sigzcenter[*,0] = omega_strct.zcenter - cmpltstr.zlim[*,0]
     omega_strct.sigzcenter[*,1] = cmpltstr.zlim[*,1] - omega_strct.zcenter
  endif else begin
     omega_strct.zcenter = mean(cmpltstr.zlim) 
     omega_strct.sigzcenter[*,0] = omega_strct.zcenter - cmpltstr.zlim[0]
     omega_strct.sigzcenter[*,1] = cmpltstr.zlim[1] - omega_strct.zcenter
  endelse 


  ;; ;;;;;;;
  ;; Un-blocked pathlength selected on EW still (because that's
  ;; best, for now) 
  ;; _extra includes /grid, /ewmax, dxwmax=
  ;; FORCE civobs_corr=0 because pathlength should already be adjusted
  ;; for part blocked by other absorbers in z and EW bin
  if keyword_set(cmpltstr.czn) then $
     cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                           alog10(ncolm), sigewobs=signcolm/(alog(10)*ncolm), $
                           sigma=sigcxw_civ, dz=0, /czw, $
                           civobs_corr=0, dxwmax=dxwmax, /silent,_extra=extra) $
  else $
     cxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                           civstr.(ewtag)[0], sigewobs=civstr.(sigewtag)[0], $
                           sigma=sigcxw_civ, dz=0, /czw, $
                           civobs_corr=0, dxwmax=dxwmax, /silent,_extra=extra)

  if keyword_set(c_ech) then begin
     omega_strct.dxtot = dxwmax ; max value (100% complete)
     omega_strct.sigdxtot = 0.
  endif else begin
     ;; Will be taking median(dX(W)) in bins 
     if keyword_set(cmpltstr.czn) then $
        dxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                              alog10(ncolm), sigewabs=signcolm/(alog(10)*ncolm), $
                              sigma=sigdxw_civ, dz=0, czw=0, $ ; <--
                              civobs_corr=0, final=final, _extra=extra) $
     else $
        dxw_civ = sdss_getdxw(cmpltstr, civstr.(ztag)[0], $
                              civstr.(ewtag)[0], sigewabs=civstr.(sigewtag)[0], $
                              sigma=sigdxw_civ, dz=0, czw=0, $ ; <--
                              civobs_corr=0, final=final, _extra=extra)
  endelse 

  
  ;; loop over redshift bins
  for zz=0,cmpltstr.nzbin-1 do begin
     zlim = omega_strct.zcenter[zz]+$
            [-omega_strct.sigzcenter[zz,0], omega_strct.sigzcenter[zz,1]]

     gd = where(civstr.(ztag)[0] ge zlim[0] and civstr.(ztag)[0] lt zlim[1],ngd)

     ;; Save all values generally
     omega_strct.numtot[zz] = ngd
     ;; _extra includes sigma=, cl=,
     omega_strct.signumtot[zz,*] = $
        sdss_calcsigpoiss(float(ngd),verbose=(not keyword_set(silent)),$
                          _extra=extra)
     if ngd eq 0 then begin
        ;; must instantiate
        if not keyword_set(c_ech) then begin
           omega_strct.dxtot[zz] = dxwmax[zz] ; max value (100% complete)
           omega_strct.sigdxtot[zz,*] = 0.
        endif 
        omega_strct.omega[zz] = 0. 
        omega_strct.sigomega[zz,0] = 0. ; there must be a limit; maybe put in 50% limit value
        omega_strct.ncolmtot[zz] = 0.
        omega_strct.signcolmtot[zz,*] = 0.
        omega_strct.ncolmmin[zz] = !values.f_nan
     endif else begin
        omega_strct.ncolmmin[zz] = alog10(min(ncolm[gd]))

        if keyword_set(c_ech) then begin
           ;; N_CIV = total( N_i / C(W_i) )
           ;; sigN_CIV^2 = sigN_obs^2 + total( (sigN_i/C(W_i))^2 + 
           ;;              (N_i*sigC(W_i) / C(W_i)^2)^2 )
           omega_strct.ncolmtot[zz] = $
              sdss_calcncolm(ncolm[gd], cxw_civ[gd], signcolm[gd,*], $
                             sigcxw_civ[gd,*],signcolm=signcolmtot)
           omega_strct.signcolmtot[zz,*] = signcolmtot
        endif else begin
           ;; N_CIV = total(N_i) 
           ;; sigN_CIV^2 = total(sigN_CIV^2) + sigN_obs^2
           omega_strct.ncolmtot[zz] = $
              sdss_calcncolm(ncolm[gd], replicate(1.,ngd), signcolm[gd,*], $
                             signcolm=signcolmtot)              
           omega_strct.signcolmtot[zz,*] = signcolmtot

           if ngd eq 1 then begin
              omega_strct.dxtot[zz] = dxw_civ[gd[0]]
              omega_strct.sigdxtot[zz,*] = sigdxw_civ[gd[0],*] ; no other choice
           endif else begin
              omega_strct.dxtot[zz] = median(dxw_civ[gd],/even)
              omega_strct.sigdxtot[zz,*] = stddev(dxw_civ[gd]) ; probably larger than sigdxw_civ
           endelse 
           
        endelse 
        
        ;; Omega = om_ceff * N_CIV / dX(z)
        ;; sigOmega = om_ceff / dX(z) * sigN_CIV
        omega_strct.omega[zz] = om_coeff * omega_strct.ncolmtot[zz] / $
                                omega_strct.dxtot[zz]
        omega_strct.sigomega[zz,0] = $
           omega_strct.omega[zz] * $
           sqrt( (omega_strct.signcolmtot[zz,0]/omega_strct.ncolmtot[zz])^2 + $
                 (omega_strct.sigdxtot[zz,0]/omega_strct.dxtot[zz])^2 )
        omega_strct.sigomega[zz,1] = $
           omega_strct.omega[zz] * $
           sqrt( (omega_strct.signcolmtot[zz,1]/omega_strct.ncolmtot[zz])^2 + $
                 (omega_strct.sigdxtot[zz,1]/omega_strct.dxtot[zz])^2 )

        ;; Log for storage
        omega_strct.signcolmtot[zz,*] = $
           omega_strct.signcolmtot[zz,*] / (alog(10)*omega_strct.ncolmtot[zz])
        omega_strct.ncolmtot[zz] = alog10(omega_strct.ncolmtot[zz])
     endelse 


     if keyword_set(afp_flg) then begin
        if keyword_set(cmpltstr.czn) then $
           stop,'sdss_calcomega() stop: accepted false-positive adjustment does not work for columns'

        ;; Correct a given omega structure for accepted false positive
        ;; dN/dX_afp (and this is just a scaling, no approximation of
        ;; what column density these contribute... because unknown)
        ;; dN/dX = dN/dX_0 - dN/dX_afp
        ;; sig(dN/dX)^2 = sqrt( sig(dN/dX_0)^2 - sig(dN/dX_afp)^2 )
        ;; dN/dX structure is organized so that there are multiple
        ;; redshift bins which must be handled indepedently

        ;; For now, going to read in the file every loop
        ;; afndndx will only be for one bin
        ;; _extra includes stuff for sdss_calcdxzw()
        afpdndx = sdss_getafpdndx(zlim, omega_strct.ewlim,_extra=extra)
        
        nwomega_strct = omega_strct
        nwomega_strct.numtot[zz] = omega_strct.numtot[zz] - afpdndx.numtot
        nwomega_strct.signumtot[zz,*] = sqrt(omega_strct.signumtot[zz,*]^2 + $
                                             afpdndx.signumtot[0,*]^2)
        ;; dxtot, sigdxtot do *not* change
        ;; scale = (1 - (dN/dX_afp)/(dN/dX_0))
        ;; Omega = scale * Omega_0
        ;; var(Omega) = var(dN/dX_afp)/(dN/dX_0)^2 * Omega_0^2 + 
        ;;             var(dN/dX_0)*( (dN/dX_afp)/(dN/dX_0)^2 * Omega_0 )^2
        ;;             + var(Omega_0) * scale^2 +
        ;;             2*sig(dN/dX_0)*sig(Omega_0)*((dN/dX_afp)/(dN/dX_0)^2*scale*Omega_0)
        ;; Last term is co-variance, dN/dX_afp acts opposite of dN/dX_0
        ;; and Omega_0 w.r.t. how Omega goes up or down.
        scale = 1. - afpdndx.dndx[0]/omega_strct.dndx[zz]
        nwomega_strct.omega[zz] = scale*omega_strct.omega[zz]
        if omega_strct.omega[zz] ne 0. then begin ; don't much with other errors
           nwomega_strct.sigomega[zz,0] = $
              (afpdndx.sigdndx[0,1]/omega_strct.dndx[zz] * omega_strct.omega[zz])^2 + $
              (omega_strct.sigdndx[zz,0]*afpdndx.dndx[0]/omega_strct.dndx[zz]^2 * omega_strct.omega[zz])^2 + $
              (omega_strct.sigomega[zz,0] * scale)^2 + $
              2*omega_strct.sigdndx[zz,0]*omega_strct.sigomega[zz,0] * $
              (afpdndx.dndx[0]/omega_strct.dndx[zz]^2 * scale * omega_strct.omega[zz])
           nwomega_strct.sigomega[zz,1] = $
              (afpdndx.sigdndx[0,0]/omega_strct.dndx[zz] * omega_strct.omega[zz])^2 + $
              (omega_strct.sigdndx[zz,1]*afpdndx.dndx[0]/omega_strct.dndx[zz]^2 * omega_strct.omega[zz])^2 + $
              (omega_strct.sigomega[zz,1] * scale)^2 + $
              2*omega_strct.sigdndx[zz,1]*omega_strct.sigomega[zz,1] * $
              (afpdndx.dndx[0]/omega_strct.dndx[zz]^2 * scale * omega_strct.omega[zz])
        endif 

        nwomega_strct.sigomega[zz,*] = sqrt(nwomega_strct.sigomega[zz,*])
        omega_strct = nwomega_strct
     endif                      ; /afp_flg

  endfor                        ; loop zz=cmpltstr.nzbin

  ;; Accepted false-positive adjustment (just like f(W)
  if keyword_set(afp_flg) then begin
     ;; Re-calculate dN/dX (self-consistent), _extra= includes ...
     dndx = sdss_calcdndx(civstr, cmpltstr, omega_strct.ewlim, final=final,$
                          /silent, afp_flg=afp_flg, _extra=extra)
     nwomega_strct.dndx = dndx.dndx
     nwomega_strct.sigdndx = dndx.sigdndx
     omega_strct = nwomega_strct
  endif                         ; /afp_flg

  if not keyword_set(silent) then $
     sdss_prntomega, omega_strct, _extra=extra ; includes /skip_null, out_fil=

  return, omega_strct
end                             ; sdss_calcomega()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_combineciv, civstrct_fil, count=count, index=index, $
                          complement=complement, icomplement=icomplement, $
                          ncomplement=ncomplement, keep=keep, ntot=ntot, $
                          dvsys=dvsys, weight=weight, debug=debug, _extra=extra
  ;; Remove duplicate candidates based on qso_name, zabs_orig, and
  ;; cflg as in sdss_chkciv
  if n_params() ne 1 then begin
     print,'Syntax - sdss_combineciv( civstrct_fil, [count=, index=,'
     print,'                          complement=, ncomplement=, '
     print,'                          icomplement=, keep=, dvsys=, weight=, /debug, /ntot, _extra=])'
  endif 
  c = 299792.458d               ; km/s
  
  if not keyword_set(dvsys) then dvsys = sdss_getspecpixscale() ; 69 km/s 
  if size(civstrct_fil,/type) eq 8 then civstrct = civstrct_fil $
  else civstrct = xmrdfits(civstrct_fil,1,/silent)
  if keyword_set(weight) then begin
     if size(weight,/type) eq 7 then wgt_opt = weight $
     else wgt_opt = 'ew'        ; default
  endif else wgt_opt = 'mean'   ; else

  civstrct = sdss_srtcivstrct(civstrct)
  if keyword_set(debug) then civstrct0 = civstrct
  nciv = (size(civstrct,/dim))[0] > 1 ; defeats singularity prob
  dvnext_all = fltarr(nciv,/nozero)

  unqqso = uniq(civstrct.qso_name) ; defined by the last per block
  nqso = (size(unqqso,/dim))[0] > 1
  if nqso eq 1 then nciv_per_qso = [unqqso[0]+1] $
  else $
     nciv_per_qso = [unqqso[0] + 1,(shift(unqqso,-1)-unqqso)[0:nqso-2]]

  keep = [-1]                      ;  must reset
  for qq=0L,nqso-1 do begin
     if qq gt 0 then rng = unqqso[qq-1] + lindgen(nciv_per_qso[qq]) + 1 $
     else rng = lindgen(nciv_per_qso[qq])

     if nciv_per_qso[qq] eq 1 then begin
        dvnext_all[rng] = !values.f_nan
        continue                ; nothing to combine
     endif 

     ;; Since sorted, only need to look one way; up the stack
     ;; will all be positive 
     dvnext = c*(shift(civstrct[rng].zabs_orig[0],-1)-civstrct[rng].zabs_orig[0])/$
              (1 + civstrct[rng].zabs_orig[0])
     dvnext = [dvnext[0:nciv_per_qso[qq]-2],!values.f_nan]
     dvnext_all[rng] = dvnext ; save for review


     ;; Find consecutive triggered 
     dup = where(dvnext lt dvsys,ndup)
     if ndup ne 0 then begin
        dup = rng[dup]          ; both consecutive
        if ndup eq 1 then begin
           dup = dup[0]
           ;; Easy case
           if keep[0] ne -1 then begin
              keep = [keep,dup]
              icomplement = [icomplement,dup+1]
           endif else begin
              keep = [dup]
              icomplement = dup+1
           endelse 

           ;; Combine information for sdss_ewciv
           idx = [dup,dup+1]
           case wgt_opt of
              'ew': begin
                 wgtI = civstrct[idx].ew_orig[0]
                 wgtII = civstrct[idx].ew_orig[1]
              end
              'ncolm': begin
                 wgtI = 10.^civstrct[idx].ncolm_orig[0]
                 wgtII = 10.^civstrct[idx].ncolm_orig[1]
              end
              else: begin       ; equal (mean)
                 wgtI = idx*0 + 1.
                 wgtII = wgtI
              endelse 
           endcase
;           civstrct[idx].zabs_orig[0] = total(civstrct[idx].zabs_orig[0]*wgtI)/total(wgtI)
;           civstrct[idx].zabs_orig[1] = total(civstrct[idx].zabs_orig[1]*wgtII)/total(wgtII)
;           civstrct[idx].rating[0] = max(civstrct[idx].rating[0])
;           if keyword_set(ntot) then begin
;              civstrct[idx].ncolm_orig[0] = alog10(total(10.^civstrct[idx].ncolm_orig[0]))
;              civstrct[idx].ncolm_orig[1] = alog10(total(10.^civstrct[idx].ncolm_orig[1]))
;           endif

        endif else begin
           ;; More complicated because have multiple systems
           gapstrt = dup[where(dup ne shift(dup,1)+1,ngap)]
           gapstop = dup[where(dup ne shift(dup,-1)-1)]
           for dd=0L,ngap-1 do begin
              idx = gapstrt[dd] + 1 + lindgen(gapstop[dd]-gapstrt[dd]+1)
              if keep[0] ne -1 then begin
                 keep = [keep,gapstrt[dd]]
                 icomplement = [icomplement,idx]
              endif else begin
                 keep = [gapstrt[dd]] ; lets 0 be option
                 icomplement = idx
              endelse 
              
              ;; Combine information for sdss_ewciv
              idx = [gapstrt[dd],idx] ; accepted one and partners
              case wgt_opt of
                 'ew': begin
                    wgtI = civstrct[idx].ew_orig[0]
                    wgtII = civstrct[idx].ew_orig[1]
                 end
                 'ncolm': begin
                    wgtI = 10.^civstrct[idx].ncolm_orig[0]
                    wgtII = 10.^civstrct[idx].ncolm_orig[1]
                 end
                 else: begin                            ; equal (mean)
                    wgtI = idx*0 + 1.
                    wgtII = wgtI
                 endelse 
              endcase

           endfor               ; loop dd=ndup

        endelse ; mutliple systems

        civstrct[idx].zabs_orig[0] = total(civstrct[idx].zabs_orig[0]*wgtI)/total(wgtI)
        civstrct[idx].zabs_orig[1] = total(civstrct[idx].zabs_orig[1]*wgtII)/total(wgtII)
        civstrct[idx].rating[0] = max(civstrct[idx].rating[0])
        if keyword_set(ntot) then begin
           civstrct[idx].ncolm_orig[0] = alog10(total(10.^civstrct[idx].ncolm_orig[0]))
           civstrct[idx].ncolm_orig[1] = alog10(total(10.^civstrct[idx].ncolm_orig[1]))
        endif
        
        civstrct[idx].wvlim_orig[0,0] = min(civstrct[idx].wvlim_orig[0,0]) ; probably useless
        civstrct[idx].wvlim_orig[0,1] = max(civstrct[idx].wvlim_orig[0,1])
        civstrct[idx].wvlim_orig[1,0] = min(civstrct[idx].wvlim_orig[1,0])
        civstrct[idx].wvlim_orig[1,1] = max(civstrct[idx].wvlim_orig[1,1])

        bd = where(civstrct[idx].wvlim_orig[0,0] gt civstrct[idx].wrest[0]*(1+civstrct[idx].zabs_orig[0]) or $
                   civstrct[idx].wvlim_orig[0,1] lt civstrct[idx].wrest[0]*(1+civstrct[idx].zabs_orig[0]) or $
                   civstrct[idx].wvlim_orig[1,0] gt civstrct[idx].wrest[1]*(1+civstrct[idx].zabs_orig[1]) or $
                   civstrct[idx].wvlim_orig[1,1] lt civstrct[idx].wrest[1]*(1+civstrct[idx].zabs_orig[1]))
        if bd[0] ne -1 then begin
           printcol,civstrct[idx].zabs_orig[0],civstrct[idx].zabs_orig[1],$
                    civstrct[idx].wvlim_orig[0,0]/civstrct[idx].wrest[0]-1.,$
                    civstrct[idx].wvlim_orig[0,1]/civstrct[idx].wrest[0]-1.,$
                    civstrct[idx].wvlim_orig[1,0]/civstrct[idx].wrest[1]-1.,$
                    civstrct[idx].wvlim_orig[1,1]/civstrct[idx].wrest[1]-1.
;           stop
        endif
          
        if keyword_set(debug) and ndup gt 1 then begin
           ;; _extra= includes e.g., dblt_name=, altdblt_name=,
           ;; /label, initials= 
           sdss_chkciv,civstrct[idx],'sdss_combineciv.fit',_extra=extra
           tmpstr = xmrdfits('sdss_combineciv.fit',1,/silent)
           spawn,'rm sdss_combineciv.fit*'
        endif                   ; /debug
 
     endif ; ndup > 0

  endfor                        ; loop qq=nqso

  ;; dvnext_all are all positive or NaN
  if keep[0] ne -1 then begin
     ncomplement = (size(icomplement,/dim))[0] > 1
     complement = civstrct[icomplement]
     mask = intarr(nciv)        ; 0 = keep; 1 = duplicate
     mask[icomplement] = 1
     index = where(mask eq 0,count)

     if keyword_set(debug) then begin
        ;; _extra= includes e.g., dblt_name=, altdblt_name=,
        ;; /label, initials= 
        strct = sdss_srtcivstrct([civstrct[keep],civstrct[icomplement]])
        sdss_chkciv,strct,'sdss_combineciv.fit',_extra=extra
        tmpstr = xmrdfits('sdss_combineciv.fit',1,/silent)
        spawn,'rm sdss_combineciv.fit*'

        ;; Print useful information
        print,'Comparing old to new redshifts'
        print,'Min/Max change ',$
              min(civstrct0[keep].zabs_orig[0]-civstrct[keep].zabs_orig[0],max=mx),mx
        print,'Median change (abs)',$
              median(abs(civstrct0[keep].zabs_orig[0]-civstrct[keep].zabs_orig[0]))

        stop,'sdss_combineciv() debug stop: pausing'

     endif                      ; /debug
  endif else begin
     ;; Must instantiate
     count = 0
     index = lindgen(nciv)      ; all remain
     icomplement = -1
     ncomplement = -1
  endelse 
  subcivstr = civstrct[index]

  return, subcivstr

end                             ; sdss_combineciv()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_runcombineciv, allciv_fil, out_fil, civstrct_fil=civstrct_fil, _extra=extra
  ;; This isn't a really generic routine but gives information
  ;; on combining systems
  ;; _extra= includes dvsys=, /debug and stuff for sdss_chkciv

  allstr = xmrdfits(allciv_fil,1,/silent)
  dupstr0 = xmrdfits(allciv_fil,2,/silent)
  
  ;; Only want to combine the good and definite ratings (including BAL LOS)
  gdall = where(allstr.rating[0] ge sdss_getrating(/good),complement=bdall)
  teststr = sdss_combineciv(allstr[gdall],index=suball,icomplement=dupall,$
                            _extra=extra)

  ;; Combine and sort
  newallstr = sdss_srtcivstrct([allstr[bdall],teststr])

  ;; Test
  civstrct0 = sdss_getcivstrct(civstrct_fil,/default)
  subcivstr = sdss_combineciv(civstrct0,index=sub,icomplement=dup,dvsys=250.)
  gd = where(newallstr.rating[0] ge sdss_getrating(/good) and $
             (newallstr.balflg and 12) ne 12)
  test = where(newallstr[gd].qso_name ne subcivstr.qso_name or $
               newallstr[gd].zabs_orig[0] ne subcivstr.zabs_orig[0] or $
               newallstr[gd].rating[0] ne subcivstr.rating[0],ntest)
  if ntest ne 0 then $
     stop,'sdss_runcombineciv stop: failed sanity check'

  mwrfits,newallstr,out_fil,/create,/silent
  mwrfits,dupstr0,out_fil,/silent
  mwrfits,allstr[gdall[dupall]],out_fil,/silent

  ;; Measure EW
;  sdss_ewciv_strct,newallstr,'test_newallstr.fit'


end                             ; sdss_runcombineciv


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_funstats, nciv, ncolm=ncolm, density=density, dqso=dqso
  if n_params() ne 1 then begin
     print,'Syntax - sdss_funstats, nciv [ncolm=, density=, dqso=]'
     return
  endif 
  if not keyword_set(ncolm) then ncolm = 1.e13 ; atom/cm^2
  if not keyword_set(density) then density = 1. ; g/cm^3
  if not keyword_set(dqso) then dqso = 1. ; pc
  pc_to_cm = 30.857e18                       ; cm
  cm3_to_km3 = 1.e-15                     ; km^3
  mcarbon = 2.e-23                             ; g per atom
  rearth = 6371.                ; km
  mearth = 5.9742e24            ; kg
  rsun = 6.96e5                 ; km
  msun = 1.98892e30             ; kg

  ;; Assumptions
  print,''
  print,''
  print,'Assume QSO diameter (pc):',dqso
  print,'Assume average column (cm^-2):',ncolm
  print,'Assume density (g cm^-3):',density
  print,''
  
  ;; mass 
  ;; system * atom /(cm^2 * system) * cm^2 * g / atom
  mass = nciv * ncolm * dqso^2 * mcarbon * pc_to_cm * pc_to_cm ; g
  print,'Total mass (g):',mass
  print,'Total mass (Mearth):',mass / (1000. * mearth)
  print,'Total mass (Msun):',mass / (1000. * msun)

  ;; volume
  volume = mass / density * cm3_to_km3 ; km^3
  print,'Volume (km^3):',volume

  ;; box size
  dbox = volume^(1./3.)
  print,'Box size (km):',dbox
  print,'Box size (Rearth):',dbox / rearth
  print,'Box size (Rsun):',dbox / rsun
  print,''

end                             ; sdss_funstats


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_geturl, qso_name, sdsstab=sdsstab, _extra=extra
  if n_params() ne 1 then begin
     print,'Syntax - sdss_geturl(qso_name,[sdsstab=, _extra=])'
     return,-1
  endif 

  url_root = 'http://cas.sdss.org/dr7/en/tools/explore/obj.asp?ra=#&dec=#'
  prs = strsplit(url_root,'#',/extract,count=nprs)

  if keyword_set(sdsstab) then begin
     if size(sdsstab,/type) eq 7 then $
        sdsssum = xmrdfits(sdsstab,1,/silent) $
     else sdsssum = sdsstab
  endif else sdsssum = sdss_getqsostrct(_extra=extra) ; incl /nobal
  file = sdss_getname(sdsssum,root=qso_name_all)

  nqso = (size(qso_name,/dim))[0] > 1
  url = strarr(nqso)

  for qq=0L,nqso-1 do begin
     mtch = where(qso_name[qq] eq qso_name_all,nmtch)
     if nmtch ne 1 then stop,'sdss_showsite stop: no match ',qso_name[qq]
     url[qq] = prs[0] + strtrim(sdsssum[mtch].ra,2) + $
               prs[1] + strtrim(sdsssum[mtch].dec,2) 
  endfor                        ; loop qq=nqso

  return, url
end                             ; sdss_geturl()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_applycivobscorr, cmpltstr_fil, civobs_corr, $
                               rec_param=rec_param, _extra=extra
  ;;;;;;;;;
  ;; Factor in the effect of pathlength blocked by actual CIV systems 
  ;; Since did *not* remove any/all during completeness tests, real
  ;; CIV systems will artificially allow any weaker simulated profile
  ;; laid on top
  if size(cmpltstr_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmpltstr_fil,1,/silent) $
  else cmpltstr = cmpltstr_fil

  typ = size(civobs_corr,/type)
  if typ eq 2L then begin
     ;; Assume want to use default and whatever else
     civstr = sdss_getcivstrct(/default,_extra=extra,count=nciv)
  endif else begin
     ;; it's a file
     ;; _extra includes /default, rating=, dvqso=, /noBAL, /final, dvem=
     ;; ...
     civstr = sdss_getcivstrct(civobs_corr,_extra=extra,count=nciv)
  endelse
  tags = tag_names(civstr)
  if keyword_set(rec_param) then begin
     ;; I'm going to force this to be true because it's fair
     ;; And hopefully _extra doesn't have /final in it
     ;; (really should have made this tag business a function long
     ;; ago) 
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endif else begin
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endelse 
  
  ;; Loop over redshift bins and get the new curve
  ewindx = lindgen(cmpltstr.newbin)
  ewabs = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  ewabs = ewabs * 1.01d
  
  civstr_lo = civstr
  civstr_lo.(ewtag)[0] = civstr.(ewtag)[0] - civstr.(sigewtag)[0]
  civstr_hi = civstr
  civstr_hi.(ewtag)[0] = civstr.(ewtag)[0] + civstr.(sigewtag)[0]

  newcmpltstr = cmpltstr
  for zz=0,cmpltstr.nzbin-1 do begin
     ;; Set the grid redshifts plus a little buffer to avoid
     ;; rounding problems
     if (size(cmpltstr.zlim,/n_dim))[0] gt 1 then $
        zabs = replicate(mean(cmpltstr.zlim[zz,*]), cmpltstr.newbin) $
     else zabs = replicate(mean(cmpltstr.zlim), cmpltstr.newbin)
     zabs = zabs * 1.01d
     ;; Using /grid will force sdss_getdxw() to return
     ;; exactly what's in cmpltstr.czw_2darr in the absence
     ;; of civstrct_fil.
     ;; _extra includes ewmax=, dxwmax=, /final, dxciv= and
     ;; whatever its _extra goes to
     indx = zz + ewindx * cmpltstr.nzbin
     newcmpltstr.czw_2darr[indx] = $
        sdss_getdxw(cmpltstr, zabs, ewabs, sigma=signewczw, $
                    /czw, civobs_corr=civstr, silent=silent, $
                    sigewabs=0, /grid, _extra=extra)
     ;; Estimate error by seeing change when modify civstr EW and
     ;; add in quadrature to the error.
     ;; The /grid option may be funny because some
     ;; civstr.(sigewtag)[0] may keep it within same grid, which I
     ;; guess is ok.
     ;; Subtracting sigEW actually increases C(z,W) and adding
     ;; sigEW decreases C(z,W) so the below names just refer to the
     ;; sigEW sign.
     newczw_lo = sdss_getdxw(cmpltstr, zabs, ewabs, $
                             /czw, civobs_corr=civstr_lo, silent=silent, $
                             sigewabs=0, /grid, _extra=extra)
     newczw_hi = sdss_getdxw(cmpltstr, zabs, ewabs, $
                             /czw, civobs_corr=civstr_hi, silent=silent, $
                             sigewabs=0, /grid, _extra=extra)
     
     ;; some useful plotting checks
;     x_splot,ewabs,cmpltstr.czw_2darr[indx],ytwo=newcmpltstr.czw_2darr[indx],psym1=10,psym2=10,/block,ythr=newczw_lo,psym3=10,yfou=newczw_hi,psym4=10,lgnd=['Original','New','-1 sigEW','+1 sigEW'],title='CIV Obs Corr',xtitle='EW',ytitle='C(W)'
;     x_splot,ewabs,sqrt(signewczw[*,0]^2 + (newcmpltstr.czw_2darr[indx] - newczw_lo)^2)/signewczw[*,0],psym1=10,ytwo=sqrt(signewczw[*,1]^2 + (newcmpltstr.czw_2darr[indx] - newczw_hi)^2)/signewczw[*,1],psym2=10,/block,lgnd=['New high err','New low err'],title='CIV Obs Corr',xtitle='EW',ytitle='Error'
;     x_splot,ewabs,newcmpltstr.czw_2darr[indx]/cmpltstr.czw_2darr[indx],psym1=10,/block,ltitle='CIV Obs Corr',xtitle='EW',ytitle='New/Old Ratio'
     
     ;; About how *_hi and *_lo are used here, see above
     newcmpltstr.sigczw_2darr[indx,0] = $
        sqrt(signewczw[*,0]^2 + $
             (newcmpltstr.czw_2darr[indx] - newczw_hi)^2)
     newcmpltstr.sigczw_2darr[indx,1] = $
        sqrt(signewczw[*,1]^2 + $
             (newczw_lo - newcmpltstr.czw_2darr[indx])^2)
  endfor                        ; loop zz=cmpltstr.nzbin
  

  return, newcmpltstr
end                             ; sdss_applycivobscorr()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_applyuserbias, cmpltstr_fil, biasuser_fil, _extra=extra
  ;;;;;;;;;
  ;; Factor in the effect of user bias as measured by
  ;; sdss_completeness_userbias 
  if size(cmpltstr_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmpltstr_fil,1,/silent) $
  else cmpltstr = cmpltstr_fil
  newcmpltstr = cmpltstr

  if size(biasuser_fil,/type) eq 7 then $
     userbiasstr = xmrdfits(biasuser_fil,1,/silent) $
  else userbiasstr = biasuser_fil
  nzdim_user = size(userbiasstr.zlim,/n_dim)
  if nzdim_user eq 1 then nzbin_user = 1 $
  else $
     nzbin_user = (size(userbiasstr.zlim[*,0],/dim))[0] 
  
  ;; Loop over redshift bins and get the new curve
  ewindx = lindgen(cmpltstr.newbin)
  ewabs = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  ewabs = ewabs * 1.01d
  
  for zz=0,cmpltstr.nzbin-1 do begin

     if cmpltstr.nzbin eq 1 then zlim = cmpltstr.zlim $
     else zlim = reform(cmpltstr.zlim[zz,*])
     user_indx = sdss_getuserbiasindex(userbiasstr,zlim,tp_coeff=coeff)

     czw_user = sdss_calcexpgrowthfunc(ewabs, coeff, sigma=sigczw_user)
     
;     ;; Don't extrapolate to larger than fit over
;     bd = where(ewabs lt userbiasstr.ewlim[0] or $
;                ewabs ge userbiasstr.ewlim[1])
;     if bd[0] ne -1 then begin
;        czw_user[bd] = 1.
;        sigczw_user[bd,*] = 0.
;     endif 
     bd = where(czw_user gt 1.)
     if bd[0] ne -1 then stop,'sdss_applyuserbias() stop: C_user > 1'

     ;; Multiply because czw_basic = Nrec / Ninput and czw_user =
     ;; Naccept / Nrec
     ;; so now newczw = Naccept / Ninput
     ;; Estimate error by basic error propagation
     indx = zz + ewindx * cmpltstr.nzbin

     newcmpltstr.czw_2darr[indx] = cmpltstr.czw_2darr[indx] * czw_user
     newcmpltstr.sigczw_2darr[indx,0] = $
        sqrt((cmpltstr.sigczw_2darr[indx,0]*czw_user)^2 + $
             (sigczw_user[*,0]*cmpltstr.czw_2darr[indx])^2)
     newcmpltstr.sigczw_2darr[indx,1] = $
        sqrt((cmpltstr.sigczw_2darr[indx,1]*czw_user)^2 + $
             (sigczw_user[*,1]*cmpltstr.czw_2darr[indx])^2)


;     ;; some useful plotting checks
;    x_splot,ewabs,cmpltstr.czw_2darr[indx],ytwo=czw_user,psym1=10,psym2=10,/block,ythr=newcmpltstr.czw_2darr[indx],psym3=10,lgnd=['Basic','User','New'],xtitle='EW',ytitle='C(W)',title='User Bias Corr'
;    x_splot,ewabs,czw_user/cmpltstr.czw_2darr[indx],psym1=10,xtitle='EW',ytitle='User/Basic Ratio',/block,title='User Bias Corr'
  endfor                        ; loop zz=cmpltstr.nzbin
  
  ;; Add information to structure
  newcmpltstr = create_struct(cmpltstr, $
                              'ZLIM_TP', userbiasstr.zlim, $
                              'COEFF_TP', userbiasstr.coeff_tp)
  return, newcmpltstr
end                             ; sdss_applyuserbias()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_adjafpfxwfit, fitstrct, errstrct, newerrstrct=newerrstrct, $
                            file=file, _extra=extra
  if (n_params() ne 2 and not keyword_set(file)) or $
     (n_params() ne 1 and keyword_set(file)) then begin
     print,'Syntax - sdss_adjafpfxwfit( fitstrct, {errstrct or /file}, [newerrstrct=, _extra])'
     return,-1
  endif 
  ;; Assumer error in dN/dX (denom) is neglible (and b/c
  ;; can't quite compute dN/dX error from fit)
  ;; scale = (1 - (dN/dX_afp)/(dN/dX))
  ;; vars(scale) = var(dN/dX_afp)/(dN/dX)^2 + 
  ;;               var(dN/dX)*((dN/dX_afp)/(dN/dX)^2)^2  --> O(0)
  ;; Then re-scaling of f(W) or k looks like:
  ;; X = scale * X_0
  ;; sig(X)[0] = sqrt( var(scale)[1] * X_0^2 + 
  ;;                   scale^2 * sig(X_0)[0]^2)
  ;; sig(X)[1] = sqrt( var(scale)[0] * X_0^2 + 
  ;;                   scale^2 * sig(X_0)[1]^2)
  ;; because scale is larger if dN/dX_afp is smaller, and a
  ;; larger scale leads to less correction (larger X)
  ;; Correct a given f(W) fit structure from sdss_fitfxw.pro
  if keyword_set(file) then begin
     fitstrct = xmrdfits(fitstrct,1,/silent)
     errstrct = xmrdfits(fitstrct,2,/silent)
  endif 
  newfitstrct = create_struct(fitstrct,'coeff0',fitstrct.coeff,$
                              'sigcoeff0',fitstrct.sigcoeff)

  ;; Calculate dN/dX just over integration limits and data limits
  ;; zrng = [min, median, max], intlim = [min, sat, max]
  ;; _extra includes /dz and stuff for sdss_calcdxzw()
  afpdndx = sdss_getafpdndx(fitstrct.zrng[[0,2]], fitstrct.intlim[[0,2]], $
                            user_indx=user_indx,_extra=extra)

  ;; Using fit dN/dX integrated over the region where f(W) fit
  scale = 1. - afpdndx.dndx[0]/fitstrct.dndx
  newfitstrct.coeff = scale * fitstrct.coeff
  newfitstrct.sigcoeff[0] = $
     (afpdndx.sigdndx[0,1]/fitstrct.dndx * fitstrct.coeff)^2 + $
     (fitstrct.sigcoeff[0]*afpdndx.dndx[0]/fitstrct.dndx^2 * fitstrct.coeff)^2 + $
     (fitstrct.sigcoeff[0] * scale)^2 + $
     2*fitstrct.sigcoeff[0]*fitstrct.sigcoeff[0] * $
     (afpdndx.dndx[0]/fitstrct.dndx^2 * scale * fitstrct.coeff)
  newfitstrct.sigcoeff[1] = $
     (afpdndx.sigdndx[0,0]/fitstrct.dndx * fitstrct.coeff)^2 + $
     (fitstrct.sigcoeff[1]*afpdndx.dndx[0]/fitstrct.dndx^2 * fitstrct.coeff)^2 + $
     (fitstrct.sigcoeff[1] * scale)^2 + $
     2*fitstrct.sigcoeff[1]*fitstrct.sigcoeff[1] * $
     (afpdndx.dndx[0]/fitstrct.dndx^2 * scale * fitstrct.coeff)
  newfitstrct.sigcoeff = sqrt(newfitstrct.sigcoeff)

  ;; Stretch the coefficient grid so that < and > the
  ;; best-fit values allow the 1-sigma surface to reach the
  ;; new errors
  ;; Basically, *knowing* it's uniform binning currently,
  ;; want the same number of bins to span the 1-sigma space
  ;; dcoeff = sig(k)/sig(k_0)*dcoeff_0
  ;; k_i = i*dcoeff + k_best
  ;;     = sig(k)/sig(k_0)*i*dcoeff_0 + scale*k_0,best
  ;;     = sig(k)/sig(k_0)*(k_i,0-k_best,0) + k_best
  ;; Slight differential stretching b/c of asymmetric errors
  newerrstrct = create_struct(errstrct,'afp_flg',1,'coeff_grid0',$
                              errstrct.coeff_grid) ; fitstrct already has

  ;; Stretching it distorts the iso-dN/dX relation but just shifting
  ;; it 
  ;; self-consistently find grid point of interest, don't make
  ;; it disjoint 
  nbin_coeff = (size(errstrct.coeff_grid,/dim))[0] > 1 
  logl_mx = max(errstrct.logl,imx)
  if logl_mx ne fitstrct.logl then $
     print,'sdss_adjafpfxwfit(): maximum likelihood values do not match',$
           fitstrct.logl,logl_mx
  icoeff_best = imx mod nbin_coeff ; nbin_coeff/2    ; by fiat
  newerrstrct.coeff_grid[0:icoeff_best-1] = $
     newfitstrct.sigcoeff[0]/fitstrct.sigcoeff[0]*$
     (errstrct.coeff_grid[0:icoeff_best-1]-fitstrct.coeff) + $
     newfitstrct.coeff
  newerrstrct.coeff_grid[icoeff_best:*] = $
     newfitstrct.sigcoeff[1]/fitstrct.sigcoeff[1]*$
     (errstrct.coeff_grid[icoeff_best:*]-fitstrct.coeff) + $
     newfitstrct.coeff
  
  return,newfitstrct
end                             ; sdss_adjafpfxwfit()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_expandcontistrct, npix, nlin=nlin
  ;; Expand sdsscontistrct__define pixels to allow for storage of the
  ;; stacking results
  if n_params() ne 1 then begin
     print,'sdss_expandcontistrct( npix, [nlin=])'
     return,-1
  endif 
  cstrct0 = { sdsscontistrct }  ; template
  tags = tag_names(cstrct0)
  ntags = (size(tags,/dim))[0] > 1
  nconti = (size(cstrct0.ncent,/dim))[0] > 1 ; # flavor conti

  ;; Know the first one is not going to change
  cstrct = create_struct(tags[0],cstrct0.(0)) ; instantiate
  for tt=1,ntags-1 do begin
     case tags[tt] of 
        'SNR_CONV': arr = fltarr(npix,nconti)
        'CONTI': arr = dblarr(npix,nconti)
        'SIGCONTI': arr = dblarr(npix,nconti)
        else: begin
           if keyword_set(nlin) then begin
              case tags[tt] of 
                 'CENTROID': arr = fltarr(nlin,nconti)
                 'WVLIM_ORIG': arr = fltarr(nlin,2)
                 'EW_ORIG': arr = fltarr(nlin)
                 'SIGEW_ORIG': arr = fltarr(nlin)
                 else: arr = cstrct0.(tt)
              endcase
           endif else arr = cstrct0.(tt)
        end
     endcase 
     cstrct = create_struct(cstrct,tags[tt],arr)
  endfor                        ; tt=1,ntags-1

  return, cstrct
end                             ; sdss_expandcontistrct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_expandqalstrct, nlin, qstrct0=qstrct0
  ;; Expand qalcharstrct__define number of lines abled to be stored
  ;; for use with high-res, higher-SNR spectra
  if n_params() ne 1 then begin
     print,'sdss_expandqalstrct(nlin)'
     return,-1
  endif
  if not keyword_set(qstrct0) then $
     qstrct0 = { qalcharstrct } ; template
  tags = tag_names(qstrct0)
  ntags = (size(tags,/dim))[0] > 1

  ;; Know the first one is not going to change
  qstrct = create_struct(tags[0],qstrct0.(0)) ; instantiate
  for tt=1,ntags-1 do begin
     case tags[tt] of
        'ZLIM2': arr = dblarr(nlin,2) ; added in sdss_fndciv
        else: begin
           if stregex(tags[tt],'DLA_Z',/boolean) then $
              arr = dblarr(nlin) $
           else begin
              if stregex(tags[tt],'DLA_',/boolean) then $
                 arr = fltarr(nlin) $
              else arr = qstrct0.(tt) ; unchaged ones
           endelse
        end
     endcase
     qstrct = create_struct(qstrct,tags[tt],arr)
  endfor                        ; tt=1,ntags-1
  return, qstrct
end

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_mkstacksumm, inp_fil, outfil=outfil, list=list, lin_fil=lin_fil, $
                           dwvtol=dwvtol, clobber=clobber
  ;; Compile all the information from a list of stacks
  if n_params() ne 1 then begin
     print,'Syntax - sdss_mkstacksumm( inp_fil, [outfil=, lin_fil=, dwvtol=,'
     print,'                          /list, /clobber])'
     return,-1
  endif 

  if not keyword_set(dwvtol) then dwvtol = 1. ; Ang
  if not keyword_set(lin_fil) then $
     lin_fil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls_stack.lst'
  linstr = x_setllst(lin_fil,0) ; lls.lst for
  nlin = (size(linstr,/dim))[0] > 1

  if keyword_set(list) then $
     readcol,inp_fil,stack_fil,format='a',/silent $
  else stack_fil = inp_fil      ; may be array of cstrct/abslin structures
  nfil = (size(stack_fil,/dim))[0] > 1

  ;; "Ave" tags are [mean,median]
  tmpltstr = {$
             STACK_FIL:'',MEDIAN:-1, PERCENTILE:fltarr(2),$
             NABS:0L,WVION:0.,$
             ZAVE:dblarr(2),ZLIM:dblarr(2),$
             EWAVE:fltarr(4),EWLIM:fltarr(4),$ ; if sdss_stackciv_jackknife output, indices 3,4 contain info of excluded sample
             ION:linstr.name,LSNR:0.,FVAL:linstr.fval,$ ; store oscillator strength
             WREST:linstr.wave,ZABS:fltarr(nlin),WVLIM:fltarr(nlin,2),$
             EW:fltarr(nlin),SIGEW:fltarr(nlin),$ ;EWMINMAX:fltarr(nlin,2),$
             NCOLM:fltarr(nlin),SIGNCOLM:fltarr(nlin) $ ; don't actually have these
             }    
  strct = replicate(tmpltstr,nfil)
  if size(stack_fil,/type) eq 7 then strct.stack_fil = stack_fil $
  else strct.stack_fil = 'Input Structure'

  for ff=0,nfil-1 do begin
     if size(stack_fil,/type) eq 7 then begin
        ;; Copy meta-information
        hdr = xheadfits(strct[ff].stack_fil)
        strct[ff].median = sxpar(hdr,'MEDIAN') ; 0: mean; 1: median
        strct[ff].percentile[0] = sxpar(hdr,'PERLOW')
        strct[ff].percentile[1] = sxpar(hdr,'PERHIGH')
        strct[ff].nabs = sxpar(hdr,'NABS')
        strct[ff].wvion = sxpar(hdr,'WVION')
        strct[ff].zave[0] = sxpar(hdr,'ZMEAN')
        strct[ff].zave[1] = sxpar(hdr,'ZMED')
        strct[ff].zlim[0] = sxpar(hdr,'ZMIN')
        strct[ff].zlim[1] = sxpar(hdr,'ZMAX')
        strct[ff].ewave[0] = sxpar(hdr,'EWMEAN')
        strct[ff].ewave[1] = sxpar(hdr,'EWMED')
        strct[ff].ewave[2] = sxpar(hdr,'EWAVE_JK') ; "EWMEAN_JK" too long
        strct[ff].ewave[3] = sxpar(hdr,'EWMED_JK')
        strct[ff].ewlim[0] = sxpar(hdr,'EWMIN')
        strct[ff].ewlim[1] = sxpar(hdr,'EWMAX')
        strct[ff].ewlim[2] = sxpar(hdr,'EWMIN_JK')
        strct[ff].ewlim[3] = sxpar(hdr,'EWMAX_JK')

        ;; Read in conti structure
        cstrct = xmrdfits(strct[ff].stack_fil,1,/silent)
;        if keyword_set(sxpar(hdr,'NITER')) then $
;           ;; for min/max EW if ran sdss_stackciv_errmc()
;           gstrct = xmrdfits(strct[ff].stack_fil,2,/silent) 
     endif else cstrct = stack_fil[ff]
     
     cindx = fix(alog(cstrct.cflg)/alog(2))
     strct[ff].lsnr = cstrct.snr_conv[cstrct.npix,cindx]

     ;; Only loop over ions possibly detectable
     rng = lindgen(cstrct.ncent[cindx])
     wvmn = min(cstrct.centroid[rng,cindx],max=wvmx)
     sub = where(strct[ff].wrest ge wvmn-dwvtol and strct[ff].wrest le wvmx+dwvtol,nsub)
     for ss=0,nsub-1 do begin
        mtch = where(abs(cstrct.centroid[rng,cindx]-strct[ff].wrest[sub[ss]]) lt dwvtol,nmtch)
        if nmtch eq 0 then continue ; no match
        if nmtch gt 1 then begin
           mn = min(cstrct.centroid[rng[mtch],cindx]-strct[ff].wrest[sub[ss]],imn,/abs)
           mtch = rng[mtch[imn]]
        endif else mtch = rng[mtch[0]]

        ;; Save results
        strct[ff].zabs[sub[ss]] = cstrct.centroid[mtch,cindx] / strct[ff].wrest[sub[ss]] - 1.
        strct[ff].wvlim[sub[ss],*] = cstrct.wvlim_orig[mtch,*] 
        strct[ff].ew[sub[ss]] = cstrct.ew_orig[mtch] / (1 + strct[ff].zabs[sub[ss]])
        strct[ff].sigew[sub[ss]] = cstrct.sigew_orig[mtch] / (1 + strct[ff].zabs[sub[ss]])
;        if keyword_set(gstrct) then begin
;           stop,'sdss_mkstacksumm() stop: currently EWMINMAX not calculable'
;        endif
     endfor                     ; ss=nsub
     
  endfor                        ; loop ff=nfil

  ;; Write (tested before)
  if keyword_set(outfil) then begin
     test = file_search(outfil+'*',count=ntest)
     if ntest ne 0 and not keyword_set(clobber) then $
        print,'sdss_mkstacksumm(): file exists; will not clobber ',outfil $
     else begin
        mwrfits,strct,outfil,/create,/silent
        spawn,'gzip -f '+outfil
        print,'sdss_mkstacksumm(): created ',outfil
     endelse
  endif 

  return, strct
end                             ; sdss_mkstacksumm()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_printratedsumm, strct_fil
  ;; Print final results of rating
  if n_params() ne 1 then begin
     print,'Syntax - sdss_printratedsumm, strct_fil'
     return
  endif
  
  if size(strct_fil,/type) eq 8 then $
     civstr = strct_fil $
  else civstr = xmrdfits(strct_fil,1,/silent) 
  print,''
  tmp = sdss_getrateddblt(civstr,/definite,count=ndef)
  tmp = sdss_getrateddblt(civstr,/good,count=ngd)
  tmp = sdss_getrateddblt(civstr,/maybe,count=nmb)
  tmp = sdss_getrateddblt(civstr,/bad,count=nbd)
  tmp = sdss_getrateddblt(civstr,/unrated,count=nun)
  tmp = where(civstr.balflg ne 0,nbal) ; sdss_getbalflg(/visual),nbal)
  tmp = where(civstr.rating[9] eq sdss_getblendflg(),nblnd) ; default
  print,'sdss_chkciv: Summary:'
  print,'Definite: ',ndef,format='(a12,1x,i6)'
  print,'Good: ',ngd,format='(a12,1x,i6)'
  print,'Maybe: ',nmb,format='(a12,1x,i6)'
  print,'Bad: ',nbd,format='(a12,1x,i6)'
  print,'Unrated: ',nun,format='(a12,1x,i6)'
  print,'BALs: ',nbal, format='(a12,1x,i6)'
  print,'Blended: ',nblnd, format='(a12,1x,i6)'
  print,'Total: ',n_elements(civstr),format='(a12,1x,i6)'
  print,''

end                             ; sdss_printratedsumm


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_getstackdat, stackstrct_fil, z_ion, ion, zrng=zrng, $
                           dztol=dztol, dwvtol=dwvtol, skip_null=skip_null, $
                           count=count, nosrt=nosrt, fnorm=fnorm
  ;; Either get all the ions at one redshift or all the redshifts for
  ;; one ion 
  if n_params() ne 3 then begin
     print,'Syntax - sdss_getstackdat( stackstrct_fil, z_ion, ion, [/zrng, '
     print,'                           dztol=, dwvtol=, /skip_null, count=, /nosrt, /fnorm])'
     return,-1
  endif 
  if size(stackstrct_fil,/type) eq 7 then $ ; from sdss_mkstacksumm()
     stackstr = xmrdfits(stackstrct_fil,1,/silent) $
  else stackstr = stackstrct_fil
  stackstr.ion = strtrim(stackstr.ion,2)
  nstack = (size(stackstr,/dim))[0] > 1
  nlin = (size(stackstr.ion,/dim))[0] > 1
  nz_ion = (size(z_ion,/dim))[0] > 1
  nion = (size(ion,/dim))[0] > 1

  if keyword_set(fnorm) then begin
     ;; Normalize all equivalent widths by oscillator strength
     for ss=0,nstack-1 do $
        stackstr[ss].ew *= 1/stackstr[ss].fval
  endif

  ;; Sanity ckeck on uniformity
  unq = uniq(stackstr.median)
  if n_elements(unq) ne 1 then begin
     print,'sdss_getstackdat(): ERROR!!! stacks not all mean or median; exiting.'
     return,-1
  endif
  if keyword_set(stackstr[0].median) then iave = 1 else iave = 0

  unqlo = uniq(stackstr.percentile[0])
  unqhi = uniq(stackstr.percentile[1])
  if n_elements(unqlo) ne 1 or n_elements(unqhi) ne 1 then begin
     print,'sdss_getstackdat(): ERROR!!! stacks not all same percentile; exiting.'
     return,-1
  endif
  percentile = stackstr[0].percentile ; conveyance of info

  if not keyword_set(dwvtol) then dwvtol = 1.e-2 ; Ang
  if not keyword_set(dztol) then dztol = 250./2.998e5 ; 250 km/s
  if nion gt 1 and nz_ion gt 1 then begin
     print,'sdss_getstackdat(): cannot retrieve all ion EWs at all z'
     return,-1
  endif
  
  if nz_ion gt 1 then begin
     ;; Return one ion EW at all redshifts
     ;; Assume all ions in same order
     if size(ion,/type) eq 7 then $
        mtch = where(stackstr[0].ion eq ion,nmtch) $
     else mtch = where(abs(stackstr[0].wrest-ion) lt dwvtol,nmtch)
     if nmtch ne 1 then $
        stop,'sdss_getstackdat() stop: less or more than one ion match ',ion,nmtch
     
     if keyword_set(zrng) then begin
        ;; All in range
        sub = where(stackstr.zave[iave] ge z_ion[0] and $
                    stackstr.zave[iave] lt z_ion[1],count)
     endif else begin
        ;; Specific redshifts
        mask = intarr(nstack)
        for zz=0,nz_ion-1 do begin
           gd = where(abs(stackstr.zave[iave]-z_ion[zz]) lt dztol)
           if gd[0] ne -1 then mask[gd]++
        endfor                  ; zz=nz_ion
        sub = where(mask ne 0,count)
     endelse

     if sub[0] eq -1 then stop,'sdss_getstackdat() stop: no matching redshifts'

     stackion = stackstr[sub].wvion ; don't lose what defines stack
     ewlim = fltarr(count,4,/nozero)
     ewlim[*,0] = stackstr[sub].ewlim[0]
     ewlim[*,1] = stackstr[sub].ewlim[1]
     ewlim[*,2] = stackstr[sub].ewlim[2]
     ewlim[*,3] = stackstr[sub].ewlim[3]
     ewabs = fltarr(count,2,/nozero) ; included and excluded
     ewabs[*,0] = stackstr[sub].ewave[iave]
     ewabs[*,1] = stackstr[sub].ewave[2+iave]
     wrest = stackstr[sub].wrest[mtch[0]]
     wvlim = stackstr[sub].wvlim[mtch[0],*]
     zbin = stackstr[sub].zave[iave]
     zabs = stackstr[sub].zabs[mtch[0]] ; more like a delta z
     xdat = zbin
     sigxdat = fltarr(count,2,/nozero)
     sigxdat[*,0] = xdat-stackstr[sub].zlim[0]
     sigxdat[*,1] = stackstr[sub].zlim[1]-xdat
     ydat = stackstr[sub].ew[mtch[0]]
     if count gt 1 then $
        sigydat = rebin(stackstr[sub].sigew[mtch[0]],count,2) $
     else sigydat = replicate(stackstr[sub].sigew[mtch[0]],1,2)
  endif else begin
     ;; Return set of ion EW at one redshift
     mtch = where(abs(stackstr.zave[iave]-z_ion) lt dztol,nmtch)
     if nmtch ne 1 then $
        stop,'sdss_getstackdat() stop: less or more than one z match ',z_ion,nmtch

     if stregex(ion,'all',/boolean,/fold_case) then begin
        sub = indgen(nlin) 
        count = nlin
     endif else begin
        mask = intarr(nlin)
        for ii=0,nion-1 do begin
           if size(ion,/type) eq 7 then $
              gd = where(stackstr[mtch[0]].ion eq ion[ii]) $
           else gd = where(abs(stackstr[mtch[0]].wrest-ion[ii]) lt dwvtol)
           if gd[0] ne -1 then mask[gd]++
        endfor 
        sub = where(mask ne 0,count)
     endelse 
     
     if sub eq -1 then stop,'sdss_getstackdat() stop: no matching ions'
     
     stackion = stackstr[mtch[0]].wvion ; don't lose what defines stack
     ewlim = fltarr(count,4,/nozero)
     ewlim[*,0] = stackstr[mtch[0]].ewlim[0]
     ewlim[*,1] = stackstr[mtch[0]].ewlim[1]
     ewlim[*,2] = stackstr[mtch[0]].ewlim[2]
     ewlim[*,3] = stackstr[mtch[0]].ewlim[3]
     ewabs = fltarr(count,2,/nozero)
     ewabs[*,0] = stackstr[mtch[0]].ewave[iave]
     ewabs[*,1] = stackstr[mtch[0]].ewave[2+iave]
     wrest = stackstr[mtch[0]].wrest[sub]
     zbin = stackstr[sub].zave[iave]
     zabs = stackstr[mtch[0]].zabs[sub] ; more like a delta z
     wvlim = stackstr[mtch[0]].wvlim[sub,*]
     xdat = wrest
     sigxdat = fltarr(count,2,/nozero)   ; no error
     ydat = stackstr[mtch[0]].ew[sub]
     sigydat = rebin(stackstr[mtch[0]].sigew[sub],nstack,2)
  endelse 
  
  if keyword_set(skip_null) then begin
     gd = where(ydat ne 0.,count)
     if count ne 0 then begin
        stackion = stackion[gd]
        ewlim = ewlim[gd,*]
        ewabs = ewabs[gd,*]
        wrest = wrest[gd]
        zbin = zbin[gd]
        wvlim = wvlim[gd,*]
        xdat = xdat[gd]
        sigxdat = sigxdat[gd,*]
        ydat = ydat[gd]
        sigydat = sigydat[gd,*]
     endif 
  endif 

  if keyword_set(nosrt) then srt = lindgen(count) $
  else srt = sort(xdat)
  ostrct = {stackion:stackion,$ ; rest wavelength
            median:iave, $      ; 0: mean; 1: median
            percentile:percentile, $ ; [low,high]
            fnorm:keyword_set(fnorm),$ ; 0: rest EWs; 1: scalled by f_osc
            zion:z_ion, ion:ion, mean:keyword_set(mean),$
            dztol:dztol, dwvtol:dwvtol, zrng:keyword_set(zrng), $
            ewlim:ewlim[srt,*], ewabs:ewabs[srt,*], wrest:wrest[srt], $
            zabs:zabs, zbin:zbin, wvlim:wvlim[srt,*], $
            xdat:xdat[srt], sigxdat:sigxdat[srt,*], $
            ydat:ydat[srt], sigydat:sigydat[srt,*]}

  return, ostrct
end                             ; sdss_getstackdat()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpZhuMenard2strct,zm13strct_fil,sdsstab,bal_fil=bal_fil,nvdisp=nvdisp,$
                                _extra=extra
  ;; Copy Zhu & Menard (2013) MgII & FeII structure to
  ;; sdsscivstrct format 
  ;; http://www.pha.jhu.edu/~gz323/Site/Download_Absorber_Catalog.html

  if n_params() ne 1 then begin
     print,'Syntax - sdss_cpZhuMenard2strct(zm13strct_fil,[sdsstab,bal_fil=,nvdisp=,_extra=])'
     return,-1
  endif
  
  if not keyword_set(nvdisp) then nvdisp = 1. ; Gaussian width

  if size(zm13strct_fil,/type) eq 7 then $
     zm13strct = xmrdfits(zm13strct_fil,1,/silent) $
  else zm13strct = zm13strct_fil
  nzm13 = (size(zm13strct,/dim))[0] > 1 

  if not keyword_set(bal_fil) then bal_fil = sdss_getqsostrct(/BAL)
  if size(balstrct,/type) eq 7 then $
     balstrct = xmrdfits(bal_fil,1,/silent) $
  else balstrct = bal_fil
  tmp = sdss_getname(balstrct,root=balqso_name)
  
  ;; Construct Zhu & Menard "qso_name" (JJJJJ-PPPP-FFF)
  mjd = strtrim(zm13strct.mjd,2) ; may not be spectro MJD

  plate = strtrim(zm13strct.plate,2)
  bd = where(strlen(plate) ne 4,nbd)
  while nbd ne 0 do begin
     plate[bd] = '0' + plate[bd]
     bd = where(strlen(plate) ne 4,nbd)
  endwhile 

  fiber = strtrim(zm13strct.fiber,2)
  bd = where(strlen(fiber) ne 3,nbd)
  while nbd ne 0 do begin
     fiber[bd] = '0' + fiber[bd]
     bd = where(strlen(fiber) ne 3,nbd)
  endwhile 

  zm13qso_name = mjd+'-'+plate+'-'+fiber
  zm13conti_fil = 'abslin/1d_26/'+plate+'/1d/spSpec-'+zm13qso_name+'-abslin.fit' ; may not exist

  ;; Sort 
  srt = sort(zm13qso_name)
  zm13strct = zm13strct[srt]
  zm13qso_name = zm13qso_name[srt]
  zm13unq = uniq(zm13qso_name)
  nzm13unq = (size(zm13unq,/dim))[0] > 1
  nabs_per_zm13 = [zm13unq[0]+1,(shift(zm13unq,-1)-zm13unq)[0:nzm13unq-2]]
  zm13tags = tag_names(zm13strct)

  ;; _extra= includes dr=, name=, /bal, /nobal
  if not keyword_set(sdsstab) then sdsstab = sdss_getqsostrct(_extra=extra)
  
  if size(sdsstab,/type) eq 7 then $
     sdsssum = xmrdfits(sdsstab,1,/silent) $
  else sdsssum = sdsstab
  spec_fil = sdss_getname(sdsssum,root=qso_name,dir=sdir)
  spec_fil = sdir+spec_Fil

  ;; Set up output
  newstrct = create_struct({sdsscivstrct},'z_qso_hw10',0.d)
  mgiistrct = replicate(newstrct,nzm13)
  mgiistrct.rating[0] = sdss_getrating(/def)
  dblt = dblt_retrieve('MgII')
  mgiistrct.wrest[0] = dblt.wvI
  iewmgii_I = (where(zm13tags eq 'REW_MGII_2796'))[0]
  ivdmgii_I = (where(zm13tags eq 'VDISP_MGII_2796'))[0]
  mgiistrct.wrest[1] = dblt.wvII
  iewmgii_II = (where(zm13tags eq 'REW_MGII_2803'))[0]
  ivdmgii_II = (where(zm13tags eq 'VDISP_MGII_2803'))[0]
  dblt = dblt_retrieve('FeII')  ; 2600, 2586
  mgiistrct.wrest[2] = dblt.wvI
  iewfeii_I = (where(zm13tags eq 'REW_FEII_2600'))[0]
  ivdfeii_I = (where(zm13tags eq 'VDISP_FEII_2600'))[0]
  mgiistrct.wrest[3] = dblt.wvII
  iewfeii_II = (where(zm13tags eq 'REW_FEII_2586'))[0]
  ivdfeii_II = (where(zm13tags eq 'VDISP_FEII_2586'))[0]
  mgiistrct.wrest[4] = 2382.765
  iewfeii_III = (where(zm13tags eq 'REW_FEII_2383'))[0]
  ivdfeii_III = (where(zm13tags eq 'VDISP_FEII_2383'))[0]
  mgiistrct.wrest[5] = 2374.4612
  iewfeii_IV = (where(zm13tags eq 'REW_FEII_2374'))[0]
  ivdfeii_IV = (where(zm13tags eq 'VDISP_FEII_2374'))[0]
  mgiistrct.wrest[6] = 2344.214 
  iewfeii_V = (where(zm13tags eq 'REW_FEII_2344'))[0]
  ivdfeii_V = (where(zm13tags eq 'VDISP_FEII_2344'))[0]
  mgiistrct.abslin_fil = zm13conti_fil

  ;; Big loop
  for qq=0L,nzm13unq-1 do begin
     if qq gt 0 then zm13rng = zm13unq[qq-1] + lindgen(nabs_per_zm13[qq]) + 1 $
     else zm13rng = lindgen(nabs_per_zm13[qq])

     ;; Try on name only
     mtch = where(zm13qso_name[zm13rng[0]] eq qso_name,nmtch)

     case nmtch of 
        0: begin
           ;; Find by distance
           gcirc,1,zm13strct[zm13rng[0]].ra*24/360.,zm13strct[zm13rng[0]].dec,$
                 sdsssum.ra*24/360.,sdsssum.dec,dist
           mindist = min(dist,mtch)
        end
        else: begin
           if nmtch eq 1 then begin
              mtch = mtch[0] 

              gcirc,1,zm13strct[zm13rng[0]].ra*24/360.,zm13strct[zm13rng[0]].dec,$
                    sdsssum[mtch].ra*24/360.,sdsssum[mtch].dec,mindist
           endif else begin
              ;; Find closest
              gcirc,1,zm13strct[zm13rng[0]].ra*24/360.,zm13strct[zm13rng[0]].dec,$
                    sdsssum[mtch].ra*24/360.,sdsssum[mtch].dec,dist
              mindist = min(dist,imn)
              mtch = mtch[imn]
           endelse

           ;; Sanity check #1
           if mindist gt 10. then begin
              ;; Find closest
              gcirc,1,zm13strct[zm13rng[0]].ra*24/360.,zm13strct[zm13rng[0]].dec,$
                    sdsssum.ra*24/360.,sdsssum.dec,dist
              mindist = min(dist,mtch)
              print,'sdss_cpZhuMenard2strct(): Over-riding qso_name match for closest ',$
                    zm13qso_name[zm13rng[0]],' '+qso_name[mtch],mindist
           endif
        end
     endcase


     ;; Sanity check #2
     if mindist gt 10. then begin
        print,'sdss_cpZhuMenard2strct(): WARNING! no quasar match within 10 arcsec ',$
             zm13qso_name[zm13rng[0]],mindist
        ;; Load up information from Zhu & Menard (2013) structure
        test = where(zm13qso_name[zm13rng[0]] eq qso_name)
        if test[0] ne -1 then stop
        mgiistrct[zm13rng].qso_name = zm13qso_name[zm13rng[0]]
        mgiistrct[zm13rng].mjd = zm13strct[zm13rng[0]].mjd
        mgiistrct[zm13rng].plate = zm13strct[zm13rng[0]].plate
        mgiistrct[zm13rng].fiber = zm13strct[zm13rng[0]].fiber
        mgiistrct[zm13rng].ra = zm13strct[zm13rng[0]].ra
        mgiistrct[zm13rng].dec = zm13strct[zm13rng[0]].dec
        mgiistrct[zm13rng].rmag = -999 ; let this be the warning flag
        ;; skipping SNR (below)
        mgiistrct[zm13rng].sdss_obs[0] = 'spectro/1d_26/'+plate[zm13rng]+'/1d/spSpec-'+$
                                     zm13qso_name[zm13rng]+'.fit'
        mgiistrct[zm13rng].z_qso = zm13strct[zm13rng[0]].zqso
        mgiistrct[zm13rng].z_qso_hw10 = zm13strct[zm13rng[0]].zqso
        ;; skipping BALFLG
        mgiistrct[zm13rng].balflg = 1024 ; Zhu & Menard Unique
     endif else begin
        ;; Load up the information from Schneider et al. (2010) structure
        mgiistrct[zm13rng].qso_name = qso_name[mtch]
        mgiistrct[zm13rng].mjd = sdsssum[mtch].smjd
        mgiistrct[zm13rng].plate = sdsssum[mtch].plate
        mgiistrct[zm13rng].fiber = sdsssum[mtch].fiber
        mgiistrct[zm13rng].ra = sdsssum[mtch].ra
        mgiistrct[zm13rng].dec = sdsssum[mtch].dec
        mgiistrct[zm13rng].rmag = sdsssum[mtch].rtmag
        ;; skipping SNR (below)
        mgiistrct[zm13rng].sdss_obs[0] = spec_fil[mtch]
        mgiistrct[zm13rng].z_qso = sdsssum[mtch].z
        mgiistrct[zm13rng].z_qso_hw10 = zm13strct[zm13rng[0]].zqso
        ;; BALFLG
        test = where(balqso_name eq qso_name[mtch[0]])
        if test[0] ne -1 then mgiistrct[zm13rng].balflg = 12 ; 8 and 4
     endelse
     mgiistrct[zm13rng].abslin_fil = 'abslin'+strmid(spec_fil[mtch],7,36)+'-abslin.fit'

     for ii=0,nabs_per_zm13[qq]-1 do begin
        mgiistrct[zm13rng[ii]].zabs_orig[0:6] = zm13strct[zm13rng[ii]].zabs
        mgiistrct[zm13rng[ii]].sigzabs_orig[0:6] = zm13strct[zm13rng[ii]].err_zabs

        mgiistrct[zm13rng[ii]].snr[0] = zm13strct[zm13rng[ii]].med_sdeviation_red
        mgiistrct[zm13rng[ii]].snr[1] = zm13strct[zm13rng[ii]].med_sdeviation_blue
        mgiistrct[zm13rng[ii]].snr[2] = zm13strct[zm13rng[ii]].spec_snr_median

        for jj=0,6 do begin
           
           case jj of
              0: begin 
                 iew = iewmgii_I
                 ivd = ivdmgii_I
              end
              1: begin
                 iew = iewmgii_II
                 ivd = ivdmgii_II
              end
              2: begin
                 iew = iewfeii_I
                 ivd = ivdfeii_I
              end
              3: begin
                 iew = iewfeii_II
                 ivd = ivdfeii_II
              end
              4: begin
                 iew = iewfeii_III
                 ivd = ivdfeii_III
              end
              5: begin
                 iew = iewfeii_IV
                 ivd = ivdfeii_IV
              end
              6: begin
                 iew = iewfeii_V
                 ivd = ivdfeii_V
              end
           end
           isigew = iew + 1

           mgiistrct[zm13rng[ii]].ew_orig[jj] = zm13strct[zm13rng[ii]].(iew)
           mgiistrct[zm13rng[ii]].sigew_orig[jj] = zm13strct[zm13rng[ii]].(isigew)

           mgiistrct[zm13rng[ii]].gwidth[jj] = zm13strct[zm13rng[ii]].(ivd)
        endfor                  ; loop jj=0,6

     endfor                     ; loop ii=nabs_per_zm13
     
  endfor                        ; loop qq=nzm13unq

  c = 299792.458                ; km/s
  cinv = 1./c
  for jj=0,6 do begin
     wvobs = mgiistrct.wrest[jj]*(1+mgiistrct.zabs_orig[jj]) 
     dwv = nvdisp*mgiistrct.gwidth[jj]*cinv*wvobs
     mgiistrct.wvlim_orig[jj,0] = wvobs - dwv
     mgiistrct.wvlim_orig[jj,1] = wvobs + dwv
  endfor                        ; loop jj=0,6
  
  ;; Polish off
  mgiistrct.cflg = sdss_getcflg(/hyb)

  return, mgiistrct
end                             ;  sdss_cpZhuMenard2strct()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpProchteretal2strct, p06tsv_fil,sdsstab,bal_fil=bal_fil,$
                                    debug=debug, _extra=extra
  ;; Copy pipe (|) seperated file from
  ;; http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJ/639/766
  ;; to sdsscivstrct format

  if n_params() ne 1 then begin
     print,'Syntax - sdss_cpProchteretal2strct(p06tsv_fil,[sdsstab,bal_fil=,/debug,_extra=])'
     return,-1
  endif 
  c = 299792.458                ; km/s

  ;; Should already be sorted 
  readcol,p06tsv_fil,id,jname,ra,dec,zqso,zmgii,ewr,delimiter='|',$
          format='(i,a,a,a,f,f,f)',/silent
  np06 = (size(id,/dim))[0] > 1
  p06unq = uniq(jname)
  np06unq = (size(p06unq,/dim))[0] > 1
  nabs_per_p06 = [p06unq[0]+1,(shift(p06unq,-1)-p06unq)[0:np06unq-2]]

  if not keyword_set(bal_fil) then bal_fil = sdss_getqsostrct(/BAL)
  if size(balstrct,/type) eq 7 then $
     balstrct = xmrdfits(bal_fil,1,/silent) $
  else balstrct = bal_fil
  tmp = sdss_getname(balstrct,root=balqso_name)
  
  ;; _extra= includes dr=, name=, /bal, /nobal
  if not keyword_set(sdsstab) then sdsstab = sdss_getqsostrct(_extra=extra)
  if size(sdsstab,/type) eq 7 then $
     sdsssum = xmrdfits(sdsstab,1,/silent) $
  else sdsssum = sdsstab
  spec_fil = sdss_getname(sdsssum,root=qso_name,dir=sdir)
  spec_fil = sdir + spec_fil
  maxnamelen = strlen(qso_name[0])

  mgiistrct = replicate({sdsscivstrct},np06)
  mgiistrct.rating[0] = sdss_getrating(/def)
  dblt = dblt_retrieve('MgII')
  
  mgiistrct.z_qso = zqso        ; may be overwritten
  mgiistrct.wrest[0] = dblt.wvI
  mgiistrct.zabs_orig[0] = zmgii
  mgiistrct.ew_orig[0] = ewr    ; EW(2803) isn't known
  mgiistrct.wrest[1] = dblt.wvII
  mgiistrct.zabs_orig[1] = zmgii

  ;; Figure out RA and DEC in degrees
  strput,ra,':',2
  strput,ra,':',5
  strput,dec,':',3
  strput,dec,':',6
  x_radec,ra,dec,rad,decd
  mgiistrct.ra = rad            ; may be overwritten
  mgiistrct.dec = decd


  ;; Big loop
  for qq=0L,np06unq-1 do begin
     if qq gt 0 then p06rng = p06unq[qq-1] + lindgen(nabs_per_p06[qq]) + 1 $
     else p06rng = lindgen(nabs_per_p06[qq])

     gcirc,1,rad[p06unq[qq]]*24/360.,decd[p06unq[qq]],$
           sdsssum.ra*24/360.,sdsssum.dec,dist
     mindist = min(dist,imn)

     if mindist lt 10. then begin
        ;; Store all the information
        mgiistrct[p06rng].qso_name = qso_name[imn]
        mgiistrct[p06rng].ra = sdsssum[imn].ra
        mgiistrct[p06rng].dec = sdsssum[imn].dec
        ;; Sanity check
        dvqso = c*(mgiistrct[p06rng[0]].z_qso-sdsssum[imn].z)
        if abs(dvqso) gt 1000. and keyword_set(debug) then $
           print,'sdss_cpProchteretal2strct(): |dvqso| > 1000 km/s: ',qso_name[imn],$
                 sdsssum[imn].z, mgiistrct[p06rng[0]].z_qso,dvqso
        mgiistrct[p06rng].z_qso = sdsssum[imn].z
        mgiistrct[p06rng].mjd = sdsssum[imn].smjd
        mgiistrct[p06rng].plate = sdsssum[imn].plate
        mgiistrct[p06rng].fiber = sdsssum[imn].fiber
        mgiistrct[p06rng].rmag = sdsssum[imn].rtmag
        mgiistrct[p06rng].sdss_obs[0] = spec_fil[imn]
        mgiistrct[p06rng].abslin_fil = 'abslin'+$
                                       strmid(mgiistrct[p06rng].sdss_obs[0],7,36)+$
                                       '-abslin.fit'
        ;; BALFLG
        test = where(balqso_name eq qso_name[imn])
        if test[0] ne -1 then mgiistrct[p06rng].balflg = 12 ; 8 and 4
     endif else begin
        print,'sdss_cpProchteretal2strct(): WARNING! no quasar match within 10 arcsec ',$
              jname[p06rng[0]],mindist
        mgiistrct[p06rng].qso_name = strmid(jname[p06rng],0,maxnamelen) ; avoids having to strtrim later
        mgiistrct[p06rng].rmag = -999 ;  let this be the warning flag
        ;; Skipping BALFLG
     endelse

  endfor                        ; loop qq=np06unq

  return, mgiistrct
  
end                             ; sdss_cpProchteretal2strct()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_cpQuideretal2strct, q11tab_fil, sdsstab, bal_fil=bal_fil, debug=debug, $
                                  _extra=extra

  ;; Penn State absorption line catalog last described in Quider et
  ;; al. (2011)
  if n_params() ne 1 then begin
     print,'Syntax - sdss_cpQuideretal2strct(q11tab_fil,[sdsstab,bal_fil=,/debug,_extra=])'
     return,-1
  endif 
  c = 299792.458                ; km/s
  
  readcol,q11tab_fil,jname,mjd,plate,fiber,zmgii,$
          ew2796,sigew2796,ew2803,sigew2803,$
          ew2586,sigew2586,ew2600,sigew2600,$
          ewmnii2576,sigewmnii2576,ewmnii2594,sigewmnii2594,$
          ewmnii2606,sigewmnii2606,ewmgi2852,sigewmgi2852,$
          format='(a,a,a,a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)',/silent
  nq11 = (size(jname,/dim))[0] > 1
  q11unq = uniq(jname)
  nq11unq = (size(q11unq,/dim))[0] > 1
  nabs_per_q11 = [q11unq[0]+1,(shift(q11unq,-1)-q11unq)[0:nq11unq-2]]
  mgiistrct = replicate({sdsscivstrct},nq11)
  mgiistrct.rating[0] = sdss_getrating(/def)
  dblt = dblt_retrieve('MgII')

  ;; Parse the name
  ratmp = strmid(jname,1,9)
  dectmp = strmid(jname,10)
  rastr = strmid(ratmp,0,2)+':'+strmid(ratmp,2,2)+':'+strmid(ratmp,4)
  decstr = strmid(dectmp,0,3)+':'+strmid(dectmp,3,2)+':'+strmid(dectmp,5)
  x_radec,rastr,decstr,ra,dec

  ;; Make one very large array
  lines = [dblt.wvI,dblt.wvII,2586.6500,2600.1729,2576.877,2594.499,2606.462,2852.9642]
  nlines = (size(lines,/dim))[0] > 1
  tmp = [ew2796,sigew2796,ew2803,sigew2803,$
         ew2586,sigew2586,ew2600,sigew2600,$
         ewmnii2576,sigewmnii2576,ewmnii2594,sigewmnii2594,$
         ewmnii2606,sigewmnii2606,ewmgi2852,sigewmgi2852]
  dat = reform(tmp,nq11,16)

  q11qso_name = mjd+'-'+plate+'-'+fiber
  q11spec_fil = 'spectro/1d_26/'+plate+'/1d/spSpec-'+q11qso_name+'.fit' ; may not exist
  q11conti_fil = 'abslin/1d_26/'+plate+'/1d/spSpec-'+q11qso_name+'-abslin.fit' ; may not exist

  if not keyword_set(bal_fil) then bal_fil = sdss_getqsostrct(/BAL)
  if size(balstrct,/type) eq 7 then $
     balstrct = xmrdfits(bal_fil,1,/silent) $
  else balstrct = bal_fil
  tmp = sdss_getname(balstrct,root=balqso_name)
  
  ;; _extra= includes dr=, name=, /bal, /nobal
  if not keyword_set(sdsstab) then sdsstab = sdss_getqsostrct(_extra=extra)
  if size(sdsstab,/type) eq 7 then $
     sdsssum = xmrdfits(sdsstab,1,/silent) $
  else sdsssum = sdsstab
  spec_fil = sdss_getname(sdsssum,root=qso_name,dir=sdir)
  spec_fil = sdir + spec_fil

  ;; Load up structure 
  for ii=0,nlines-1 do begin
     mgiistrct.zabs_orig[2*ii] = zmgii

     mgiistrct.wrest[ii] = lines[ii]

     mgiistrct.ew_orig[ii] = dat[*,2*ii]
     mgiistrct.sigew_orig[ii] = dat[*,2*ii+1]
  endfor                        ; loop ii=nlines

  ;; ID and map sightlines
  for qq=0L,nq11unq-1 do begin
     if qq gt 0 then q11rng = q11unq[qq-1] + lindgen(nabs_per_q11[qq]) + 1 $
     else q11rng = lindgen(nabs_per_q11[qq])

     ;; Try on name only
     mtch = where(q11qso_name[q11rng[0]] eq qso_name,nmtch)

     case nmtch of 
        0: begin
           ;; Find by distance
           gcirc,1,ra[q11rng[0]]*24/360.,dec[q11rng[0]],$
                 sdsssum.ra*24/360.,sdsssum.dec,dist
           mindist = min(dist,mtch)
        end
        else: begin
           if nmtch eq 1 then begin
              mtch = mtch[0] 

              gcirc,1,ra[q11rng[0]]*24/360.,dec[q11rng[0]],$
                    sdsssum[mtch].ra*24/360.,sdsssum[mtch].dec,mindist
           endif else begin
              ;; Find closest
              gcirc,1,ra[q11rng[0]]*24/360.,dec[q11rng[0]],$
                    sdsssum[mtch].ra*24/360.,sdsssum[mtch].dec,dist
              mindist = min(dist,imn)
              mtch = mtch[imn]
           endelse

           ;; Sanity check #1
           if mindist gt 10. then begin
              ;; Find closest
              gcirc,1,ra[q11rng[0]]*24/360.,dec[q11rng[0]],$
                    sdsssum.ra*24/360.,sdsssum.dec,dist
              mindist = min(dist,mtch)
              print,'sdss_cpQuideretal2strct(): Over-riding qso_name match for closest ',$
                    q11qso_name[q11rng[0]],' '+qso_name[mtch],mindist
           endif
        end
     endcase

     ;; Sanity check #2
     if mindist gt 10. then begin
        print,'sdss_cpQuideretal2strct(): WARNING! no quasar match within 10 arcsec ',$
             q11qso_name[q11rng[0]],mindist
        ;; Load up information from Quider et al (2011)
        mgiistrct[q11rng].qso_name = q11qso_name[q11rng]
        mgiistrct[q11rng].mjd = mjd[q11rng]
        mgiistrct[q11rng].plate = plate[q11rng]
        mgiistrct[q11rng].fiber = fiber[q11rng]
        mgiistrct[q11rng].ra = ra[q11rng]
        mgiistrct[q11rng].dec = dec[q11rng]
        mgiistrct[q11rng].rmag = -999 ; let this be the warning flag
        ;; skipping SNR (below)
        mgiistrct[q11rng].sdss_obs[0] = q11spec_fil[q11rng]
        mgiistrct[q11rng].z_qso = -999
        ;; skipping BALFLG
     endif else begin
        ;; Load up the information from Schneider et al. (2010) structure
        mgiistrct[q11rng].qso_name = qso_name[mtch]
        mgiistrct[q11rng].mjd = sdsssum[mtch].smjd
        mgiistrct[q11rng].plate = sdsssum[mtch].plate
        mgiistrct[q11rng].fiber = sdsssum[mtch].fiber
        mgiistrct[q11rng].ra = sdsssum[mtch].ra
        mgiistrct[q11rng].dec = sdsssum[mtch].dec
        mgiistrct[q11rng].rmag = sdsssum[mtch].rtmag
        ;; skipping SNR (below)
        mgiistrct[q11rng].sdss_obs[0] = spec_fil[mtch]
        mgiistrct[q11rng].z_qso = sdsssum[mtch].z
        ;; BALFLG
        test = where(balqso_name eq qso_name[mtch[0]])
        if test[0] ne -1 then mgiistrct[q11rng].balflg = 12 ; 8 and 4
        mgiistrct[q11rng].abslin_fil = 'abslin'+strmid(spec_fil[mtch],7,36)+'-abslin.fit'
     endelse
     
  endfor                        ; loop qq=nq11unq
  
  return, mgiistrct
end                             ; sdss_cpQuideretal2strct()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_flux2ugriz, spec_fil, sdss3=sdss3, option=option, AB=AB
  ;; From a typical SDSS spectrum, estimate the magnitude
  if n_params() ne 1 then begin
     print,'Syntax - sdss_flux2ugriz( spec_fil, [/sdss3, option=, /AB]'
     return,-1
  endif

  if not keyword_set(option) then option = 'median'
  c =  299792.458d * 1e5        ; cm/s

  ;; http://www.sdss3.org/instruments/camera.php#Filters
;  mag_lbl = ['u', 'g', 'r', 'i', 'z']
;  mag_wave = transpose([[2980.,3551.,4130.], $ ; u
;                        [3630.,4686.,5830.], $ ; g
;                        [5380.,6166.,7230.], $ ; r 
;                        [6430.,7480.,8630.], $ ; i 
;                        [7730.,8932.,11230.]]) ; z
  nmag = 5 ; (size(mag_lbl,/dim))[0] > 1

  ;; One structure per extension per filter
  filter_fil = getenv('XIDL_DIR')+'/SDSS/General/SDSS_filter_curves.fits'
  
  ;; Both files store flux in units of 10^-17 erg/s/cm^2/Ang
  if keyword_set(sdss3) then begin
     dat = xmrdfits(spec_fil,1,/silent)
     wave = 10.^dat.loglam
     flux = dat.flux
     bd = where(dat.ivar le 0.,nbd)
     if nbd ne 0 then dat[bd].ivar = !values.f_nan
     sigma = sqrt(1./dat.ivar)  ; bad = NaN still
  endif else $
     parse_sdss,spec_fil,flux,wave,sig=sigma
  bdpix = where(sigma le 0 or finite(sigma,/nan) eq 1,nbdpix)
  if nbdpix ne 0 then sigma[bdpix] = -999.

  flux_nu = flux*1.e-17*wave^2/(c*1.e8) ; erg/cm^2/s/Hz
  sigflux_nu = sigma*1.e-17*wave^2/(c*1.e8) ; erg/cm^2/s/Hz

  
  magstrct = {spec_fil:spec_fil, $
              option:option, $ ; median, mean, wgtmean
              AB:keyword_set(AB), $ ; 1 = AB corrected, 0 = SDSS
              umag:0., umagerr:0., $
              gmag:0., gmagerr:0., $
              rmag:0., rmagerr:0., $
              imag:0., imagerr:0., $
              zmag:0., zmagerr:0}
  tags = tag_names(magstrct)
  i_umag = (where(tags eq 'UMAG'))[0]
  
  for mm=0,nmag-1 do begin
     filter = xmrdfits(filter_fil,mm+1,/silent)

     ;; Compare to variance-weighted and median values in bandpass
     gd_wv = where(wave ge min(filter.wavelength,max=mx) and $
                   wave le mx and sigflux_nu gt 0.,ngd_wv)
     case option of
        'median': begin
           f_ave = median(flux_nu[gd_wv],/even)
           sigf_ave = median(abs(flux_nu[gd_wv] - f_ave),/even)
        end
        'mean': begin
           f_ave = mean(flux_nu[gd_wv])
           sigf_ave = stddev(flux_nu[gd_wv])
        end
        'wgtmean': begin
           wgt = 1/sigflux_nu[gd_wv]^2
           f_ave = total(flux_nu[gd_wv]*wgt)/total(wgt)
           sigf_ave = sqrt(1./total(wgt))
        end
        'respt': stop,'sdss_flux2ugriz() stop: option=respt not working'
        'resbig': stop,'sdss_flux2ugriz() stop: option=resbig not working'
        'resnoa': stop,'sdss_flux2ugriz() stop: option=resnoa not working'
     endcase 

     ;; Convert to AB magnitudes
     ;; u_AB = u_SDSS - 0.04 mag [2980, 4130Ang] 
     ;; z_AB = z_SDSS + 0.02 mag [7730, 11230Ang]
     ;; g [3630, 5830Ang], r [5380, 7230Ang], i [6430, 8630Ang] match AB
     ;; ydat = -2.5*alog10(flux_nu*1e-27) - 48.60 ; mag
     ;; sigydat = sigflux_nu/flux_nu * 2.5/alog(10.)
     magstrct.(i_umag+2*mm) = -2.5*alog10(f_ave) - 48.60 ; mag
     magstrct.(i_umag+(2*mm+1)) = sigf_ave/f_ave * 2.5/alog(10.)
  endfor ; loop mm=nmag

  if not keyword_set(AB) then begin
     ;; Adjust
     ;; u_AB = u_SDSS - 0.04 mag [2980, 4130Ang] 
     ;; z_AB = z_SDSS + 0.02 mag [7730, 11230Ang]
     magstrct.umag = magstrct.umag + 0.04
     magstrct.umagerr = magstrct.umagerr + 0.04
     magstrct.zmag = magstrct.zmag - 0.02
     magstrct.zmagerr = magstrct.zmagerr - 0.02
  endif

  return,magstrct
end ; sdss_flux2ugriz()



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_functions, compile=compile
  ;; Print everything in this package
  ;; if /compile set, then run silent
  if keyword_set(compile) then return

  print,'SDSS_GETSDSSDIR() -- return concatenated $SDSSPATH/$SDSSDIR/'
  print,'SDSS_GETNAME() -- return spSpec-jjjjj-pppp-fff* name by convention.'
  print,'SDSS_GETQSOSTRCT() -- return DR X catalog structure or file name.'
  print,'SDSS_GETSNRSTRCT() -- return DR X S/N structure or file name.'
  print,'SDSS_GETQSOLIST() -- return DR X catalog spectra list or file name.'
  print,'SDSS_WRQSOLIST -- write formatted QSO list from structure.'
  print,'SDSS_GETQSOINLIST -- find and save subset of SDSS structure based on list.'
  print,'SDSS_EXCLBAL -- compare lists to exclude files and write new list.'
  print,'SDSS_GETCFLG() -- return continuum flag 0 (spline) or 1 (eigen; default).'
  print,'SDSS_GETRATING() -- return convenction for rating tag (-1, 0, 1, 2, 3).'
  print,'SDSS_SRTCIVSTRCT() -- sort elements in sdsscivstrct by QSO name and zabs.'
  print,'SDSS_RMCIVDUPLICATES() -- remove duplicate qso_name, zabs_orig[0], cflg.'
  print,'SDSS_GETRATEDDBLT() -- return sample of doublets with desired rating.'
  print,'SDSS_GETLIMFLG() -- report the analyze, upper, or lower limit flags (default: 1).'
  print,'SDSS_SETLIMFLG() -- combine upper and/or lower limit flags to input.'
  print,'SDSS_GETEWFLG() -- report the EW type flags (default: 1 or orig).'
  print,'SDSS_SETEWFLG() -- combine EW type flags to input.'
  print,'SDSS_GETBLENDFLG() -- report the blend flags (default: 1).'
  print,'SDSS_SETBLENDFLG() -- combine blend flags to input.'
  print,'SDSS_GETBALFLG() -- return the BAL flag (default: 4).'
  print,'SDSS_FNDBALQSO() -- return indices of BAL QSOs in input structure.'
  print,'SDSS_INSTANTBALFLG -- instantiate BALFLG tag in sdsscontistrct.'
  print,'SDSS_GETSPECWAVE() -- return SDSS wavelength range [3820, 9200].'
  print,'SDSS_GETSKYLINWAVE() -- return strong sky line wavelengths [5579, 6302].'
  print,'SDSS_SETCOSMOLOGY() -- instantiate XIDL cosm_common and return cosmology array [default WMAP5].'
  print,'SDSS_CALCDXDZ() -- return dX/dz given z and cosmology[3].'
  print,'SDSS_GETSPECPIXSCALE() -- return SDSS standard pixel scale [69 km/s].'
  print,'SDSS_GETRANDQSOSMPL -- select random QSOs from DR X catalog.'
  print,'SDSS_CALCNORMERR() -- calculate the total error including continuum error.'
  print,'SDSS_NORMSPEC() -- normalize SDSS spectrum and write to new standard SDSS spectrum file.'
  print,'SDSS_PLTCONTI -- plot spectrum and eigen or spline continuum.'
  print,'SDSS_PLTSNR -- plot convolved S/N from abslin structure.'
  print,'SDSS_PLTABSLIN -- plot spectrum with auto-detected lines.'
  print,'SDSS_CONTIQUAL() -- measure eigen or spline continuum fit quality.'
  print,'SDSS_MEASURESNR() -- measure S/N for input spectra and/or compute dblt bounds.'
  print,'SDSS_MKSNRTAB -- create structure of all S/N and observed wvlim per doublet.'
  print,'SDSS_PRNTSNRTAB -- print out a variety of S/N-based tables/lists.'
  print,'SDSS_REORGANIZE -- move old EIGCONTI/, SPLCONTI/, ABSLIN/ files to new directories.'
  print,'SDSS_SRTBALQSO -- divide input into BAL sightlines and else.'
  print,'SDSS_PRNTDR7SUMM -- print some doublet cuts and statistics.'
  print,'SDSS_CALCPARALLELJOB() -- divide input for given processor.'
  print,'SDSS_CATPARALLELJOB -- concatenate sdss_civsearch-in-parallel outputs.'
  print,'SDSS_HISTOGRAM() -- bin un-evenly to have roughly constant numbers.'
  print,'SDSS_BINTOLOC() -- bin to un-even location from sdss_histogram().'
  print,'SDSS_MTCHCIVSTRCT -- match elements in two sdsscivstrct files.'
  print,'SDSS_PRNTCANDSUMM -- print basic information from sdsscivstrct.'
  print,'SDSS_RDNOTEFIL -- read notes file from sdss_chkciv.'
  print,'SDSS_WRNOTEFIL -- write notes file.'
  print,'SDSS_CPSTRCT() -- copy one structure (array) to another.'
  print,'SDSS_SRTVPSTRCT() -- sort Voigt Profile structure by id_sys and id_comp.'
  print,'SDSS_CPVP2MCSTRCT() -- copy Voigt Profile structure to MC structure.'
  print,'SDSS_CPMC2CIVSTRCT() -- copy MC structure to CIV structure.'
  print,'SDSS_GETQSOSUBSET -- extract one QSO per z and S/N bin.'
  print,'SDSS_ITERQSOSUBSET -- iteratively concatenate sdss_getqsosubset results.'
  print,'SDSS_CALCFRACPIX() -- return array including fractional contribution.'
  print,'SDSS_CALCDZLOS() -- return dz and dX arrays including frac. pixels.'
  print,'SDSS_FNDQSOPAIRS() -- return structure of closest QSOs.'
  print,'SDSS_MKZARR() -- set redshift array based on z limits and binsize.'
  print,'SDSS_MKEWARR() -- set EW array based on limits and binsize.'
  print,'SDSS_CALCSIGPOISS() -- calculate Poisson error on number array.'
  print,'SDSS_CALCSIGBINOM() -- calculate Binomial confidence interval for input, recovered numbers.'
  print,'SDSS_MEDIANW() -- calculate the weighted median.'
  print,'SDSS_MINWAD() -- minimize weighted average deviation (by Leslie A. Young).'
  print,'SDSS_CALCNCMPLT() -- calculate completeness-corrected number.'
  print,'SDSS_CALCNCOLM() -- calculate completeness-corrected total column.'
  print,'SDSS_SETNCOLMFLG() -- better analysis of satured and upper limits.'
  print,'SDSS_CALCNCOLMDBLT() -- combine column measures and flags.'
  print,'SDSS_GETEWCLIM() -- return the EW or column at requested completeness.'
  print,'SDSS_CATCMPLTSTR() -- return a concatenated completeness structure.'
  print,'SDSS_GETDXW() -- query completeness array for given z, EW.'
  print,'SDSS_PRNTDNDX -- print structure of all sorts of dN/dX or dN/dz stuff.'
  print,'SDSS_PLTDNDX() -- plot dN/dX with optional extra dN/dX structure.'
  print,'SDSS_GETUSERBIASINDEX() -- match biasuser_fil redshift binning to desired redshift limits.'
  print,'SDSS_CALCDNDX() -- return structure of all sorts of dN/dX or dN/dz stuff.'
  print,'SDSS_PRNTFXW() -- print structure of all sorts of f(X,W) or f(z,W) stuff.'
  print,'SDSS_CALCFXWFIT() -- return fit of given form for data.'
  print,'SDSS_CALCNCOLMFIT() -- sample column densities from f(N) given limits.'
  print,'SDSS_CALCDNDXFIT() -- return dN/dX from fit of given form.'
  print,'SDSS_EXPGROWTHFUNC -- fitting function with correct params for A0(1-exp(A1*x)).'
  print,'SDSS_CALCEXPGROWTHFUNC() -- returns the above fit for given coefficients.'
  print,'SDSS_GETCIVSTRCT() -- return subset of given or default CIV systems.'
  print,'SDSS_ESTCZN() -- return an estimated C(N) structure given a C(W) structure and actual data.'
  print,'SDSS_PLTFXW -- plot f(W) with optional fit for given sdss_calcfxw() structure.'
  print,'SDSS_CALCFXW() -- return structure of all sorts of f(X,W) or f(z,W) stuff.'
  print,'SDSS_MKFXWSTACKS -- make a stack of each f(W) bin.'
  print,'SDSS_CALCFXWCHISQ() -- return chi^2 of fit and binned f(W).'
  print,'SDSS_PRNTOMEGA() -- print structure of all sorts of Omega stuff.'
  print,'SDSS_CALCOMEGA() -- return structure of all sorts of Omega stuff.'
  print,'SDSS_COMBINECIV() -- combine CIV detections smartly.'
  print,'SDSS_FUNSTATS -- print total mass, volume, and box size of survey.'
  print,'SDSS_GETURL() -- return array of SDSS quick view URLs.'
  print,'SDSS_APPLYCIVOBSCORR() -- modify input completeness structure for pathlength blocked by detections.'
  print,'SDSS_APPLYUSERBIAS() -- modify input completeness structure for user bias fit.'
  print,'SDSS_EXPANDCONTISTRCT() -- grow sdsscontistrct to large npix, (nlin).'
  print,'SDSS_EXPANDQALSTRCT() -- grow qalcharstrct to large nlin.'
  print,'SDSS_MKSTACKSUMM -- take in list of stacked spectra and compile info.'
  print,'SDSS_PRINTRATEDSUMM -- print summary of ratings, BAL flags, etc.'
  print,'SDSS_GETSTACKDAT() -- compile all ions at one z or one ion at all z.'
  print,'SDSS_CPZHUMENARD2STRCT() -- copy Zhu & Menard (2013) structure to sdsscivstrct.'
  print,'SDSS_CPPROCHTERETAL2STRCT() -- copy Prochter et al. (2006) table to sdsscivstrct.'
  print,'SDSS_CPQUIDERETAL2STRCT() -- copy Quider et al. (2011) table to sdsscivstrct.'
  print,'SDSS_FLUX2UGRIZ() -- estimate magnitudes from standard SDSS spectrum.'
  return
end 
