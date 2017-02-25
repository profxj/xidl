;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_find.pro                
; Author: Kathy Cooksey                      Date: 18 Feb 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Blind search (e.g. no Lya necessary) for 
;              doublets CIV and other likely lines (OVI, NV, SiIV,
;              Lya, Lyb, CIII, SiIII, NIII) in searchspec structure
; Input: 
;   srchstrct - searchspec structure (either file name or structure)
; Optional:
;   zlim - redshift limits to search
;   instrfil - instrument file name (to set upper limits to 
;              CIV 1550, Lya, CIII)
;   debug - stop occasionally
;   inspec - instead of instrfil, just input structure of one spectrum
;            with tags: wave, flux, error
;   instr_flg - binary instrument flag number for inspec
;   silent - suppress printing messages
; Output: 
;   root + '.fits' = civcandstrct (only valid
;            info is ion, wrest, wv_lim due to modifications to
;            srchstrct information)
;       [+ '_noCIV.fits' = null file b/c no CIV]
;        + '_orig.fits' = civcandstrct with original srchstrct info
;        + '_allID.fits' = civcandstrct mapping srchstrct to possible
;           multiple identifications within zlim
;        + '_allID.tab' = formatted table of above info
;   ostrct - return the final structure (or -1 if nothing found)
; Calls:
;   civcandstrct__define
;   civ_instantstrct
;   dblt_retrieve
;   klc_dblt_check
;   measureSN
; Example:
; History:
;   18 Feb 2008  created by KLC, adapted from dblt_find
;   20 Feb 2008  instantiate civcandstrct; correct storing of nwlim
;   26 Feb 2008  don't force flg_colm = 5
;   21 May 2008  correct instrument file access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@measureSN

pro civ_find_prntallid,allidstrct,outtab
  ;; Prints explicitly formatted table of pairs ID
  if size(allidstrct,/type) eq 7 then $
     allid = xmrdfits(allidstrct,1,/silent) $
  else allid = allidstrct

  ;; Print table of 
  nauto = n_elements(allid)
  dv = 2.998e5*(allid.wv_lim[0,1]-allid.wv_lim[0,0])/double(allid.ion[0])

  wrest = allid[uniq(allid.wrest,sort(allid.wrest))].ion
  if not keyword_set(outtab) then outtab = 'civ_find_prntallid.tab'
  openw,1,outtab
  printf,1,'Wobs','dv_obs','Ion', 'zabs',format='(a9,2x,a6,2x,a9,": ",a8)'
  printf,1,'(Ang)','(km/s)',format='(a9,2x,a6)'
  for ii=0,nauto-1 do begin
     printf,1,allid[ii].ion[0],round(dv[ii]),format='(a9,2x,i6,$)'
     gd = where(allid[ii].wrest gt 0,ngd)
     if ngd ne 0 then for jj=0,ngd-1 do printf,1,allid[ii].wrest[gd[jj]],$
       double(allid[ii].ion[0])/allid[ii].wrest[gd[jj]]-1.,$
       format='(2x,f9.4,": ",f8.5,$)'
     printf,1,''
  endfor                        ;loop nauto
  close,1
  print,'civ_find_allid: created ',outtab
end                             ;civ_find_prntallid


function civ_find_cp2civcand,indx,fusestrct,civstrct,$
                             create=create,calcrest=calcrest
  ;; Copy fuselinstrct to civstrct
  if keyword_set(create) then $
                   nwcivstrct = civ_instantstrct(n_elements(fusestrct)) $
  else nwcivstrct = civstrct
  nwcivstrct.ion[indx] = fusestrct.ion
  nwcivstrct.wrest[indx] = fusestrct.wrest
  nwcivstrct.wv_lim[indx,0] = fusestrct.wv_lim[0]
  nwcivstrct.wv_lim[indx,1] = fusestrct.wv_lim[1]
  nwcivstrct.flg_colm[indx] = fusestrct.flg
  nwcivstrct.instr[indx] = fusestrct.instr
  if keyword_set(calcrest) then begin
     nwcivstrct.ew[indx] = fusestrct.ew[0]/(1+(fusestrct.zabs>0.))   ;EWr
     nwcivstrct.sigew[indx] = fusestrct.sigew[0]/(1+(fusestrct.zabs>0.)) ;EWr
  endif else begin
     nwcivstrct.ew[indx] = fusestrct.ew[0]
     nwcivstrct.sigew[indx] = fusestrct.sigew[0]
  endelse 
  nwcivstrct.ncolm[indx] = fusestrct.ncolm
  nwcivstrct.signcolm[indx] = fusestrct.signcolm
  nwcivstrct.zabs[indx] = fusestrct.zabs
  nwcivstrct.zsig[indx] = fusestrct.zsig
  return,nwcivstrct
end                             ;civ_find_cp2civcand


function civ_find_cp2fuselin,indx,civstrct,fusestrct,$
                             create=create,calcrest=calcrest
  ;; Copy civcandstrct to fuselinstrct
  if keyword_set(create) then $
     nwfusestrct = replicate({ fuselinstrct },n_elements(civstrct)) $
  else nwfusestrct = fusestrct
  nwfusestrct.ion = civstrct.ion[indx]
  nwfusestrct.wrest = civstrct.wrest[indx]
  nwfusestrct.wv_lim[0] = civstrct.wv_lim[indx,0]
  nwfusestrct.wv_lim[1] = civstrct.wv_lim[indx,1]
  nwfusestrct.flg = civstrct.flg_colm[indx]
  nwfusestrct.instr = civstrct.instr[indx]
  if keyword_set(calcrest) then begin
     nwfusestrct.ew[0] = civstrct.ew[indx]/(1+(civstrct.zabs[indx]>0.))    ;EWr
     nwfusestrct.sigew[0] = civstrct.sigew[indx]/(1+(civstrct.zabs[indx]>0.))
  endif else begin
     nwfusestrct.ew[0] = civstrct.ew[indx]
     nwfusestrct.sigew[0] = civstrct.sigew[indx]
  endelse 
  nwfusestrct.ncolm = civstrct.ncolm[indx]
  nwfusestrct.signcolm = civstrct.signcolm[indx]
  nwfusestrct.zabs = civstrct.zabs[indx]
  nwfusestrct.zsig = civstrct.zsig[indx]
  return,nwfusestrct
end                             ;civ_find_cp2fuselin


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_find,srchstrct,siiv=siiv,root=root,zlim=zlim,instrfil=instrfil,$
             debug=debug, ostrct=ostrct, inspec=inspec, instr_flg=instr_flg,$
             silent=silent
  if size(srchstrct,/type) eq 7 then auto = xmrdfits(srchstrct,1,/silent) $
  else auto = srchstrct
  nauto = n_elements(auto)

  ;; Return some data if not automatically called for
  if not keyword_set(root) and not arg_present(ostrct) then root = 'civ_find'
  if not keyword_set(zlim) then zlim = [0.,10.]

  ;;Define doublets

  if keyword_set(siiv) then begin
     ;;;;;;;;;;;;;;;;;;
     ;; WARNING!!!!!
     ;;;;;;;;;;;;;;;;;;
     ;; Yes, this is going to be really confusing but it's easiest
     civ = dblt_retrieve('siiv')
     siivstrct = dblt_retrieve('civ')
     dblt_name = 'SiIV'
  endif else begin
     civ = dblt_retrieve('civ')
     siivstrct = dblt_retrieve('siiv')
     dblt_name = 'CIV'
  endelse 
  ovi = dblt_retrieve('ovi')
  nv = dblt_retrieve('nv')
  lya = dblt_retrieve('lya')    ;Lya (I) and Lyb (II)

  ;; CIV 1550
  civ1550 = dblt_retrieve('')
  civ1550.ion = civ.ion
  civ1550.wvI = civ.wvII
  civ1550.fI = civ.fII
  ;; Lyb
  lyb = dblt_retrieve('')
  lyb.ion = lya.ion
  lyb.wvI = lya.wvII
  lyb.fI = lya.fII
  ;; CIII (I) 
  ciii = dblt_retrieve('CIII')  ;CIII (I) and Lya (II)
  ;; NIII (I) 
  niii = dblt_retrieve('')
  niii.ion = 'NIII'
  niii.wvI = 989.799
  niii.fI = 0.1066
  ;; SiIII (I) 
  siiii = dblt_retrieve('')
  siiii.ion = 'SiIII'
  siiii.wvI = 1206.500
  siiii.fI = 1.66
  ;; SiII (I)
  siii = dblt_retrieve('')
  siii.ion = 'SiII'
  siii.wvI = 1260.4221
  siii.fI = 1.0070


  getfnam,[civ.wvI,civ.wvII],fval,civnam
  pair = [ovi,nv,siivstrct]
  getfnam,pair.wvI,fval,pairInam
  getfnam,pair.wvII,fval,pairIInam
  sngl = [lya,lyb,ciii,niii,siiii,siii]
  getfnam,sngl.wvI,fval,snglnam

  npair = n_elements(pair)
  nsngl = n_elements(sngl)

  ;;;;;;;;;;;;;;;;;;;
  ;; Load all CIV 1548 w/in zlim into civcandstrct
  ;;;;;;;;;;;;;;;;;;;
  zciv = auto.wrest/civ.wvI - 1.
  dzciv = (civ.wvII-civ.wvI)/civ.wvI ; leave space to cover CIV 1550
  gd = where(zciv ge zlim[0] and (zciv+dzciv) le zlim[1],nciv)
  if nciv eq 0 then begin
     if keyword_set(root) then begin
        if not keyword_set(silent) then $
           print,'civ_find: nothing found for ',root
        ofit = root+'_no'+dblt_name+'.fits'
        mwrfits,[-1],ofit,/create,/silent
        if not keyword_set(silent) then $
           print,'civ_find: created ',ofit
     endif 
     ostrct = -1
     return
  endif else begin
     if not keyword_set(silent) then begin
        print,'civ_find: found '+dblt_name+' '+strtrim(round(civ.wvI),2),nciv
     endif 
     civcand = civ_instantstrct(nciv)
     mask1550 = replicate(0,nauto)  

     ;; Store all possible identifications per automatically detected
     ;; feature 
     allid = civ_instantstrct(nauto)
     allid.ion[0] = auto.ion    ;technically wrest
     allid.wv_lim[0,0] = auto.wv_lim[0]
     allid.wv_lim[0,1] = auto.wv_lim[1]
     allid.ew[0] = auto.ew[0]
     allid.sigew[0] = auto.sigew[0]
     allid.instr[0] = auto.instr

     iciv = 0                   ;index for civcand system line
     ;; Load line information
     sauto = auto[gd]
     sauto.ion = civnam[0]
     sauto.wrest = civ.wvI
     sauto.zabs = zciv[gd]
     sauto.flg = 1              ;analyze
     civcand = civ_find_cp2civcand(iciv,sauto,civcand,/calcrest)
     civcand_nwlim = civcand    ; only wv_lim ever changed
     ;; Redshift range spanned by CIV 1548
     zcivlim = dblarr(nciv,2)
     zcivlim[*,0] = civcand.wv_lim[iciv,0]/civcand.wrest[iciv]-1.
     zcivlim[*,1] = civcand.wv_lim[iciv,1]/civcand.wrest[iciv]-1. 
     allid[gd].wrest[iciv] = civ.wvI

     iciv = iciv+1              ;increment civcand system line index
  endelse 


  ;;;;;;;;;;;;;;;;;;;
  ;; Add CIV 1550 
  ;;;;;;;;;;;;;;;;;;;
  zciv1550 = auto.wrest/civ1550.wvI - 1.
  gd = where(zciv1550 ge zlim[0] and zciv1550 le zlim[1],ngd) ;num possible
  if ngd ne 0 then begin
     allid[gd].wrest[iciv] = civ1550.wvI

     ;; Match to what's stored in civcand
     sauto = auto[gd]
     sauto.ion = civnam[1]
     sauto.zabs = zciv1550[gd]
     sauto.wrest = civ1550.wvI
     sauto.ew[0] = sauto.ew[0]/(1+(sauto.zabs > 0.))
     sauto.sigew[0] = sauto.sigew[0]/(1+(sauto.zabs > 0.))
     sauto.flg = 1
     szlim = dblarr(ngd,2)
     szlim[*,0] = sauto.wv_lim[0]/sauto.wrest-1.
     szlim[*,1] = sauto.wv_lim[1]/sauto.wrest-1.
     for kk=0,nciv-1 do begin
        ;; Singlet must be within redshift limits of CIV 
        ;; OR doublet must be within singlet bounds
        mtch = where((sauto.zabs ge zcivlim[kk,0] and $
                      sauto.zabs le zcivlim[kk,1]) or $
                     (szlim[*,0] le civcand[kk].zabs[0] and $
                      szlim[*,1] ge civcand[kk].zabs[0]),nmtch)
        if nmtch eq 0 then continue
        if nmtch gt 1 then begin
           mn = min(sauto[mtch].zabs-civcand[kk].zabs[0],imn,/absolute)
           mtch = mtch[imn]
        endif 

        ;; Load line information
        civcand[kk] = civ_find_cp2civcand(iciv,sauto[mtch],civcand[kk])

        ;; Adjust limits (only for CIV doublet)
        mn = min([civcand[kk].wv_lim[0,*]/civcand[kk].wrest[0],$
                  civcand[kk].wv_lim[iciv,*]/civcand[kk].wrest[iciv]],max=mx)
        civcand_nwlim[kk].wv_lim[0,*] = civcand[kk].wrest[0]*[mn,mx]
        civcand_nwlim[kk].wv_lim[iciv,*] = civcand[kk].wrest[iciv]*[mn,mx]
        mask1550[gd[mtch]] = mask1550[gd[mtch]]+1 ;feature used as CIV 1550
     endfor                                       ; loop nciv
     mtch = where(civcand.wv_lim[iciv,0] gt 0.,nmtch)
     if not keyword_set(silent) then begin
        print,'civ_find: found '+strtrim(round(civ.wvII),2)+' w/ '+$
              strtrim(round(civ.wvI),2),nmtch
     endif 
  endif                         ; CIV 1550 included

  ;;;;;;;;;;;;;;;;;;;
  ;; Search for CIV 1550 seperately 
  ;;;;;;;;;;;;;;;;;;;
  unmtch = where(mask1550 eq 0 and (zciv1550-dzciv) ge zlim[0] and $
                 zciv1550 le zlim[1],nunmtch) ; leave space to cover CIV 1548
  if nunmtch ne 0 then begin
     if not keyword_set(silent) then begin
        print,'civ_find: found '+dblt_name+' '+strtrim(round(civ.wvII),2)+$
              ' w/o '+strtrim(round(civ.wvI),2),nunmtch
     endif 
     civxtr = civ_instantstrct(nunmtch)
     
     ;; Load line information
     sauto = auto[unmtch]
     sauto.ion = civnam[1]
     sauto.zabs = auto[unmtch].wrest/civ.wvII - 1.
     sauto.wrest = civ.wvII
     sauto.flg = 1              ;analyize
     civxtr = civ_find_cp2civcand(iciv,sauto,civxtr,/calcrest)

     ;; Redshift range spanned by CIV 1550
     zcivlimxtr = dblarr(nunmtch,2)
     zcivlimxtr[*,0] = civxtr.wv_lim[iciv,0]/civxtr.wrest[iciv]-1.
     zcivlimxtr[*,1] = civxtr.wv_lim[iciv,1]/civxtr.wrest[iciv]-1. 

     ;; Load info for CIV 1548 limit
     civxtr.zabs[iciv-1] = civxtr.zabs[iciv] ;CIV 1548 (driver)
     civxtr.wv_lim[iciv-1,0] = civ.wvI*(1.+zcivlimxtr[*,0])
     civxtr.wv_lim[iciv-1,1] = civ.wvI*(1.+zcivlimxtr[*,1])
     civxtr_nwlim = civxtr      ;only wv_lim ever changed

     allid[unmtch].wrest[iciv] = civ.wvII ;iciv = 1

     ;; Append with enough dummy information to proceed
     civcand = [civcand,civxtr]
     civcand_nwlim = [civcand_nwlim,civxtr_nwlim]
     zcivlim = [zcivlim,zcivlimxtr]
     zciv = [zciv,civxtr.zabs[iciv]]
     nciv = nciv + nunmtch
  endif                         ;end add CIV 1550 w/o 1548
  iciv = iciv + 1               ;now can move on

  if iciv lt 2 then stop,'civ_find: indexing off for doublet'

  ;;;;;;;;;;;;;;;;;;;
  ;; Add single lines
  ;;;;;;;;;;;;;;;;;;;
  for jj=0,nsngl-1 do begin
     if not keyword_set(silent) then $
        print,'civ_find: finding '+snglnam[jj]+'...'
     
     zsngl = auto.wrest/sngl[jj].wvI - 1.
     gd = where(zsngl ge zlim[0] and zsngl le zlim[1],ngd)
;     print,'civ_find: ... found ',ngd
     if ngd eq 0 then begin
        if not keyword_set(silent) then $
           print,'civ_find:  ... none'
        continue
     endif 
     allid[gd].wrest[iciv] = sngl[jj].wvI

     ;; Match to what's stored in civcand
     sauto = auto[gd]
     sauto.ion = snglnam[jj]
     sauto.zabs = zsngl[gd]
     sauto.wrest = sngl[jj].wvI
     sauto.ew[0] = sauto.ew[0]/(1+(sauto.zabs > 0.))
     sauto.sigew[0] = sauto.sigew[0]/(1+(sauto.zabs > 0.))
     sauto.flg = 1
     szlim = dblarr(ngd,2)
     szlim[*,0] = sauto.wv_lim[0]/sauto.wrest-1.
     szlim[*,1] = sauto.wv_lim[1]/sauto.wrest-1.
     for kk=0,nciv-1 do begin
        ;; Singlet must be within redshift limits of CIV 
        ;; OR doublet must be within singlet bounds         
        mtch = where((sauto.zabs ge zcivlim[kk,0] and $
                      sauto.zabs le zcivlim[kk,1]) or $
                     (szlim[*,0] le civcand[kk].zabs[0] and $
                      szlim[*,1] ge civcand[kk].zabs[0]),nmtch)
        if nmtch eq 0 then continue
        if nmtch gt 1 then begin
           mn = min(sauto[mtch].zabs-civcand[kk].zabs[0],imn,/absolute)
           mtch = mtch[imn]
        endif 

        ;; Load line information
        civcand[kk] = civ_find_cp2civcand(iciv,sauto[mtch],civcand[kk])

        civcand_nwlim[kk].wv_lim[iciv,*] = civcand[kk].wv_lim[iciv,*]
     endfor                     ; loop nciv
     mtch = where(civcand.wv_lim[iciv,0] gt 0.,nmtch)
     if not keyword_set(silent) then $
        print,'civ_find: ... found ',nmtch
     iciv = iciv + 1
  endfor                        ; loop nsngl


  ;;;;;;;;;;;;;;;;;;;
  ;; Find all other pairs
  ;;;;;;;;;;;;;;;;;;;
  for jj=0,npair-1 do begin
     if not keyword_set(silent) then $
        print,'civ_find: finding '+pair[jj].ion+'...'

     ;; tmp and auto will have 1:1 indices
     tmp = klc_dblt_check(auto,pair[jj],zgal=zlim[0],zqso=zlim[1],debug=debug)
     gd = where(tmp.flag ne 0,ngd)
;     print,'civ_find: ... found ',ngd
     if ngd eq 0 then begin
        if not keyword_set(silent) then $
           print,'civ_find:  ... none'
        continue
     endif 
     allid[gd].wrest[iciv] = pair[jj].wvI
     allid[tmp[gd].flag].wrest[iciv+1] = pair[jj].wvII
     
     ;; Match to what's stored in civcand
     sautoI = auto[gd]
     sautoI.ion = pairInam[jj]
     sautoI.zabs = sautoI.wrest/pair[jj].wvI - 1.
     sautoI.wrest = pair[jj].wvI
     sautoI.ew[0] = sautoI.ew[0]/(1+(sautoI.zabs > 0.))
     sautoI.sigew[0] = sautoI.sigew[0]/(1+(sautoI.zabs > 0.))
     sautoI.flg = 1
     
     sautoII = auto[tmp[gd].flag]
     sautoII.ion = pairIInam[jj]
     sautoII.zabs = sautoII.wrest/pair[jj].wvII - 1.
     sautoII.wrest = pair[jj].wvII
     sautoII.ew[0] = sautoII.ew[0]/(1+(sautoII.zabs > 0.))
     sautoII.sigew[0] = sautoII.sigew[0]/(1+(sautoII.zabs > 0.))
     sautoII.flg = 1

     ;; Significance-weighted redshift for comparison
     nsigI = abs(sautoI.ew[0]/sautoI.sigew[0])
     nsigII = abs(sautoII.ew[0]/sautoII.sigew[0])
     zpair = (sautoI.zabs*nsigI + sautoII.zabs*nsigII)/(nsigI+nsigII)
     for kk=0,nciv-1 do begin
        ;; Incorporate pair[jj] if EW-significance-weighted average
        ;; redshift withing CIV 1548 bound 
        mtch = where(zpair ge zcivlim[kk,0] and zpair le zcivlim[kk,1],nmtch)
        if nmtch eq 0 then continue
        if nmtch gt 1 then begin
           mn = min(zpair[mtch]-civcand[kk].zabs[0],imn,/absolute)
           mtch = mtch[imn]
        endif

        ;; Load line information
        civcand[kk] = civ_find_cp2civcand(iciv,sautoI[mtch],civcand[kk])
        civcand[kk] = civ_find_cp2civcand(iciv+1,sautoII[mtch],civcand[kk])

        ;; Adjust doublet wv_lim
        mn = min([civcand[kk].wv_lim[iciv,*]/civcand[kk].wrest[iciv],$
                  civcand[kk].wv_lim[iciv+1,*]/civcand[kk].wrest[iciv+1]],$
                 max=mx)
        civcand_nwlim[kk].wv_lim[iciv,*] = civcand[kk].wrest[iciv]*[mn,mx]
        civcand_nwlim[kk].wv_lim[iciv+1,*] = civcand[kk].wrest[iciv+1]*[mn,mx]

     endfor                     ; loop nciv
     mtch = where(civcand.wv_lim[iciv,0] gt 0.,nmtch)
     if not keyword_set(silent) then $
        print,'civ_find: ... found ',nmtch
     iciv = iciv + 2
  endfor                        ; loop over pairs


  ;;;;;;;;;;;;;;;;;;;
  ;; Include upper limits for CIV 1548, 1550, CIII, Lya
  ;;;;;;;;;;;;;;;;;;;
  if keyword_set(instrfil) or keyword_set(inspec) then begin
     if keyword_set(instrfil) then begin
        civcand.instr_fil = strmid(instrfil,strpos(instrfil,'MLSS/')+5)
        
        ;; Store spectra information
        readcol,instrfil,spec,dw,w0,format='a,f,f',/silent
        spec = getenv('MLSS_DIR')+'/'+spec
        ninstr = n_elements(spec)
        wvlim = fltarr(ninstr,2)
        for ii=0,ninstr-1 do begin
           test = file_search(spec[ii],count=ntest)
           if ntest eq 0 then begin
              wvlim[ii,*] = [-1,-1]
              continue
           endif 

           if ii le 6 then $    ;FUSE
              fx = x_readspec(spec[ii],wav=wv,inflg=3) $
           else fx = x_readspec(spec[ii],wav=wv) ;stis
           gd = where(fx ne 0,ngd)
           wv = wv + dw[ii]*w0[ii]/wv
           
           wvlim[ii,*] = [min(wv[gd],max=mx),mx]
        endfor
     endif else begin
        lgv = alog(instr_flg)/alog(2)
        ninstr = fix(lgv+0.000001) 
        wvlim = replicate(-1.,ninstr+1,2)
        wvlim[ninstr,*] = [min(inspec.wave,max=mx),mx]
     endelse 

     ;; Loop CIV and see if other lines there
     ;; Only include even if full feature can NOT be included
     ii = iciv
     for kk=0,nciv-1 do begin
        ;; CIV 1548  (assume zabs[1] info from CIV 1550) 
        ;; load in 0th location
        mtch = where(abs(civcand[kk].wrest-civ.wvI) lt 1e-4,nmtch)
        if nmtch eq 0 then begin
           gd = where(wvlim[*,0] le civ.wvI*(1+civcand[kk].zabs[1]) and $
                      wvlim[*,1] ge civ.wvI*(1+civcand[kk].zabs[1]),ngd)
           if ngd ne 0 then begin
              civcand[kk].ion[0] = civnam[0] ;CIV 1548
              civcand[kk].wrest[0] = civ.wvI
              civcand[kk].zabs[0] = civcand[kk].zabs[1]
              civcand[kk].wv_lim[0,*] = [civ.wvI*(1+zcivlim[kk,0]),$
                                         civ.wvI*(1+zcivlim[kk,1])]
              civcand[kk].instr[0] = total(2^gd)
              civcand[kk].flg_colm[0] = 1
              civcand_nwlim[kk].wv_lim[0,*] = civcand[kk].wv_lim[0,*]
           endif                ;ngd ne 0
        endif                   ; CIV 1548

        ;; CIV 1550 
        ;; load 1st location
        mtch = where(abs(civcand[kk].wrest-civ1550.wvI) lt 1e-4,nmtch)
        if nmtch eq 0 then begin
           gd = where(wvlim[*,0] le civ1550.wvI*(1+civcand[kk].zabs[0]) and $
                      wvlim[*,1] ge civ1550.wvI*(1+civcand[kk].zabs[0]),ngd)
           if ngd ne 0 then begin
              civcand[kk].ion[1] = civnam[1] ;CIV 1550
              civcand[kk].wrest[1] = civ1550.wvI
              civcand[kk].zabs[1] = civcand[kk].zabs[0]
              civcand[kk].wv_lim[1,*] = [civ1550.wvI*(1+zcivlim[kk,0]),$
                                         civ1550.wvI*(1+zcivlim[kk,1])]
              civcand[kk].instr[1] = total(2^gd)
              civcand[kk].flg_colm[1] = 1
              civcand_nwlim[kk].wv_lim[1,*] = civcand[kk].wv_lim[1,*]
           endif                ;ngd ne 0
        endif                   ; CIV 1550

        ;; CIII
        mtch = where(abs(civcand[kk].wrest-ciii.wvI) lt 1e-4,nmtch)
        if nmtch eq 0 then begin
           gd = where(wvlim[*,0] le ciii.wvI*(1+civcand[kk].zabs[0]) and $
                      wvlim[*,1] ge ciii.wvI*(1+civcand[kk].zabs[0]),ngd)
           if ngd ne 0 then begin
              getfnam,ciii.wvI,fval,nam
              civcand[kk].ion[ii] = nam
              civcand[kk].wrest[ii] = ciii.wvI
              civcand[kk].zabs[ii] = civcand[kk].zabs[0]
              civcand[kk].wv_lim[ii,*] = [ciii.wvI*(1+zcivlim[kk,0]),$
                                          ciii.wvI*(1+zcivlim[kk,1])]
              civcand[kk].instr[ii] = total(2^gd) ;could isolate
              civcand[kk].flg_colm[ii] = 1
              civcand_nwlim[kk].wv_lim[ii,*] = civcand[kk].wv_lim[ii,*]
              ii = ii + 1
           endif                ;ngd ne 0
        endif                   ; CIII

        ;; Lya
        mtch = where(abs(civcand[kk].wrest-lya.wvI) lt 1e-4,nmtch)
        if nmtch eq 0 then begin
           gd = where(wvlim[*,0] le lya.wvI*(1+civcand[kk].zabs[0]) and $
                      wvlim[*,1] ge lya.wvI*(1+civcand[kk].zabs[0]),ngd)
           if ngd ne 0 then begin
              getfnam,lya.wvI,fval,nam
              civcand[kk].ion[ii] = nam
              civcand[kk].wrest[ii] = lya.wvI
              civcand[kk].zabs[ii] = civcand[kk].zabs[0]
              civcand[kk].wv_lim[ii,*] = [lya.wvI*(1+zcivlim[kk,0]),$
                                          lya.wvI*(1+zcivlim[kk,1])]
              civcand[kk].instr[ii] = total(2^gd) ;could isolate...
              civcand[kk].flg_colm[ii] = 1
              civcand_nwlim[kk].wv_lim[ii,*] = civcand[kk].wv_lim[ii,*]
              ii = ii + 1
           endif                ;ngd ne 0
        endif                   ; Lya
        ii = iciv               ;reset
     endfor                     ; loop nciv
  endif                         ;keyword instrfil set

  ;;;;;;;;;;;;;;;;;;;
  ;; Check 
  ;; (everything must have at least CIV 1548, 1550)
  ;;;;;;;;;;;;;;;;;;;
  gd = where(abs(civcand.wrest[0]-civ.wvI) lt 1e-4 and $
             abs(civcand.wrest[1]-civ.wvII) lt 1e-4,ngd)
  if ngd eq 0 then begin
     ;; Will be freaky for this to happen, but I guess GHRS ECH could
     ;; just have one feature and no room for a limit...
     if not keyword_set(silent) then $
        print,'civ_find: ultimately no '+dblt_name+' in ',root
     if keyword_set(root) then begin
        ofit = root+'_no'+dblt_name+'.fits'
        mwrfits,[-1],ofit,/create,/silent
        if not keyword_set(silent) then $
           print,'civ_find: created ',ofit 
     endif 
     ostrct = -1
     return
  endif else begin
     ;; Really just excluding edge cases (no wavelength coverage)
     civcand = civcand[gd]
     civcand_nwlim = civcand_nwlim[gd]
  endelse  

  ;;;;;;;;;;;;;;;;;;;
  ;; Output information
  ;;;;;;;;;;;;;;;;;;;
  srt = sort(civcand.zabs[0])   ;sort by redshift
  civcand = civcand[srt]
  civcand_nwlim = civcand_nwlim[srt]

  if size(srchstrct,/type) eq 7 then civcand.search_fil = srchstrct

  if keyword_set(root) then begin
     ofit = root+'_allID.fits'
     mwrfits,allid,ofit,/create,/silent
     if not keyword_set(silent) then $
        print,'civ_find: created ',ofit
     civ_find_prntallid,allid,root+'_allID.tab'
;  save,allid,filename=root+'_allID.idl'
     
     ofit = root+'_orig.fits'
     mwrfits,civcand,ofit,/create,/silent
     if not keyword_set(silent) then $
        print,'civ_find: created ',ofit
;  save,civcand,filename=root+'_orig.idl'
  endif 

  ;; Adjust wavelength limits (and erase no-longer pertinant data)
  civcand.wv_lim = civcand_nwlim.wv_lim
  civcand.ew = 0.
  civcand.sigew = 0.
  civcand.zabs = 0.
  civcand.zsig = 0.

  ;; Write final files
  if keyword_set(root) then begin
     ofit = root+'.fits'
     mwrfits,civcand,ofit,/create,/silent 
     if not keyword_set(silent) then $
        print,'civ_find: created ',ofit
;  save,civcand,filename=root+'.idl'
  endif
  ostrct = civcand

  if keyword_set(debug) then stop,'civ_find debug: finished'

end
