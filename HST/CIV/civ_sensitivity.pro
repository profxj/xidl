;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_sensitivity.pro
; Author: Kathy Cooksey     Date: 27 Feb 2008
; Project: HST Metal-line System Survey (CIV focus)
; Description: Measure either single or doublet sensitivity 
;              function of spectra
; Inputs:
;   strct_fil - civcandstrct format file (need instr_fil)
;   wrest - scalar or 2-element array of rest wavelengths
;   dopb - (km/s) defines width of line (= b/c*wrest Ang)
;   ewlim - (mA) minimum EW of feature, determines sensitivity curve
;   nsigma - desired significance of detection
; Optional:
;   outroot - root of filname for out put (otherwise use base of 
;             instrument name)
; Output: 
;   [root]_gzIONbkmsEWmA.fits - structure saved
; Optional Output: 
; History:
;   27 Feb 2008  created by KLC
;    3 Mar 2008  super effor to make fast, compact, right  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@spec_ewlin

function civ_sensitivity_fndcolm,wrest,dopb,ewlim,tol=tol
;; Converge on logN based on EW and b (within tolerance)
if not keyword_set(tol) then tol = 0.1 ;mA


logn = 13.5
ew = spec_ewlin(wrest,nb=[logn,dopb],/silent)
if abs(ewlim-ew) le tol then done = 1 else done = 0
cnt = 0
dlogn = 0.5
lst = 0                         ;last change
while not done do begin
    if ew ge ewlim then begin
        ;; Decrease logN
        if lst gt 0 then dlogn = 0.5*dlogn
        logn = logn - dlogn
        ew = spec_ewlin(wrest,nb=[logn,dopb],/silent)
        lst = -1
    endif else begin
        ;; Increase logN 
        if lst lt 0 then dlogn = 0.5*dlogn
        logn = logn + dlogn
        ew = spec_ewlin(wrest,nb=[logn,dopb],/silent)
        lst = 1
    endelse 
    if abs(ew-ewlim) le tol then done = 1 else cnt = cnt + 1
;    print,cnt,dlogn,logn,ew
endwhile 

return,logn
end                             ;civ_sensitivity_fndcolm


function civ_sensitivity_dvqso,strct_fil,dvqso=dvqso,zlim=zlim,zpath=zpath,$
                               dvgal=dvgal,debug=debug
  if size(strct_fil,/type) eq 7 then strct = xmrdfits(strct_fil,1,/silent) $
  else strct = strct_fil
  ngz = n_elements(strct.gz[0,*])
  gz = strct.gz[1,*]
  zabs = strct.gz[0,*] 
  
  ;; Adjust upper bound (also can just account for zqso)
  civcand = xmrdfits(strtrim(strct.civcand_fil,2),1,/silent)

  ;;;;;;;
  ;; ZQSO
  ;;;;;;;
  if not keyword_set(dvqso) then dvqso = 1
  if dvqso eq 1 then zqso = civcand[0].zqso $
  else zqso = civcand[0].zqso-dvqso/2.998e5*(1+civcand[0].zqso)

  if strct.gz[0,ngz-1] lt zqso then begin
     ;; All good, let it stand
     ;; dz remains the same
  endif else begin
     if strct.gz[0,0] gt zqso then begin
        ;; All bad
        gz = [0,0]
        zabs = [strct.gz[0,0],strct.gz[0,ngz-1]]
     endif else begin
        ;; Difficult, must divide
        mx = min(zqso-strct.gz[0,*],imx,/absolute)
        case imx of
           0: begin
              gz = [strct.gz[1,imx],strct.gz[1,imx],0,0]
              zabs = [strct.gz[0,imx],zqso,zqso,strct.gz[0,ngz-1]]
           end 
           else: begin
              if (imx mod 2) eq 0 then imx = imx + 1 ;use upper side of bin (odd)
              nnwgz = imx + 3
              gz = intarr(nnwgz)
              zabs = dblarr(nnwgz)
              ;; Preserve original
              for ii=0,imx-1 do begin
                 gz[ii] = strct.gz[1,ii]
                 zabs[ii] = strct.gz[0,ii]
              endfor 
              ;; Alter ending
              gz[imx:imx+2] = [strct.gz[1,imx-1],0,0]
              zabs[imx:imx+2]= [zqso,zqso,strct.gz[0,ngz-1]]
           end 
        endcase                 ;imx
        if keyword_set(debug) then $
           stop,'civ_sensitivity_dvqso debug: difficult case for dvqso'
     endelse                    ;difficult case
  endelse                       ;must alter g(z)
  nwgz = fltarr(2,n_elements(gz))
  nwgz[0,*] = zabs
  nwgz[1,*] = gz
  ngz = n_elements(nwgz[0,*])
  if keyword_set(debug) then stop,'civ_sensitivity_dvqso: finished dvqso'

  ;;;;;;;
  ;; ZGAL
  ;;;;;;;
  if not keyword_set(dvgal) then dvgal = 1
  if dvgal eq 1 then zgal = 0. $
  else zgal = dvgal/2.998e5 
  if nwgz[0,0] gt zgal then begin
     ;; All good, let it stand
     ;; dz remains the same
  endif else begin
     if nwgz[0,ngz-1] lt zgal then begin
        ;; All bad
        gz = [0,0]
        zabs = [nwgz[0,0],nwgz[0,ngz-1]]
     endif else begin
        ;; Difficult, must divide
        mn = min(zgal-nwgz[0,*],imn,/absolute)
        case imn of
           0: begin
              gz = [0,0,nwgz[1,imn],nwgz[1,imn]]
              zabs = [nwgz[0,imn],zgal,zgal,nwgz[0,ngz-1]]
           end
           else: begin
              if (imn mod 2) ne 0 then imn = imn - 1 ; use lower side of bin (even)
              nnwgz = ngz - imn + 2
              gz = intarr(nnwgz)
              zabs = dblarr(nnwgz)
              ;; Alter beginning
              gz[0:3] = [0,0,nwgz[1,imn],nwgz[1,imn+1]]
              zabs[0:3]= [nwgz[0,0],zgal,zgal,nwgz[0,imn+1]]
              ;; Preserve remainder
              jj = 4
              for ii=imn+2,ngz-1 do begin
                 gz[jj] = nwgz[1,ii]
                 zabs[jj] = nwgz[0,ii]
                 jj = jj+1
              endfor 
           end
        endcase                 ;imn
        if keyword_set(debug) then $
           stop,'civ_sensitivity_dvqso debug: difficult case for dvgal'
     endelse                    ;difficult case
  endelse                       ;must alter g(z)

  dz = zabs - shift(zabs,1)
  dz[0] = dz[1]
  gd = where(gz gt 0,ngd)
  if ngd ne 0 then $
     zlim = [min(zabs[gd],max=mx),mx] $
  else zlim = [-1,-1]
  zpath = total(gz*dz)

  nwgz = fltarr(2,n_elements(gz))
  nwgz[0,*] = zabs
  nwgz[1,*] = gz
  if keyword_set(debug) then stop,'civ_sensitivity_dvqso: about to return'
return,nwgz
end


pro civ_sensitivity_plot,strct_fil,cumulative=cumulative,zbin=zbin,$
                         single=single,dvqso=dvqso,dvgal=dvgal,$
                         cumfil=cumfil,debug=debug
  angstrom = STRING("305B)   
  nstrct = n_elements(strct_fil)

;;;;;;;;;;;;;;;;;;;;;
;; Cumulative g(z)
;;;;;;;;;;;;;;;;;;;;;
  if keyword_set(cumulative) then begin
     if size(cumulative,/type) ne 7 then $
        stop,'civ_sensitivity_plot: cumulative must be name of postscript'

     ;; Set cumulative info
     if not keyword_set(zbin) then zbin = 1000./2.998e5
     zlim = [0.,10.]            ;sufficiently large
     ncum = ceil(zlim[1]/zbin)
     cumz = dindgen(ncum)*zbin
     cumgz = dblarr(ncum)
     zpathtot = 0.

     if keyword_set(cumfil) then dat = xmrdfits(cumfil,1,/silent) $
     else begin
        for ii=0L,nstrct-1 do begin
           ;; Concatenate information
           strct = xmrdfits(strct_fil[ii],1,/silent)
           npix = n_elements(strct.gz[0,*])
           gz = strct.gz

           if keyword_set(dvqso) or keyword_set(dvgal) then begin
              nwgz = civ_sensitivity_dvqso(strct,dvqso=dvqso,dvgal=dvgal,$
                                           zlim=zlim,zpath=zpath,debug=debug)
              if zpath lt 0. then $
                 stop,'civ_sensitivity_plot: zpath lt 0'
              gz = nwgz
              strct.zlim = zlim
              strct.zpath = zpath
              if strct.zpath ne zpath then stop
           endif                ;dvqso set
           dz = gz[0,*] - shift(gz[0,*],1)
           dz[0] = 0.
           zpathtot = zpathtot + strct.zpath

           ;; Should be no double counting of a quasar spectrum 
           ;; b/c everything in largest bins possible
           ;; mtch is lower bownd of positive g(z)
           mtch = where(gz[1,*] eq 1 and shift(gz[1,*],-1) eq 1 and $
                        gz[0,*] ge 0.,nmtch)
           if nmtch eq 0 then $
              print,'civ_sensitivity_plot: no contribution ',strct_fil[ii] $
           else begin
              for jj=0,nmtch-1 do begin
                 ilo = floor(gz[0,mtch[jj]]/zbin) ; lower bound, left side
                 ihi = floor(gz[0,mtch[jj]+1]/zbin) ; upper bound, left side
                 ;; in same bin, add full value
                 if ihi eq ilo then $
                    cumgz[ilo] = cumgz[ilo] + $
                                 (gz[0,mtch[jj]+1]-gz[0,mtch[jj]])/zbin $
                 else begin
                    ;; Individual g(z) bin shared across bins
                    cumgz[ilo] = cumgz[ilo] + (cumz[ilo+1]-gz[0,mtch[jj]])/zbin
                    cumgz[ihi] = cumgz[ihi] + $
                                 (gz[0,mtch[jj]+1]-cumz[ihi])/zbin
                    if ihi-ilo gt 1 then for kk=ilo+1,ihi-1 do $
                       cumgz[kk] = cumgz[kk] + 1
                 endelse        ;account for consectuve bins
              endfor            ;loop nmtch
           endelse              ;LOS contributes
;           print,ii,zpathtot,total(cumgz*zbin)
       endfor                  ;loop nstrct

        ;; Trim
        if keyword_set(debug) then $
           stop,'civ_sensitivity_plot debug: about to trim cum g(z)'
        gd = where(cumgz gt 0,ngd)
        ncum = max(gd) + 1
        cumgz = cumgz[0:ncum-1]
        cumz = cumz[0:ncum-1]

        ;; Test
        if abs(zpathtot-total(cumz*cumgz))/zpathtot gt 0.05 then $
          stop,'civ_sensitivity_plot: 5% difference in total zpath measured'

        ;; Label info
        getfnam,strct.wrest[0],fval,nam
        nam = strsplit(strtrim(nam,2),/extract)

        if min(strct.ncolm) le 0 then begin
           strct.ncolm[0] = civ_sensitivity_fndcolm(strct.wrest[0],$
                                                    strct.dopb,strct.ewlim)
           if n_elements(strct.ncolm) gt 1 then $
              strct.ncolm[1] = civ_sensitivity_fndcolm(strct.wrest[1],$
                                                       strct.dopb,strct.ewlim)
        endif                   ;calc ncolm

        ;; Save info
        dat = {ion:nam[0],cumz:cumz,cumgz:cumgz,zpath:total(cumgz*zbin),$
               dopb:strct.dopb,ewlim:strct.ewlim,ncolm:min(strct.ncolm)}
        ofit = strmid(cumulative,0,strpos(cumulative,'.',/reverse_search))+$
               '.fits'
        mwrfits,dat,ofit,/create,/silent
        print,'civ_sensitivity_plot: created ',ofit
     endelse                    ;read in info

     ;; Plot
     x_psopen,cumulative,/maxs
     clr = getcolor(/load)
     !x.margin = [7,1.5]
     !y.margin = [4,1]
     xrng = [min(dat.cumz,max=mx),mx]
     yrng = [-0.1,max(dat.cumgz)*1.1] ; plus 10%
     
     zstr = 'z!D'+dat.ion+'!N'
     plot,xrng,yrng,/device,/nodata,background=clr.white,color=clr.black,$
          /ystyle,/xstyle,xtitle=zstr,ytitle='g('+zstr+')',charsize=2
     oplot,dat.cumz,dat.cumgz,psym=10,thick=6,color=clr.black

     ;; Label (upper right)
     xyouts,xrng[1]-0.3*(xrng[1]-xrng[0]),0.93*(yrng[1]-yrng[0])+yrng[0],$
            '!9D!3'+zstr+' = '+strtrim(string(dat.zpath,format='(f6.2)'),2),$
            charsize=1.7
     xyouts,xrng[1]-0.3*(xrng[1]-xrng[0]),0.88*(yrng[1]-yrng[0])+yrng[0],$
            'b = '+strtrim(round(dat.dopb),2)+' km s!E-1!N',$
            charsize=1.7
     xyouts,xrng[1]-0.3*(xrng[1]-xrng[0]),0.83*(yrng[1]-yrng[0])+yrng[0],$
            'W!Dr,lim!N = '+strtrim(round(dat.ewlim),2)+' m'+angstrom,$
            charsize=1.7
     xyouts,xrng[1]-0.3*(xrng[1]-xrng[0]),0.78*(yrng[1]-yrng[0])+yrng[0],$
            'log N!DAOD,lim!N = '+string(dat.ncolm,format='(f5.2)'),$
            charsize=1.7
     x_psclose
     print,'civ_sensitivity_plot: created ',cumulative
     if keyword_set(debug) then $
        stop,'civ_sensitivity_plot debug: about to exit cumulative'
  endif                         ;cumulative


  ;;;;;;;;;;;;;;;;;;;;;
  ;; Single EW sensitivity
  ;;;;;;;;;;;;;;;;;;;;;
  if keyword_set(single) then begin
     if size(single,/type) ne 7 then $
        psfil = strmid(strct_fil,0,strpos(strct_fil,'.',/reverse_search))+$
                '.ps' $
     else psfil = single

     if size(strct_fil,/type) eq 7 then strct = xmrdfits(strct_fil,1,/silent) $
     else strct = xmrdfits(strct_fil,1,/silent)

     zabs = strct.wave/strct.wrest[0] - 1.
     nlin = n_elements(strct.wrest)

     x_psopen,psfil,/maxs
     clr = getcolor(/load)
     !x.margin = [7,1.5]
     !y.margin = [4,1]
     xrng = [min(strct.wave,max=mx)/max(strct.wrest,min=mn)-1.,mx/mn-1.]
     yrng = [-0.1,strct.ewlim*1.2] ;plus 20%
     
     getfnam,strct.wrest[0],fval,nam
     nam = strsplit(strtrim(nam,2),/extract)
     zstr = 'z!D'+nam[0]+'!N'

     ;; Data
     plot,xrng,yrng,/device,/nodata,background=clr.white,color=clr.black,$
          /ystyle,/xstyle,xtitle=zstr,ytitle='W!Dr,lim!N (m'+angstrom+')',$
          charsize=2
     bd = where(strct.sigew[0,*] le 0,nbd)    ;gaps
     if nbd ne 0 then strct.sigew[0,bd] = 1e6 ;very large
     oplot,strct.wave/strct.wrest[0]-1,strct.nsigma*strct.sigew[0,*],$
           psym=10,thick=3,color=clr.black
     if nlin ne 1 then begin
        bd = where(strct.sigew[1,*] le 0,nbd)    ;gaps
        if nbd ne 0 then strct.sigew[1,bd] = 1e6 ;very large
        oplot,strct.wave/strct.wrest[1]-1,strct.nsigma*strct.sigew[1,*],$
              psym=10,thick=2,color=clr.red
     endif 
     oplot,[-100.,100.],[strct.ewlim,strct.ewlim],linestyle=2,$
           color=clr.black,thick=3

     dwv = strct.wave-shift(strct.wave,1)
     dwv[0] = dwv[1]
     civcand = xmrdfits(strtrim(strct.civcand_fil,2),1,/silent)
     if keyword_set(dvqso) then begin
        nwgz = civ_sensitivity_dvqso(strct,dvqso=dvqso,dvgal=dvgal,zlim=zlim,$
                                     zpath=zpath,debug=debug)
        strct.zlim = zlim
        strct.zpath = zpath
     endif                      ;dvqso set

     ;; Labels
     xyouts,0.4*(xrng[1]-xrng[0])+xrng[0],0.93*(yrng[1]-yrng[0])+yrng[0],$
            strtrim(civcand.qso,2)+$
            ' (z!DQSO!N = '+strtrim(string(civcand.zqso,$
                                           format='(f9.5)'),2)+')',$
            charsize=1.7
     xyouts,0.4*(xrng[1]-xrng[0])+xrng[0],0.88*(yrng[1]-yrng[0])+yrng[0],$
            '!9D!3'+zstr+' = '+strtrim(string(strct.zpath,format='(f6.3)'),2)+$
            ', b = '+strtrim(round(strct.dopb),2)+' km s!E-1!N',$
            charsize=1.7
     x_psclose
     print,'civ_sensitivity_plot: created ',psfil
  endif                         ;single

end                             ;civ_sensitivity_plot



pro civ_sensitivity_print,strct_fil,outtab,dvqso=dvqso,dvgal=dvgal,filnam=filnam,_EXTRA=extra
  ;; Print table for given dobp, ewlim, and nsigma for all LOS
  nstrct = n_elements(strct_fil)
  if nstrct le 1 or size(strct_fil,/type) ne 7 then $
     stop,'civ_sensitivity_print: strct_fil must be array of structure names'

  openw,1,outtab
  for ii=0,nstrct-1 do begin
     strct = xmrdfits(strct_fil[ii],1,/silent)
     if ii eq 0 then begin
         ;; Table header
         dopb = strct.dopb
         printf,1,strct.dopb,format='("# b = ",i3," km/s")'
         printf,1,strct.ewlim,format='("# EWlim = ",i3," mA")'
         printf,1,strct.nsigma,format='("# nsigma = ",i3)'
         printf,1,'QSO','zqso','RA','Dec','zmin','zmax','zpath',$
           format='(a20,1x,a9,1x,a12,1x,a12,2x,a9,1x,a9,1x,a9)'
     endif else if strct.dopb ne dopb then $
       stop,'civ_sensitivity_print: dopb changes'

     ;; Need QSO info
     civcand = xmrdfits(strtrim(strct.civcand_fil,2),1,/silent)
     unq = uniq(strtrim(civcand.qso,2),sort(civcand.qso))
     if n_elements(unq) ne 1 then $
        stop,'civ_sensitivity_print: confused civcandstrct with more than one QSO'
     civcand = civcand[0]
     x_radec,raqso,decqso,civcand.ra,civcand.dec,/flip

     if keyword_set(dvqso) then begin
        ;; Re-evaluate bounds (also sneaky way of taking into account
         ;; zqso)
        nwgz = civ_sensitivity_dvqso(strct,dvqso=dvqso,dvgal=dvgal,$
                                     zlim=zlim,zpath=zpath,_extra=extra)
        strct.zlim = zlim
        strct.zpath = zpath
     endif

     ;; Print line
     if keyword_set(filnam) then begin
        ;; Instead of QSO use strct_fil, which has info about QSO and grating
        prs = strsplit(strct_fil[ii],'/',count=nprs)
        str = strmid(strct_fil[ii],prs[nprs-1])
        prs = strsplit(str,'_',count=nprs)
        fil = strmid(str,0,prs[nprs-1]-1)
        if strct.zpath gt 0. then $
           printf,1,strtrim(fil,2),civcand.zqso,raqso,decqso,$
                  strct.zlim[0],strct.zlim[1],strct.zpath,$
       format='(a20,1x,f9.5,1x,a12,1x,a12,2x,f9.5,1x,f9.5,1x,f9.5)' $
        else printf,1,strtrim(fil,2),civcand.zqso,raqso,decqso,$
                    '...','...',0.,$
                    format='(a20,1x,f9.5,1x,a12,1x,a12,2x,a9,1x,a9,1x,f9.5)' 
     endif else begin
        if strct.zpath gt 0. then $
           printf,1,strtrim(civcand.qso,2),civcand.zqso,raqso,decqso,$
                  strct.zlim[0],strct.zlim[1],strct.zpath,$
       format='(a20,1x,f9.5,1x,a12,1x,a12,2x,f9.5,1x,f9.5,1x,f9.5)' $
        else printf,1,strtrim(civcand.qso,2),civcand.zqso,raqso,decqso,$
                    '...','...',0.,$
                    format='(a20,1x,f9.5,1x,a12,1x,a12,2x,a9,1x,a9,1x,f9.5)' 
     endelse 
     if ii gt 0 then zpathtot = zpathtot + strct.zpath $
     else zpathtot = strct.zpath 
  endfor                        ;loop nstrct
  printf,1,'Total unblocked redshift pathlength: ',$
    string(zpathtot,format='(f9.5)')

  close,1
  print,'civ_sensitivity_print: created ',outtab
end                             ;civ_sensitivity_concat



function civ_sensitivity_gz,wave,sigew,wrest,nsigma,ewlim,$
                            zpath=zpath,zlim=zlim
;; Evaluate sigew (binned error)
  npix = n_elements(wave)
  nlin = n_elements(wrest)

  gz = replicate(0L,npix)
  zabs = wave/wrest[0]-1. 
  dwv = wave-shift(wave,1)
  dwv[0] = dwv[1]
  dz = (shift(wave,-1)-wave)/wrest[0] ;assumption!!
  dz[npix-1] = dz[npix-2]

  if nlin eq 2 then begin
     ;; Doublet
     zabsII = wave/wrest[1]-1. 
     ;; See where doublet can be detected at >= 3*sigew and where
     ;; partner can also be detected (first must be in spectrum)
     gd = where(sigew[0,*] gt 0. and $ ;gaps
                sigew[0,*]*nsigma le ewlim and $
                zabsII ge min(zabs) and zabsII le max(zabs) and $
                zabs ge 0.,ngd)
     if ngd ne 0 then begin
        gz[gd] = 1              ;temporary yes
        ;; See if partner also detectable
        for ii=0L,ngd-1 do begin
           mn = min(zabs[gd[ii]]-zabsII,imn,/absolute)
           if sigew[1,imn]*nsigma gt ewlim or $
              sigew[1,imn] le 0. then gz[gd[ii]] = 0 ;no
        endfor                                       ;loop ngd
     endif 
  endif else begin
     ;; Single
     gd = where(sigew[0,*] gt 0. and $
                sigew[0,*]*nsigma le ewlim $
                and zabs ge 0.,ngd)
     if ngd ne 0 then gz[ii,gd] = 1
  endelse 

;; Trim g(z) to basic necessity
;; Discontinuities
  jmpu=where((gz-shift(gz,1)) ne 0,nu)
  jmpd=where((gz-shift(gz,-1)) ne 0,nd)
  if nu eq nd and nu ne 0 then begin
     jmp=[0,jmpu,jmpd,npix-1]
     jmp = jmp[sort(jmp)]       ;weave
     njmp = n_elements(jmp)
     gzcrs= gz[jmp]        
     zabscrs = zabs[jmp]
  endif else begin
     gzcrs = [gz[0],gz[npix-1]]
     zabscrs = [zabs[0],zabs[npix-1]]
  endelse 

;; Redshift bounds and unblocked redshift pathlength
  dzcrs = zabscrs - shift(zabscrs,1)
  dzcrs[0] = 0.
  gd = where(zabscrs ge 0 and gzcrs eq 1,ngd)
  if ngd ne 0 then $
     zlim = [min(zabscrs[gd],max=mx),mx] $
  else zlim = [-1.,-1.]
  zpath = total(dzcrs*gzcrs)
  if zlim[0] lt 0. then zlim[0] = 0. ;can occasionally be -0.
  
;; TESTS
  if abs(max(zabs*gz)-zlim[1]) gt max(dz) and zlim[1] ne -1 then $
     stop,'civ_sensitivity_gz: max redshift limit not right'
  if abs(zpath-total(dz*gz)) gt max(dz) then $
     stop,'civ_sensitivity_gz: unblocked pathlength not right'

;; Return array
  ngz = n_elements(gzcrs)
  gz = fltarr(2,ngz)
  gz[0,*] = zabscrs
  gz[1,*] = gzcrs
;  stop,'civ_sensitivity_gz: about to return'
  return,gz                     ;2-element array
end                             ; civ_sensitivity_gz


pro civ_sensitivity_mc, cmplt_fil, ostrct_fil, siiv=siiv,$
                        list=list, zlim=zlim, lyaforest=lyaforest,$
                        binncolm=binncolm, binew=binew, $
                        gallin_fil=gallin_fil, dvgal=dvgal,$
                        dvqso=dvqso,_extra=extra, bingz=bingz
  maxno = 9999.d                  ; some large number
  
  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'
  civ = dblt_retrieve(dblt_name)

  ;; Read file
  if keyword_set(list) then begin
     if keyword_set(dvqso) then $
        readcol,cmplt_fil,fil,zqso,format='a,f',comment=';' $
     else readcol,cmplt_fil,fil,format='a',comment=';' 
  endif else fil = cmplt_fil
  fil = fil[sort(fil)]
  nfil = n_elements(fil)

  ;; Truncate to just Lya-forest coverage
  if keyword_set(lyaforest) then begin
     if not keyword_set(zqso) then $
        stop,'civ_sensitivity_mc: QSO redshifts unknown'
     ;; Reset zqso and dvqso
     lya = dblt_retrieve('Lya')
     zqso = lya.wvI*(1.+zqso) / civ.wvI - 1.
     gd = where(zqso gt 0.,ngd)
     if ngd ne 0 then begin
        zqso = zqso[gd]
        fil = fil[gd]
        nfil = ngd
     endif 
  endif                         ; /lyaforest

  if not keyword_set(bingz) then bingz = 0.05

  for ff=0,nfil-1 do begin
     ;; Read in
     strct = xmrdfits(fil[ff],2,/silent) ; completeness structure
     ;; Eliminate masked out pixels immediately
     gd = where(strct.cmplt95[*,0] lt maxno and strct.cmplt95[*,1] lt maxno,$
                ngd,complement=bd,ncomplement=nbd)
     if ngd eq 0 then begin
        print,'civ_sensitivity_mc: No good data ',fil[ff]
        continue
     endif
     if nbd ne 0 then begin
        strct = {locz:strct.locz[gd], dz:strct.dz[gd],$
                 cmplt90:strct.cmplt90[gd,*], cmplt95:strct.cmplt95[gd,*]}   
        print,'civ_sensitivity_mc: Blacked out dz = ',total(strct.dz[bd]),$
              ' for ',fil[ff]
     endif 
     if keyword_set(zlim) then begin
        gd = where(strct.locz ge zlim[0] and strct.locz lt zlim[1],ngd)
        if ngd eq 0 then continue
        strct = {locz:strct.locz[gd], dz:strct.dz[gd],$
                 cmplt90:strct.cmplt90[gd,*], cmplt95:strct.cmplt95[gd,*]}
     endif                      ; zlim=
     if keyword_set(dvqso) then begin
        if not keyword_set(zqso) then $
           stop,'civ_sensitivity_mc: QSO redshifts unknown'
        if dvqso eq 1 then $
           gd = where(strct.locz lt zqso[ff],ngd) $
        else $
           gd = where(strct.locz lt zqso[ff]-dvqso/2.998e5*(1+zqso[ff]),ngd)
        if ngd eq 0 then continue
        strct = {locz:strct.locz[gd], dz:strct.dz[gd],$
                 cmplt90:strct.cmplt90[gd,*], cmplt95:strct.cmplt95[gd,*]}
     endif                      ; zlim=
     nbinz = n_elements(strct.locz)

     ;; Shift bin centers 
     strct.locz = strct.locz + 0.5*strct.dz
     ;; Compute X
     locx = cosm_xz(strct.locz,/silent,_extra=extra) ; cosmology
     dx = cosm_xz(strct.locz + 0.5*strct.dz,zmin=strct.locz - 0.5*strct.dz,$
                  /silent,_extra=extra)

     ;; Concatenate results
     if not keyword_set(ncmplt90_all) then begin
        ncmplt90_all = strct.cmplt90[*,0]
        ewcmplt90_all = strct.cmplt90[*,1]
        ncmplt95_all = strct.cmplt95[*,0]
        ewcmplt95_all = strct.cmplt95[*,1]
        locz_all = strct.locz
        dz_all = strct.dz
        locx_all = locx
        dx_all = dx
        indx_all = replicate(ff,nbinz)
     endif else begin
        ncmplt90_all = [ncmplt90_all,strct.cmplt90[*,0]]
        ewcmplt90_all = [ewcmplt90_all,strct.cmplt90[*,1]]
        ncmplt95_all = [ncmplt95_all,strct.cmplt95[*,0]]
        ewcmplt95_all = [ewcmplt95_all,strct.cmplt95[*,1]]
        locz_all = [locz_all,strct.locz]
        dz_all = [dz_all,strct.dz]
        locx_all = [locx_all,locx]
        dx_all = [dx_all,dx]
        indx_all = [indx_all,replicate(ff,nbinz)]
     endelse 

  endfor                        ; loop nfil
  n_all = n_elements(ncmplt90_all)

  if keyword_set(gallin_fil) then begin
     ;; Black out galactic lines
     if not keyword_set(dvgal) then dvgal = 0.

     ;; Re-sample data
     nsmpl = 10L                ; looks sufficient
     n_tmp = n_all * nsmpl
     ncmplt90_tmp = dblarr(n_tmp)
     ewcmplt90_tmp = dblarr(n_tmp)
     ncmplt95_tmp = dblarr(n_tmp)
     ewcmplt95_tmp = dblarr(n_tmp)
     locz_tmp = dblarr(n_tmp)
     dz_tmp = dblarr(n_tmp)
     indx_tmp = lonarr(n_tmp)

     for ii=0L,nsmpl-1 do begin
        rng = ii + lindgen(n_all)*nsmpl
        ncmplt90_tmp[rng] = ncmplt90_all
        ewcmplt90_tmp[rng] = ewcmplt90_all
        ncmplt95_tmp[rng] = ncmplt95_all
        ewcmplt95_tmp[rng] = ewcmplt95_all
        dz_tmp[rng] = dz_all/double(nsmpl)
        locz_tmp[rng] = locz_all - 0.5*dz_all + (ii+0.5)*dz_tmp[rng]
        indx_tmp[rng] = indx_all
     endfor 
     ;; Compute X
     locx_tmp = cosm_xz(locz_tmp,/silent,_extra=extra) ; cosmology
     dx_tmp = cosm_xz(locz_tmp + 0.5*dz_tmp,zmin=locz_tmp - 0.5*dz_tmp,$
                      /silent,_extra=extra)

     ;; Overwrite
     n_all = n_tmp
     ncmplt90_all = ncmplt90_tmp 
     ewcmplt90_all = ewcmplt90_tmp 
     ncmplt95_all = ncmplt95_tmp 
     ewcmplt95_all = ewcmplt95_tmp 
     locz_all = locz_tmp 
     dz_all = dz_tmp 
     locx_all = locz_tmp
     dx_all = dx_tmp
     indx_all = indx_tmp

     ;; Read in
     if size(gallin_fil,/type) eq 7 then $
        gallin = xmrdfits(gallin_fil,1,/silent) $
     else gallin = gallin_fil
     nlin = n_elements(gallin)
     gd = where(gallin.wv_lim[0]/civ.wvI-1. ge min(locz_all,max=mx) and $
                gallin.wv_lim[1]/civ.wvII-1. le mx,nlin)
     gallin = gallin[gd]

     mask = replicate(0L,n_all)
     ;; Set 1548 and 1550 effects (including dvgal)
     zlimI = dblarr(nlin,2)
     zlimII = dblarr(nlin,2)
     if dvgal lt 1. and dvgal gt 0. then begin
        ;; Use fraction
        dzI = (gallin.wv_lim[1]-gallin.wv_lim[0])/civ.wvI
        dzII = (gallin.wv_lim[1]-gallin.wv_lim[0])/civ.wvII
        zlimI[*,0] = gallin.wv_lim[0]/civ.wvI - 1. - dvgal*dzI
        zlimI[*,1] = gallin.wv_lim[1]/civ.wvI - 1. + dvgal*dzI
        zlimII[*,0] = gallin.wv_lim[0]/civ.wvII - 1. - dvgal*dzII
        zlimII[*,1] = gallin.wv_lim[1]/civ.wvII - 1. + dvgal*dzII
     endif else begin
        ;; Flat cut
        zlimI[*,0] = (-dvgal/3.e5 + 1.) * gallin.wv_lim[0]/civ.wvI - 1.
        zlimI[*,1] = (dvgal/3.e5 + 1.) * gallin.wv_lim[1]/civ.wvI - 1.
        zlimII[*,0] = (-dvgal/3.e5 + 1.) * gallin.wv_lim[0]/civ.wvII - 1.
        zlimII[*,1] = (dvgal/3.e5 + 1.) * gallin.wv_lim[1]/civ.wvII - 1.
     endelse 

     for ii=0,nlin-1 do begin
        bd = where((locz_all ge zlimI[ii,0] and locz_all le zlimI[ii,1]) or $
                   (locz_all ge zlimII[ii,0] and locz_all le zlimII[ii,1]),nbd)
        if nbd ne 0 then mask[bd] = mask[bd] + 1
     endfor                     ; loop nlin

     gz = replicate(0.d,n_all,2)
     gz[*,0] = locz_all
     gd = where(mask eq 0,n_all,complement=bd)
     if n_all ne 0 then begin
        ncmplt90_all = ncmplt90_all[gd]
        ewcmplt90_all = ewcmplt90_all[gd]
        ncmplt95_all = ncmplt95_all[gd]
        ewcmplt95_all = ewcmplt95_all[gd]
        locz_all = locz_all[gd]
        dz_all = dz_all[gd]
        locx_all = locx_all[gd]
        dx_all = dx_all[gd]
        indx_all = indx_all[gd]
        gz[gd,1] = 1            ; open pathlength
     endif else $
        stop,'civ_sensitivity_mc: everything blacked out by Galactic lines'
     srt = sort(gz[*,0])
     gz[*,0] = gz[srt,0]
     gz[*,1] = gz[srt,1]
  endif else begin                        ; gallin_fil=
     gz = replicate(0.d,n_all,2)
     gz[*,0] = locz_all[sort(locz_all)]
     gz[*,1] = 1 
  endelse 

  ;; Proper g(z) + buffer
  nbingz = ceil((max(locz_all,min=mn)-mn)/bingz) + 2
  cumgz = dblarr(nbingz,2)
  cumgz[*,0] = mn + (dindgen(nbingz)-1)*bingz
  for gg=0,nbingz-1 do begin
     gd = where(locz_all ge cumgz[gg,0] and $
                locz_all lt cumgz[gg,0] + bingz,ngd)
     if ngd eq 0 then continue
     unq = uniq(indx_all[gd],sort(indx_all[gd]))
     cumgz[gg,1] = n_elements(unq)
  endfor                        ; loop nbingz
  

  ;; Column density
  if not keyword_set(binncolm) then binncolm = 0.01
  bd = where(ncmplt90_all eq 0. or ncmplt95_all eq 0. or $
            ncmplt90_all eq maxno or ncmplt95_all eq maxno)
  if bd[0] ne -1 then $
     print,'civ_sensitivity_mc: column density dX blacked out = ',$
           total(dx_all[bd])
  tmp = [ncmplt90_all,ncmplt95_all]
  gd = where(tmp gt 0. and tmp lt maxno)
  nlim = [min(tmp[gd],max=mx),mx]
  nbinncolm = ceil(nlim[1]-nlim[0])/binncolm
  locncolm = nlim[0] + lindgen(nbinncolm)*binncolm
  x_ncolm = dblarr(nbinncolm,2)
  z_ncolm = dblarr(nbinncolm,2)
  for ii=0,nbinncolm-1 do begin
     gd = where(ncmplt90_all ge locncolm[ii] and $
                ncmplt90_all lt locncolm[ii]+binncolm,ngd)
     if ngd ne 0 then begin
        x_ncolm[ii,0] = total(dx_all[gd])
        z_ncolm[ii,0] = total(dz_all[gd])
     endif 
     gd = where(ncmplt95_all ge locncolm[ii] and $
                ncmplt95_all lt locncolm[ii]+binncolm,ngd)
     if ngd ne 0 then begin
        x_ncolm[ii,1] = total(dx_all[gd])
        z_ncolm[ii,1] = total(dz_all[gd])
     endif 
  endfor                        ; loop nbinncolm
  cumx_ncolm = dblarr(nbinncolm,2)
  cumx_ncolm[*,0] = total(x_ncolm[*,0],/cum)
  cumx_ncolm[*,1] = total(x_ncolm[*,1],/cum)
  cumz_ncolm = dblarr(nbinncolm,2)
  cumz_ncolm[*,0] = total(z_ncolm[*,0],/cum)
  cumz_ncolm[*,1] = total(z_ncolm[*,1],/cum)

  ;; Equivalent width
  if not keyword_set(binew) then binew = 1. ; mA
  bd = where(ewcmplt90_all eq 0. or ewcmplt95_all eq 0.)
  if bd[0] ne -1 then $
     print,'civ_sensitivity_mc: EW dX blacked out = ',$
           total(dx_all[bd])
  tmp = [ewcmplt90_all,ewcmplt95_all]
  gd = where(tmp gt 0. and tmp lt maxno)
  ewlim = [min(tmp[gd],max=mx),mx]
  nbinew = ceil(ewlim[1]-ewlim[0])/binew
  locew = ewlim[0] + lindgen(nbinew)*binew
  x_ew = dblarr(nbinew,2)
  z_ew = dblarr(nbinew,2)
  for ii=0,nbinew-1 do begin
     gd = where(ewcmplt90_all ge locew[ii] and $
                ewcmplt90_all lt locew[ii]+binew,ngd)
     if ngd ne 0 then begin
        x_ew[ii,0] = total(dx_all[gd])
        z_ew[ii,0] = total(dz_all[gd])
     endif 
     gd = where(ewcmplt95_all ge locew[ii] and $
                ewcmplt95_all lt locew[ii]+binew,ngd)
     if ngd ne 0 then begin
        x_ew[ii,1] = total(dx_all[gd])
        z_ew[ii,1] = total(dz_all[gd])
     endif 
  endfor                        ; loop nbinew
  cumx_ew = dblarr(nbinew,2)
  cumx_ew[*,0] = total(x_ew[*,0],/cum)
  cumx_ew[*,1] = total(x_ew[*,1],/cum)
  cumz_ew = dblarr(nbinew,2)
  cumz_ew[*,0] = total(z_ew[*,0],/cum)
  cumz_ew[*,1] = total(z_ew[*,1],/cum)


  ;; Write out results
  strct = {locn:locncolm, cumx_ncolm:cumx_ncolm, cumz_ncolm:cumz_ncolm, $
           locew:locew, cumx_ew:cumx_ew, cumz_ew:cumz_ew, $
           gz:gz, cumgz:cumgz}
  mwrfits,strct,ostrct_fil,/create,/silent
  print,'civ_sensitivity_mc: created ', ostrct_fil

end                             ; civ_sensitivity_mc


function civ_sensitivity_x, sens_fil, sigx=sigx, ncolm=ncolm, signcolm=signcolm, $
                            ew=ew, sigew=sigew, z=z, sigz=sigz
  ;; 95% confidence limit

  if size(sens_fil,/type) eq 7 then sens=xmrdfits(sens_fil,1,/silent) $
  else sens = sens_fil
  
  if not keyword_set(ncolm) and not keyword_set(ew) then begin
     print,'civ_sensitivity_x(): column density nor EW set'
     return, -1
  endif 

  if keyword_set(ncolm) then num = n_elements(ncolm) $
  else num = n_elements(ew)
  nbinn = n_elements(sens.locn)
  nbinew = n_elements(sens.locew)

  x = dblarr(num)
  sigx = dblarr(num,2)
  z = dblarr(num)
  sigz = dblarr(num,2)

  ;; Ease the looping (for now, just greatly help simulations)
  if keyword_set(ncolm) then begin
     ;; Find where value plateaus
     gd = where(sens.cumx_ncolm[0,1] eq sens.cumx_ncolm[*,1],ngd)
     if ngd eq n_elements(sens.cumx_ncolm[*,1]) then begin
        x = sens.cumx_ncolm[0,1]
        sigx[*,0] = 0.
        sigx[*,1] = 0.
        z = sens.cumz_ncolm[0,1]
        sigz[*,0] = 0.
        sigz[*,1] = 0.
        num = 0                 ; skip loop
     endif 
  endif else begin
     ;; Find where value plateaus
     gd = where(sens.cumx_ew[0,1] eq sens.cumx_ew[*,1],ngd)
     if ngd eq n_elements(sens.cumx_ncolm[*,1]) then begin
        x = sens.cumx_ew[0,1]
        sigx[*,0] = 0.
        sigx[*,1] = 0.
        z = sens.cumz_ew[0,1]
        sigz[*,0] = 0.
        sigz[*,1] = 0.
        num = 0                 ; skip loop
     endif 
  endelse 

  for ii=0,num-1 do begin
     ;; Interpolate 
     if keyword_set(ncolm) then begin
        ;; Column Density
        if ncolm[ii] lt sens.locn[0] then begin
           ;; Floor
           x[ii] = 0. 
           z[ii] = 0. 
           ;; Error
           sigx[ii,1] = sens.cumx_ncolm[0,1] ; max upper; no lower
           sigz[ii,1] = sens.cumz_ncolm[0,1] ; max upper; no lower
        endif else begin
           if ncolm[ii] gt sens.locn[nbinn-1] then begin
              ;; Ceiling
              x[ii] = sens.cumx_ncolm[nbinn-1,1] 
              z[ii] = sens.cumz_ncolm[nbinn-1,1] 
              ;; Error
              gd = where(sens.cumx_ncolm[*,1] lt sens.cumx_ncolm[nbinn-1,1],ngd)
              if ngd ne 0 then begin
                 ;; Function may be flat
                 sigx[ii,0] = x[ii] - sens.cumx_ncolm[gd[ngd-1],1] ; max lower; no upper
                 sigz[ii,0] = z[ii] - sens.cumz_ncolm[gd[ngd-1],1] ; max lower; no upper
              endif 
           endif else begin
              ;; Interpolate
              rslt = interpol(sens.cumx_ncolm[*,1], sens.locn, ncolm[ii])
              rsltz = interpol(sens.cumz_ncolm[*,1], sens.locn, ncolm[ii])
              x[ii] = rslt
              z[ii] = rsltz
              ;; Error
              mn = min(sens.cumx_ncolm[*,1]-x[ii],imn,/absolute)
              if x[ii] lt sens.cumx_ncolm[imn,1] then begin
                 if imn gt 0 then begin
                    sigx[ii,0] = x[ii] - sens.cumx_ncolm[imn-1,1]
                    sigz[ii,1] = sens.cumz_ncolm[imn,1] - z[ii]
                 endif 
                 sigx[ii,1] = sens.cumx_ncolm[imn,1] - x[ii]
                 sigz[ii,1] = sens.cumz_ncolm[imn,1] - z[ii]
              endif else begin
                 if imn lt nbinn-1 then begin
                    sigx[ii,0] = sens.cumx_ncolm[imn+1,1] - x[ii]
                    sigz[ii,0] = sens.cumz_ncolm[imn+1,1] - z[ii]
                 endif 
                 sigx[ii,1] = x[ii] - sens.cumx_ncolm[imn,1]
                 sigz[ii,1] = z[ii] - sens.cumz_ncolm[imn,1]
              endelse 
          endelse              ; interpolate
        endelse                 ; not < min locn 
     endif else begin           ; end ncolm=
        ;; Equivalent width
        if ew[ii] lt sens.locew[0] then begin
           ;; Floor
           x[ii] = 0. 
           z[ii] = 0. 
           ;; Error
           sigx[ii,1] = sens.cumx_ew[0,1] ; max upper; no lower
           sigz[ii,1] = sens.cumz_ew[0,1] ; max upper; no lower
        endif else begin
           if ew[ii] gt sens.locew[nbinew-1] then begin
              ;; Ceiling
              x[ii] = sens.cumx_ew[nbinew-1,1] 
              z[ii] = sens.cumz_ew[nbinew-1,1] 
              ;; Error
              gd = where(sens.cumx_ew[*,1] lt sens.cumx_ew[nbinew-1,1],ngd)
              if ngd ne 0 then begin
                 ;; Function may be flat
                 sigx[ii,0] = x[ii] - sens.cumx_ew[gd[ngd-1],1] ; max lower; no upper
                 sigz[ii,0] = z[ii] - sens.cumz_ew[gd[ngd-1],1] ; max lower; no upper
              endif 
           endif else begin
              ;; Interpolate
              rslt = interpol(sens.cumx_ew[*,1], sens.locew, ew[ii])
              rsltz = interpol(sens.cumz_ew[*,1], sens.locew, ew[ii])
              x[ii] = rslt
              z[ii] = rsltz
              mn = min(sens.cumx_ew[*,1]-x[ii],imn,/absolute)
              if x[ii] lt sens.cumx_ew[imn,1] then begin
                 if imn gt 0 then begin
                    sigx[ii,0] = x[ii] - sens.cumx_ew[imn-1,1]
                    sigz[ii,0] = z[ii] - sens.cumz_ew[imn-1,1]
                 endif 
                 sigx[ii,1] = sens.cumx_ew[imn,1] - x[ii]
                 sigz[ii,1] = sens.cumz_ew[imn,1] - z[ii]
              endif else begin
                 if imn lt nbinew-1 then begin
                    sigx[ii,0] = sens.cumx_ew[imn+1,1] - x[ii]
                    sigz[ii,0] = sens.cumz_ew[imn+1,1] - z[ii]
                 endif 
                 sigx[ii,1] = x[ii] - sens.cumx_ew[imn,1]
                 sigz[ii,1] = z[ii] - sens.cumz_ew[imn,1]
              endelse 
           endelse              ; interpolate
        endelse                 ; not < min locew   
     endelse                    ; end ew=
  endfor                        ; loop num

  ;; Account for error in measurement
  if keyword_set(signcolm) then begin
     ;; Treat ncolm as log, always
     xhi = civ_sensitivity_x(sens,$
                             ncolm=(ncolm+alog10(abs(1.+alog(10.)*signcolm))), $
                             z=zhi)
     xlo = civ_sensitivity_x(sens,$
                             ncolm=(ncolm+alog10(abs(1.-alog(10.)*signcolm))), $
                             z=zlo)
     ;; Add in quadrature
     sigx[*,0] = sqrt(sigx[*,0]^2 + (x - xlo)^2)
     sigx[*,1] = sqrt(sigx[*,1]^2 + (xhi - x)^2)
     sigz[*,0] = sqrt(sigz[*,0]^2 + (z - zlo)^2)
     sigz[*,1] = sqrt(sigz[*,1]^2 + (zhi - z)^2)
  endif 
  if keyword_set(sigew) then begin
     xhi = civ_sensitivity_x(sens,ew=(ew+sigew), z=zhi)
     xlo = civ_sensitivity_x(sens,ew=(ew-sigew), z=zlo)
     ;; Add in quadrature
     sigx[*,0] = sqrt(sigx[*,0]^2 + (x - xlo)^2)
     sigx[*,1] = sqrt(sigx[*,1]^2 + (xhi - x)^2)
     sigz[*,0] = sqrt(sigz[*,0]^2 + (z - zlo)^2)
     sigz[*,1] = sqrt(sigz[*,1]^2 + (zhi - z)^2)
  endif 

  return, x
  
end                             ; civ_sensitivity_x()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_sensitivity, strct_fil, wrest, dopb, nsigma, ewlim, $
                     outroot=outroot

if (n_params() lt 5) then begin
    print,'Syntax -'+$
      'civ_sensitivity, strct_fil, wrest, dopb, nsigma, ewlim '
    print,'     [outroot=,ndopb=]'
    return
endif 

;; Determie root directory
mlss_dir = getenv('MLSS_DIR')
if strtrim(mlss_dir,2) eq '' then $
   mlss_dir = strtrim(getenv('HOME'),2)+'/MLSS'
mlss_dir = mlss_dir + '/'

;; Open structure
if size(strct_fil,/type) ne 7 then $
  stop,'civ_sensitivity: strct_fil must be name'

strct = xmrdfits(strct_fil, 1, /silent) 
nstrct = n_elements(strct)

unq = uniq(strtrim(strct.instr_fil,2),sort(strct.instr_fil))
instr_list = mlss_dir+strtrim(strct[unq].instr_fil,2)
nlist = n_elements(instr_list)

nlin = n_elements(wrest)
if nlin gt 2 then $
  stop,'civ_sensitivity: only works for single line or doublet',wrest
ncolm = fltarr(nlin)

;; Loop different LOS (instrument files)
for ll=0,nlist-1 do begin
    readcol,instr_list[ll],instr_fil, inst_dw, inst_w0,format='a,f,f'
    instr_fil = mlss_dir + instr_fil
    ninst = n_elements(instr_fil)
    if ninst lt 7 then stop,'civ_sensitivity: too few spectra'

    ;; Loop instruments
    for ii=0,ninst-1 do begin
        ;; If spectrum doesn't exist continue
        test = file_search(instr_fil[ii],count=ntest)
        if ntest eq 0 then goto,next_instr ;nothing to do

        if ii LE 6 then $
          fx = x_readspec(instr_fil[ii], sig=sig, wav=wave, $
                          NPIX=npix, inflg=3)$
        else begin              ; HST
            spos = strpos(instr_fil[ii], 'f.fits')
            sig_fil = strmid(instr_fil[ii], 0, spos)+'e.fits'
            fx = x_readspec(instr_fil[ii], SIG=sig, wav=wave, $
                            NPIX=npix,fil_sig=sig_fil, inflg=0)
        endelse
        ;; Shift wave
        wave = wave + inst_dw[ii]*inst_w0[ii]/wave
        
        ;; Test if instr_fil covers
        if max(wave)/min(wrest)-1. le 0. then goto,next_instr

        dwv = wave - shift(wave,1)
        dwv[0] = dwv[1]


        ;; Cut the overhead 
        bd = where(sig le 0.,nbd)
        bdsig = replicate(0,npix)
        if nbd ne 0 then begin
           bdsig[bd] = 1        ;num pix
           sig[bd] = 0.         ;make calc right
        endif 
        sigdwv2 = (sig*dwv)^2
        cumsigdwv2 = total((sig*dwv)^2,/cumulative)
        cumbdsig = total(bdsig,/cumulative)
        sigew = fltarr(nlin,npix)

        for jj=0L,nlin-1 do begin
            wlo = wave - 0.5*dopb/2.998e5*wrest[jj]
            whi = wave + 0.5*dopb/2.998e5*wrest[jj]
            
            ;; Flag regions not broad enough
            bd = where(wlo lt min(wave),nbd)
            if nbd ne 0 then sigew[jj,bd] = -1
            bd = where(whi gt max(wave),nbd)
            if nbd ne 0 then sigew[jj,bd] = -1
            zabs = wave/wrest[jj] - 1. > 0. ;floor
            gdz = where(sigew[jj,*] ge 0,ngdz)
            for kk=0L,ngdz-1 do begin
                gg = gdz[kk]
                mn = min(wave-wlo[gg],ilo,/absolute)
                mn = min(wave-whi[gg],ihi,/absolute)
                if ihi eq ilo then stop,'civ_sensitivity: no data in window'
                if (cumbdsig[ihi]-cumbdsig[ilo])/double(ihi-ilo+1) gt 0.1 then $
                   sigew[jj,gg] = -1 $ ;>10% bad 
                else sigew[jj,gg] = sqrt(cumsigdwv2[ihi]-cumsigdwv2[ilo])/$
                                    (1+zabs[gg])*1000. ;rest EW sig (mA)
            endfor              ;loop npix

            ncolm[jj] = civ_sensitivity_fndcolm(wrest[jj], dopb,ewlim)
        endfor                  ; loop nlin (wrest)

        ;; Determine g(z)
        gz = civ_sensitivity_gz(wave,sigew,wrest,nsigma,ewlim,zlim=zlim,$
                               zpath=zpath)

        getfnam,wrest[0],fval,nam
        nam = strsplit(strtrim(nam,2),/extract)
        if keyword_set(outroot) then begin
           ;; e.g. root_128gzHI80mA.fits
           ofil = outroot+strtrim(2^ii,2)+'_gz'+nam[0]+$
             strtrim(round(dopb),2)+'kms'+$
                  strtrim(round(ewlim),2)+'mA.fits'
        endif else begin
           ;; e.g. qso_E140M_gzHI80mA.fits
            if ii le 6 then $   ; FUSE
              ofil = strmid(instr_fil[ii],0,$
                            strpos(instr_fil[ii],'.',/reverse_search)) $
            else $ ; HST
              ofil = strmid(instr_fil[ii],0,$
                            strpos(instr_fil[ii],'_f.',/reverse_search))
            ofil = ofil+'_gz'+nam[0]+strtrim(round(dopb),2)+'kms'+$
              strtrim(round(ewlim),2)+'mA.fits'
        endelse

        ;; Create structure
        sigstrct = {civcand_fil:strct_fil, $
                    wrest:wrest, $ ;rest wavelength (Ang)
                    dopb:dopb, $ ;Doppler parameter/width (km/s)
                    nsigma:nsigma, $ ;EW/sigEW significance limit
                    ewlim:ewlim, $ ;EW limit
                    ncolm:ncolm, $ ;from dobp and ewlim
                    wave:wave, $           ; --> zabs, dz
                    sigew:sigew, $ ;error in window
                    gz:gz, $ ;sensitivity redshift & curve
                    zlim:zlim, $  ; bounds >0 up to max
                    zpath:zpath $ ;unblocked pathlength
                   } 
        mwrfits,sigstrct,ofil,/create,/silent
        print,'civ_sensitivity: created ',ofil

        ;; If instr_fil doesn't exist or doesn't cover redshift
        next_instr:
    endfor                      ;loop ninst
endfor                          ;loop nlist
end
