;+ 
; NAME:
; civ_aodm
;   Version 1.1
;
; PURPOSE:
;    Calculates an AODM column density profile given a civcandstrct (with
;   instr_fil, zabs, wrest value set).  
;
; CALLING SEQUENCE:
;   
;   civ_aodm, strct, /LOG
;
; INPUTS:
;   strct_fil -- FITS file of civcandstrct, IDL save file of civcand,
;                or actual civcandstrct
;   [instr_list -- list of instruments]
;
; RETURNS:
;
; OUTPUTS:
;   strct_fil -- overwrites with aodm_civ, sigaodm_civ, aodm_vel,
;                flg_sat in whatever format supplied
;
; OPTIONAL KEYWORDS:
;   /LOG - Return log answers
;   nbin - number of pixels overwhich to smooth (should be odd)
;
; OPTIONAL OUTPUTS:
;   FLG_SAT=  Keyword indicating whether the transtition was saturated
;             (also already written to strct_fil upon output)
;   flg_mtch= fraction profile agreement (1-sigma agreement)
;   aodm_vel= in addition to writing to strcture, return results
;   aodm_civ= in addition to writing to structure, return results 
;   sigaodm_civ= in addition to writing to structure, return results
;
; COMMENTS:
;
; EXAMPLES:
;   civ_aodm, strct, /log, nbin=3
;
; PROCEDURES CALLED:
;   smooth
;   dblt_retrieve
;   x_readspec
;
; REVISION HISTORY:
;   21-Feb-2008  created by KLC, copied from XIDL x_aodm
;   26-Feb-2008  each aodm_vel locked to measured centroid of feature
;    4-mar-2008  added civ_aodm_mtch for consistency
;-
;------------------------------------------------------------------------------
@hdf5_linlst                    ; resolve hdf5_linlst_ion()
function civ_aodm_mtch,strct
  ;; will fail if strct more than one element
  gd = where(strtrim(strct.ion,2) ne '')
  prs = strsplit(strct.ion[gd[0]],/extract,count=nprs) 
  dblt_name = prs[0]
  
  civ = where(stregex(strct.ion,dblt_name,/boolean),nciv)
  if nciv ne 2 then stop,'civ_aodm_mtch: no '+dblt_name+' doublet'

  civ = civ[sort(strct.wrest[civ])]

  ;; Handle hdf5_fil, which may not have error array
  test = where(strct.sigaodm_civ eq 0.,ntest)
  if ntest eq n_elements(strct.sigaodm_civ) then begin
;     print,'civ_aodm_mtch(): no AOD profile error'
     return,1.
  endif 

  for jj=0,nciv-1 do begin
     vlim = 2.998e5*(strct.wv_lim[civ[jj],*]/$
                     (strct.wrest[civ[jj]]*$
                      (1+strct.zabs[civ[0]]))-1.)
     if jj eq 0 then vmnx = vlim $
     else begin
        tmp = [vmnx[0] < vlim[0],vmnx[0] > vlim[1]] ;largest range
        vmnx = tmp
     endelse 
  endfor                        ;loop doublet
  mn = min(strct.aodm_vel[*,0]-vmnx[0],imn1548,/absolute)
  mx = min(strct.aodm_vel[*,0]-vmnx[1],imx1548,/absolute)
  mn = min(strct.aodm_vel[*,1]-vmnx[0],imn1550,/absolute)
  mx = min(strct.aodm_vel[*,1]-vmnx[1],imx1550,/absolute)
  cmp = strct.aodm_civ[imn1550:imx1550,1] - $
        strct.aodm_civ[imn1548:imx1548,0]
  sigcmp = sqrt(strct.sigaodm_civ[imn1550:imx1550,1]^2 + $
                strct.sigaodm_civ[imn1548:imx1548,0]^2)
  gd = where(abs(cmp) le sigcmp,ngd)
  mx = imx1548-imn1548 > imx1550-imn1550
  return,ngd/double(mx+1)
end





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_aodm, strct_fil, instr_list, savfil=savfil, log=log, nbin=nbin,$
              flg_sat=flg_sat, flg_mtch=flg_mtch, view=view, inspec=inspec, $
              aodm_vel=aodm_vel, aodm_civ=aodm_civ, sigaodm_civ=sigaodm_civ,$
              dblt_name=dblt_name, hdf5_fil=hdf5_fil, debug=debug,_extra=extra

  if (N_params() LT 1) then begin 
     print,'Syntax - ' + $
           'civ_aodm, strct_fil, FLG_SAT=,' + $
           ' VELO= [v1.1]'
     return
  endif 

  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_aodm: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent)
  endif else strct = strct_fil
  strct.instr_fil = strtrim(strct.instr_fil,2)
  nstrct = n_elements(strct)
  flg_mtch = fltarr(nstrct)
  flg_sat = fltarr(nstrct,2)

  if not keyword_set(dblt_name) then begin
     ;; Assume first filled-in element is structure focus
     gd = where(strtrim(strct[0].ion,2) ne '',ngd)
     if ngd lt 2 then stop,'civ_aodm: no doublet'
     prs = strsplit(strct[0].ion[gd[0]],/extract,count=nprs) 
     dblt_name = prs[0]
  endif

  ;; Read instrument file list (same for all)
  if keyword_set(instr_list) then begin
     readcol, getenv('MLSS_DIR')+'/'+$
              instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f' 
     instr_fil = getenv('MLSS_DIR')+'/'+instr_fil
     nlist = n_elements(instr_fil)
     if nlist LT 7 then stop,'civ_aodm: too few spectra'
  endif 

  ;; Size as used in civcandstrct
  sza = size(strct[0].aodm_vel, /dim)
  nsvpx = sza[0]
  fxaod = fltarr(nsvpx,2)
  sigaod = fltarr(nsvpx,2)
  aodm_vel = fltarr(nstrct,nsvpx,2)
  aodm_civ = dblarr(nstrct,nsvpx,2)
  sigaodm_civ = dblarr(nstrct,nsvpx,2)

  ;; Get the info
  civinfo = dblt_retrieve(dblt_name)
  
  ;; Constant
  cst = (10.d^14.5761)/([civinfo.fI,civinfo.fII]*[civinfo.wvI,civinfo.wvII])
  fxflr = [0.05,exp(alog(0.05)*cst[0]/cst[1])] ;flux floor
  svinstr_list = ''

  if keyword_set(inspec) then begin
     fx = inspec.flux
     wave = inspec.wave
     sig = inspec.error
     npix = n_elements(fx)

     ;; Trim flux (reduces effect of spurious pixels)
     ;; Un-physical, indicates redux problem or echelle gap
     bd = where(sig le 0.,nbd,complement=gd,ncomplement=ngd) 
     if nbd ne 0 then begin
        sig[bd] = 0.
        fx[bd] = -1.
     endif 
     if ngd eq 0 then stop,'civ_aodm: no data with sig > 0'
     bd = where(fx[gd] lt -abs(sig[gd]),nbd)
     if nbd ne 0 then fx[gd[bd]] = -abs(sig[gd[bd]]) ;floor
     bd = where(fx[gd] gt 1+sig[gd],nbd)
     if nbd ne 0 then fx[gd[bd]] = 1+sig[gd[bd]] ;ceil
  endif 

  for ii=0,nstrct-1 do begin
     ;; Read instrument file list (when different)
     if not keyword_set(instr_list) and not keyword_set(inspec) and $
        strct[ii].instr_fil ne svinstr_list and $
        not keyword_set(hdf5_fil) then begin 
        readcol, getenv('MLSS_DIR')+'/'+$
                 strct[ii].instr_fil, instr_fil, inst_dw, inst_w0, $
                 format='a,f,f'
        instr_fil = getenv('MLSS_DIR')+'/'+instr_fil
        nlist = n_elements(instr_fil)
        if nlist LT 7 then stop,'civ_aodm: too few spectra'
        svinstr_list = strct[ii].instr_fil
     endif 

     ;; Find CIV 
     civ = where(stregex(strct[ii].ion,dblt_name,/boolean),nciv)
     if nciv ne 2 then stop,'civ_aodm: no '+dblt_name+' doublet to evaluate'
     srt = sort(strct[ii].wrest[civ])
     civ = civ[srt]

     ;; Prep arrays
     if keyword_set(hdf5_fil) then instr = strct[ii].instr[civ] $
     else begin
        lgv = alog(strct[ii].instr[civ])/alog(2)
        instr = fix(lgv+0.00001)
     endelse 
     svinst = -1
     for jj=0,nciv-1 do begin   ;loop nciv 
        if instr[jj] ne svinst and not keyword_set(inspec) then begin
           if keyword_set(hdf5_fil) then begin
              fx = hdf5_readspec(getenv('MLSS_DIR')+'/'+strct[ii].instr_fil,$
                                 instr[jj],wav=wave,sig=sig,npix=npix,$
                                 _extra=extra) ; /clean, snr=, ion_name=, /ion_only, etc
              gd = lindgen(npix)
              ngd = npix
              test = where(sig eq 0.,ntest)
              if ntest eq npix then clean = 1      ; flag for later
           endif else begin
              if instr[jj] LE 6 then $ ;FUSE
                 fx = x_readspec(instr_fil[instr[jj]],SIG=sig,wav=wave,$
                                 NPIX=npix,inflg=3)$
              else begin        ; HST
                 spos = strpos(instr_fil[instr[jj]], 'f.fits')
                 sig_fil = strmid(instr_fil[instr[jj]], 0, spos)+'e.fits'
                 fx = x_readspec(instr_fil[instr[jj]], SIG=sig, wav=wave, $
                                 NPIX=npix, fil_sig=sig_fil, inflg=0)
              endelse
              ;; Shift wave
              wave = wave + inst_dw[instr[jj]]*inst_w0[instr[jj]]/wave

              ;; Trim flux (reduces effect of spurious pixels)
              ;; Un-physical, indicates redux problem or echelle gap
              bd = where(sig le 0.,nbd,complement=gd,ncomplement=ngd) 
              if nbd ne 0 and nbd ne npix then begin
                 sig[bd] = 0.
                 fx[bd] = -1.
              endif 
              if ngd eq 0 then stop,'civ_aodm: no data with sig > 0'
           endelse 

           ;; Still trim to treat like data
           bd = where(fx[gd] lt -abs(sig[gd]),nbd)
           if nbd ne 0 then fx[gd[bd]] = -abs(sig[gd[bd]]) ;floor
           bd = where(fx[gd] gt 1+sig[gd],nbd)
           if nbd ne 0 then fx[gd[bd]] = 1+sig[gd[bd]] ;ceil
        endif                                          ; new spec to read

        ;; Align arrays (will have to adjust later to align)
        wcent = min(wave-strct[ii].wrest[civ[jj]]*(1+strct[ii].zabs[civ[jj]]),$
                    pcent,/absolute)
        pmn = pcent - nsvpx/2 > 0
        pmx = pcent + nsvpx/2 + 1 < npix-1
        imnx = nsvpx/2+ [-(pcent-pmn),(pmx-pcent)]
        aodm_vel[ii,imnx[0]:imnx[1],jj] = $
           2.998e5*(wave[pmn:pmx]-wave[pcent])/wave[pcent]
        ;; Pad (if necessary)
        if imnx[0] ne 0 then $
           aodm_vel[ii,0:imnx[0],jj] = aodm_vel[ii,imnx[0],jj]-$
                                       reverse(lindgen(imnx[0]+1))*$
                                       (aodm_vel[ii,imnx[0]+1,jj]-$
                                        aodm_vel[ii,imnx[0],jj])
        if imnx[1] ne nsvpx-1 then $
           aodm_vel[ii,imnx[1]+1:nsvpx-1,jj] = aodm_vel[ii,imnx[1],jj]+$
                                               lindgen(nsvpx-imnx[1]-1)*$
                                               (aodm_vel[ii,imnx[1],jj]-$
                                                aodm_vel[ii,imnx[1]-1,jj])

        fxaod[imnx[0]:imnx[1],jj] = fx[pmn:pmx]
        sigaod[imnx[0]:imnx[1],jj] = sig[pmn:pmx]
        svinst = instr[jj]      ;for next iteration comparison
     endfor                     ;loop nciv I

     if keyword_set(debug) then stop

     ;; Error catching
     if n_elements(fxaod) NE n_elements(sigaod) then $
        message, 'civ_aodm: flux and error arrays mis-match'


     ;; Calculate AODM profiles
     for jj=0,nciv-1 do begin   ;loop nciv II
        ;; Set nndt = ln(1/fx)*cst
        nndt = dblarr(nsvpx)

        ;; Deal with saturated pixels (f < 0.05 or 0.2*sig)
        if keyword_set(clean) then begin
           sat = -1
           nsat = 0
           unsat = lindgen(nsvpx)
           nunsat = nsvpx
        endif else $
           sat = where(fxaod[*,jj] LE sigaod[*,jj]/5. OR fxaod[*,jj] LT fxflr[0],$
                       nsat, complement=unsat, ncomplement=nunsat)
        if nsat NE 0 then begin
           ;; Floor flux to the corresponding limit (tau1550=tau1548)
           if keyword_set(clean) then gd = lindgen(nsat) $
           else gd = where(sigaod[sat,jj] GT 0., ngd)
           if ngd NE 0 then begin
              sub = fxflr[jj] > sigaod[sat[gd],jj]/5.d 
              nndt[sat[gd]] = alog(1./sub)*cst[jj]
           endif 
        endif

        ;; Unsaturated pixels
        if nunsat NE 0 then nndt[unsat] = alog(1.d/fxaod[unsat,jj])*cst[jj]

        ;; Arrays to return
        delv = aodm_vel[ii,*,jj] - shift(aodm_vel[ii,*,jj],1)
        delv[0] = delv[1]
        aodm_civ[ii,*,jj] = nndt*delv
        sigaodm_civ[ii,*,jj] = abs(delv*cst[jj]*sigaod[*,jj]/fxaod[*,jj])

        if keyword_set(debug) then stop

        ;; Perhaps due to padded profile
        bd = where(finite(aodm_civ[ii,*,jj]) eq 0 or $
                   finite(sigaodm_civ[ii,*,jj]) eq 0,nbd)
        if nbd ne 0 then begin
           aodm_civ[ii,bd,jj] = 0.
           sigaodm_civ[ii,bd,jj] = 0.
        endif 

        ;; Smooth (npix should be odd) and can make error be < 0.
        if keyword_set(nbin) then begin
           if nsvpx le nbin then stop,'civ_aodm: unable to smooth'
           tmpn = smooth(aodm_civ[ii,*,jj],nbin,/edge_truncate)
           tmps = smooth(sigaodm_civ[ii,*,jj],nbin,/edge_truncate)

           ;; Since can smoothing can make error < 0, reset
           bd = where(tmps le 0.,nbd)
           if nbd ne 0 and not keyword_set(clean) then begin
              tmpn[bd] = aodm_civ[ii,bd,jj]
              tmps[bd] = sigaodm_civ[ii,bd,jj]
           endif 
           
           bd = where(finite(tmpn) eq 0 or finite(tmps) eq 0,nbd)
           if nbd ne 0 then begin
              tmpn[bd] = 0.
              tmps[bd] = 0.
           endif 
           aodm_civ[ii,*,jj] = tmpn
           sigaodm_civ[ii,*,jj] = tmps
        endif                   ;nbin set

        if keyword_set(debug) then stop

        if keyword_set(log) then begin
           tmpn = aodm_civ[ii,*,jj]
           tmps = sigaodm_civ[ii,*,jj]
           ;; No need to make this sigaodm_civ cut, I believe
           if keyword_set(clean) then $ 
              bd = where(aodm_civ[ii,*,jj] lt 1.,nbd,$
                         complement=gd,ncomplement=ngd) $
           else $
              bd = where(aodm_civ[ii,*,jj] lt 1. or sigaodm_civ[ii,*,jj] lt 1.,$
                         nbd,complement=gd,ncomplement=ngd)
           if nbd ne 0 then begin
              tmpn[bd] = 0.
              tmps[bd] = 0.
           endif 
           if ngd ne 0 then begin
              tmpn[gd] = alog10(aodm_civ[ii,gd,jj])
              tmps[gd] = sigaodm_civ[ii,gd,jj]/abs(alog(10.)*aodm_civ[ii,gd,jj])
           endif 

           bd = where(finite(tmpn) eq 0 or finite(tmps) eq 0,nbd)
           if nbd ne 0 then begin
              tmpn[bd] = 0.
              tmps[bd] = 0.
           endif 
           aodm_civ[ii,*,jj] = tmpn
           sigaodm_civ[ii,*,jj] = tmps
        endif                   ;/log


        ;; Comparison information (largest range)
        if jj eq 0 then $
           vmnx = 2.998e5*(strct[ii].wv_lim[civ[jj],*]/$
                           (strct[ii].wrest[civ[jj]]*$
                            (1+strct[ii].zabs[civ[jj]]))-1.) $
        else begin
           tmp = 2.998e5*(strct[ii].wv_lim[civ[jj],*]/$
                          (strct[ii].wrest[civ[jj]]*$
                           (1+strct[ii].zabs[civ[jj]]))-1.)
           vmnx = [vmnx[0] < tmp[0], vmnx[1] > tmp[1]]
        endelse 
        flg_sat[ii,jj] = nsat
     endfor                     ;loop nciv II


     ;; Store data
     strct[ii].aodm_vel[*,0] = aodm_vel[ii,*,0]
     strct[ii].aodm_vel[*,1] = aodm_vel[ii,*,1]
     strct[ii].aodm_civ[*,0] = aodm_civ[ii,*,0]
     strct[ii].aodm_civ[*,1] = aodm_civ[ii,*,1]
     strct[ii].sigaodm_civ[*,0] = sigaodm_civ[ii,*,0]
     strct[ii].sigaodm_civ[*,1] = sigaodm_civ[ii,*,1]

     ;; Comparison information for binned data (flg_mtch)
     mn = min(strct[ii].aodm_vel[*,0]-vmnx[0],imn1548,/absolute)
     mx = min(strct[ii].aodm_vel[*,0]-vmnx[1],imx1548,/absolute)
     flg_sat[ii,0] = flg_sat[ii,0]/double(imx1548-imn1548+1)
     mn = min(strct[ii].aodm_vel[*,1]-vmnx[0],imn1550,/absolute)
     mx = min(strct[ii].aodm_vel[*,1]-vmnx[1],imx1550,/absolute)
     flg_mtch[ii] = ngd/double(mx+1)
     flg_mtch[ii] = civ_aodm_mtch(strct[ii])
     
     ;; View
     if keyword_set(view) then begin
        print,'civ_aodm: view '+dblt_name+' profiles z=',strct[ii].zabs[civ[0]]
        tmp = fx                ;b/c x_velplt modifies
        x_velplt,tmp,strct[ii].zabs[civ[0]],wave=wave,/sublls
        print,'civ_aodm: flg_sat 1548,1550: ',flg_sat[ii,0],flg_sat[ii,1]
        print,'civ_aodm: flg_mtch 1548-1550: ',flg_mtch[ii],$
              ngd,vmnx[1]-vmnx[0]+1
        ;; Plot profiles and indicate max bounds for comparison
        x_splot,strct[ii].aodm_vel[*,0],strct[ii].aodm_civ[*,0],psym1=10,$
                xtwo=strct[ii].aodm_vel[*,1],ytwo=strct[ii].aodm_civ[*,1],$
                psym2=10,xthr=[vmnx[0],vmnx[0],vmnx[1],vmnx[1]],$
                ythr=[0.,20.,20.,0.],psym3=10,/block
        stop,'civ_aodm: ... stop after viewing profiles'
     endif                      ; stop after both loaded
     
     ;; Reset
;     aodm_vel[ii,*,*] = 0.
     fxaod[*,*] = 0.
     sigaod[*,*] = 0.
;     aodm_civ[ii,*,*] = 0.
;     sigaodm_civ[ii,*,*] = 0.
  endfor                        ;loop nstrct

  ;; Overwrite
  if size(strct_fil,/type) eq 7 then begin
     spawn,'cp '+strct_fil+' '+strtrim(strct_fil,2)+'.bkp'
     if keyword_set(savfil) then begin
        civcand = strct
        save,civcand,filename=strct_fil
     endif else mwrfits, strct, strct_fil, /create,/silent
     print,'civ_aodm: overwrote ',strct_fil
  endif else strct_fil = strct

  return
end
