;+ 
; NAME:
; civ_calcewn
;  V1.1
;
; PURPOSE:
;    Given a CIV structure file, the code measures EW, N(X) and fills
;  up the CIV structure.  It takes the highest significance EW data
;  when multiple detections.  .
;
; CALLING SEQUENCE:
;   civ_calcewn, strct_fil, [instr_list]
;
; INPUTS:
;  strct_fil  -- IDL CIV file
;
; RETURNS:
;
; OUTPUTS:
;  strct_fil  -- The EW, ncolm, z, and errors as well as flg_colm 
;                are filled up by this routine
;
; OPTIONAL KEYWORDS:
;  savfil -- strct_fil is IDL save file
;  inspec -- single spectrum structure {wave, flux, error}
;  instr_flg -- instrument flag number (default to use all)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   civ_calcewn, struct
;
;
; PROCEDURES CALLED:
;   x_aodm
;
; REVISION HISTORY:
;   18-Feb-2008 created by KLC, adapted from XIDL fuse_calcewn
;   16-Dec-2008 added civ_calcewn_ndblt()
;    8-Apr-2009 streamline use with inspec=, instr_flg=
;-
;------------------------------------------------------------------------------

function civ_calcewn_ndblt,strct_fil,dblt_name,signcolm=signcolm,log=log,$
                           subscript=subscript,savfil=savfil, $
                           eq_weight=eq_weight, flg_colm=flg_colm,silent=silent
  ;; Return error-weighted average of two doublet values
  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_calcewn: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent,structyp='civcandstrct') 
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  ;; Doublet
  if keyword_set(dblt_name) then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_retrieve('CIV')
  if strtrim(dblt.ion,2) eq '' then begin
     print,'civ_caclewn_ndblt: no such doublet ',dblt_name
     signcolm = -1
     return,-1
  endif 

  ;; Error-weighted average of AODM column density
  ncolm = dblarr(nstrct)
  signcolm = dblarr(nstrct)
  flg_colm = replicate(1L,nstrct)
  subscript = replicate(-1,nstrct,2)

  for ii=0, nstrct-1 do begin
     ;; Find doublet
     gd = where(stregex(strct[ii].ion,dblt.ion,/boolean),ngd)
     
     case ngd of
        0: goto,skip_loop
        1: begin
           if dblt.ion eq 'HI' then begin
              ncolm[ii] = 10.^strct[ii].ncolm[gd]
              signcolm[ii] = alog(10.)*strct[ii].signcolm[gd]*ncolm[ii]
              flg_colm[ii] = strct[ii].flg_colm[gd]
              subscript[ii,*] = [gd[0],-1]
              goto,skip_loop
           endif 
        end
        2:                      ; process 
        else: stop,'civ_calcewn_ndblt: ambiguous doublet'
     endcase 

     ;; Make sure ordered right
     mn = min(abs(strct[ii].wrest[gd]-dblt.wvI),imn,subscript_max=imx)
     subscript[ii,*] = [gd[imn],gd[imx]]

     ;; Error-weighted average
     nI = 10^strct[ii].ncolm[subscript[ii,0]]
     nII = 10^strct[ii].ncolm[subscript[ii,1]]
     if keyword_set(eq_weight) or $
        (strct[ii].signcolm[subscript[ii,0]] eq 0. or $
         strct[ii].signcolm[subscript[ii,1]] eq 0.) then begin
        signI = 1.d
        signII = 1.d
     endif else begin
        signI = strct[ii].signcolm[subscript[ii,0]]*alog(10.)*nI
        signII = strct[ii].signcolm[subscript[ii,1]]*alog(10.)*nII
     endelse 

     ;; Check flags (1 = analyze; 2 = lower limit; 4 = upper limit; 
     ;; 8 = COG value)
     if strct[ii].flg_colm[subscript[ii,0]] eq 1 and $
        strct[ii].flg_colm[subscript[ii,1]] eq 1 then begin

        ;; Take error-weighted value
        wgt = (1./signI^2 + 1./signII^2) 
        ncolm[ii] = (nI/signI^2 + nII/signII^2)/wgt
        signcolm[ii] = sqrt(1./wgt)
        flg_colm[ii] = flg_colm[ii] or 1
     endif else begin        ; both good measurements
        ;; Look for one good measurement
        if strct[ii].flg_colm[subscript[ii,0]] eq 1 or $
           strct[ii].flg_colm[subscript[ii,1]] eq 1 then begin
           ;; Take one _or_ the other
           if strct[ii].flg_colm[subscript[ii,0]] eq 1 then begin
              ncolm[ii] = nI
              signcolm[ii] = signI
              
              ;; Checks
              if (strct[ii].flg_colm[subscript[ii,1]] and 2) eq 2 and $
                 nI lt nII and not keyword_set(silent) then begin
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvII,2),$
                       ' lower limit > ',$
                       strtrim(dblt.wvI,2),' measurement (',$
                       strtrim((nII-nI)/signI,2),' sig)'
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
              endif 
              if (strct[ii].flg_colm[subscript[ii,1]] and 4) eq 4 and $
                 nI gt nII and not keyword_set(silent) then begin
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvII,2),$
                       ' upper limit < ',$
                       strtrim(dblt.wvI,2),' measurement (',$
                       strtrim((nI-nII)/signI,2),' sig)'
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
              endif

           endif else begin
              ncolm[ii] = nII
              signcolm[ii] = signII

              ;; Checks
              if (strct[ii].flg_colm[subscript[ii,0]] and 2) eq 2 and $
                 nII lt nI and not keyword_set(silent) then begin
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvI,2),$
                       ' lower limit > ',$
                       strtrim(dblt.wvII,2),' measurement (',$
                       strtrim((nI-nII)/signII,2),' sig)'
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
              endif 
              if (strct[ii].flg_colm[subscript[ii,0]] and 4) eq 4 and $
                 nII gt nI and not keyword_set(silent) then begin
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvI,2),$
                         ' upper limit <  ',$
                         strtrim(dblt.wvII,2),' measurement (',$
                          strtrim((nII-nI)/signII,2),' sig)'
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
              endif 

           endelse               
           flg_colm[ii] = flg_colm[ii] or 1
        endif else begin        ; one good measurement
           ;; No precise measurement 

           ;; Check for ambiguity (simultaneously lower and upper limit)
           if (strct[ii].flg_colm[subscript[ii,0]] and 6) eq 6 then begin
              if not keyword_set(silent) then $
                 print,'civ_calcewn_ndblt(): upper and lower limit ambiguity: ',$
                       strct[ii].qso,strct[ii].zabs[0],$
                       strct[ii].ion[subscript[ii,0]]
              if strct[ii].signcolm[subscript[ii,0]] lt 1./(3*alog(10.)) then $
                 begin 
                 print,'civ_calcewn_ndblt(): > 3 sigma detection; use upper limit'
                 strct[ii].flg_colm[subscript[ii,0]] = $
                    strct[ii].flg_colm[subscript[ii,0]] - $
                    (strct[ii].flg_colm[subscript[ii,0]] and 2) ; 5
              endif else $
                 stop,'civ_calcewn_ndblt(): < 3 sigma detection; so why lower limit?'
           endif  
           if (strct[ii].flg_colm[subscript[ii,1]] and 6) eq 6 then begin
             if not keyword_set(silent) then $
                 print,'civ_calcewn_ndblt(): upper and lower limit ambiguity: ',$
                       strct[ii].qso,strct[ii].zabs[0],$
                       strct[ii].ion[subscript[ii,1]]
              if strct[ii].signcolm[subscript[ii,1]] lt 1./(3*alog(10.)) then $
                 begin 
                 print,'civ_calcewn_ndblt(): > 3 sigma detection; use upper limit'
                 strct[ii].flg_colm[subscript[ii,1]] = $
                    strct[ii].flg_colm[subscript[ii,1]] - $
                    (strct[ii].flg_colm[subscript[ii,1]] and 2) ; 5
              endif else $
                 stop,'civ_calcewn_ndblt(): < 3 sigma detection; so why lower limit?'
           endif  
           
           if (strct[ii].flg_colm[subscript[ii,0]] and 2) eq 2 and $
              (strct[ii].flg_colm[subscript[ii,1]] and 2) eq 2 then begin
              ;; Both lower limits; take highest
              if nI gt nII then begin
                 ncolm[ii] = nI
                 signcolm[ii] = signI
              endif else begin
                 ncolm[ii] = nII
                 signcolm[ii] = signII
              endelse 
              flg_colm[ii] = flg_colm[ii] or 2
           endif 

           if (strct[ii].flg_colm[subscript[ii,0]] and 4) eq 4 and $
              (strct[ii].flg_colm[subscript[ii,1]] and 4) eq 4 then begin
              ;; Both upper limits; take lowest
              if nI lt nII then begin
                 ncolm[ii] = nI 
                 signcolm[ii] = signI
              endif else begin
                 ncolm[ii] = nII
                 signcolm[ii] = signII
              endelse 
              flg_colm[ii] = flg_colm[ii] or 4
           endif 

           ;; Case where limits may disagree
           if (strct[ii].flg_colm[subscript[ii,0]] and 2) eq 2 and $
              (strct[ii].flg_colm[subscript[ii,1]] and 4) eq 4 then begin
              if not keyword_set(silent) then $
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvI,2),$
                       ' is lower limit; ',strtrim(dblt.wvII,2),$
                       ' is upper limit'
              if nI gt nII then begin
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                          strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
                 ;; Try making bounds more extreme
                 if nI-signI gt nII+signII then begin
                    stop,'civ_calcewn_ndblt(): constraints disagree'
                 endif else begin
                    ncolm[ii] = 0.5*(nI-signI + nII+signII)
                    signcolm[ii] = 0.5*(nII+signII - (nI-signI)) 
                    flg_colm[ii] = flg_colm[ii] or 16 ; new
                    print,'civ_calcewn_ndblt(): constraints disagree so increased allowance'
                 endelse 
              endif else begin
                 ncolm[ii] = 0.5*(nI+nII)
                 signcolm[ii] = 0.5*(nII-nI) 
              endelse 
              flg_colm[ii] = flg_colm[ii] or 6
           endif 

           if (strct[ii].flg_colm[subscript[ii,0]] and 4) eq 4 and $
              (strct[ii].flg_colm[subscript[ii,1]] and 2) eq 2 then begin
              if not keyword_set(silent) then $
                 print,'civ_calcewn_ndblt(): ',strtrim(dblt.wvI,2),$
                       ' is upper limit; ',strtrim(dblt.wvII,2),$
                       ' is lower limit'
              if nI lt nII then begin
                 print,strct[ii].qso,strct[ii].zabs[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,0]],$
                       strct[ii].flg_colm[subscript[ii,1]]
                 ;; Try making bounds more extreme
                 if nI+signI lt nII-signII then begin
                    stop,'civ_calcewn_ndblt(): constraints disagree'
                 endif else begin
                    ncolm[ii] = 0.5*(nI+signI + nII-signII)
                    signcolm[ii] = 0.5*(nI+signI - (nII-signII)) 
                    flg_colm[ii] = flg_colm[ii] or 16 ; new
                    print,'civ_calcewn_ndblt(): constraints disagree so increased allowance'
                 endelse 
              endif else begin
                 ncolm[ii] = 0.5*(nI+nII)
                 signcolm[ii] = 0.5*(nI-nII) 
              endelse 
              flg_colm[ii] = flg_colm[ii] or 6
           endif 
        endelse                 ; neither good measurement
     endelse                    ; both not good measurements
     skip_loop:
  endfor                         ; loop nstrct

  ;; Log 
  if keyword_set(log) then begin
     gd = where(ncolm gt 0.,ngd)
     if ngd ne 0 then begin
        signcolm[gd] = signcolm[gd]/(alog(10.)*ncolm[gd])
        ncolm[gd] = alog10(ncolm[gd])
     endif 
  endif 

  return,ncolm
end                             ; civ_calcewn_ndblt()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_calcewn, strct_fil, instr_list, inspec=inspec, instr_flg=instr_flg, $
                 savfil=savfil,hdf5_fil=hdf5_fil,_extra=extra

  if (N_params() LT 1) then begin 
     print,'Syntax - ' + $
           'civ_calcewn, strct_fil [instr_lst, /savfil] [v1.1]' 
     return
  endif 

  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_calcewn: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent)
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  ;; Determie root directory
  mlss_dir = getenv('MLSS_DIR')
  if strtrim(mlss_dir,2) eq '' then $
     mlss_dir = strtrim(getenv('HOME'),2)+'/MLSS'
  mlss_dir = mlss_dir + '/'

  ;;;;;;;;;;;;;;;;;;;;
  ;; Multiple or single spectra
  ;;;;;;;;;;;;;;;;;;;;
  ;; Reset EW values (or flag to)
  strct.zsig = 0.
  
  if keyword_set(inspec) then begin
     ;; Only one spectrum
     nlist = 1
  endif else begin
     if keyword_set(hdf5_fil) then begin
        ion_indx = where(strct[0].wrest gt 0.)
        unq = uniq(strct.instr[ion_indx[0]],sort(strct.instr[ion_indx[0]]))
        instr = strct[unq].instr[ion_indx[0]]
        nlist = n_elements(unq)
     endif else begin
        ;; Read instrument file list
        if keyword_set(instr_list) then $
           readcol, mlss_dir+instr_list, instr_fil, inst_dw, inst_w0, $
                    format='a,f,f' $
        else $                  ;assume whole structure from one instrument
           readcol, mlss_dir+strtrim(strct[0].instr_fil,2), $
                    instr_fil, inst_dw, inst_w0, $
                    format='a,f,f'
        instr_fil = mlss_dir+instr_fil
        nlist = n_elements(instr_fil)
        if nlist LT 7 then stop,'civ_calcewn: too few spectra'
     endelse                    ; variety of instr_fil
  endelse

  ;; Loop on instruments
  for qq=0L,nlist-1 do begin
     ;; Find lines
     if keyword_set(inspec) then begin
        wave = inspec.wave
        fx = inspec.flux
        sig = inspec.error
        if keyword_set(instr_flg) then begin
           lgv = alog(double(instr_flg))/alog(2)
           nn = fix(lgv+0.00001)
           lin = where(strct.instr MOD 2^(nn+1) GT (2^nn-1), nlin)
        endif else begin
           nn = 1
           lin = where(strct.instr ge nn,nlin)
        endelse 
     endif else begin
        if keyword_set(hdf5_fil) then begin
           fx = hdf5_readspec(getenv('MLSS_DIR')+'/'+strct[0].instr_fil,$
                              instr[qq],wav=wave,sig=sig,npix=npix,$
                              _extra=extra) ; /clean, snr=,/ion_only, etc
           lin = where(strct.instr eq instr[qq],nlin)
           test = where(sig eq 0.,ntest)
           if ntest eq npix then clean = 1 ; flag for later
        endif else begin
           lin = where(strct.instr MOD 2^(qq+1) GT (2^qq-1) and $
                       strct.wrest gt 0., nlin)
           if nlin eq 0 then continue
           test = file_search(instr_fil[qq],count=ntest)
           if ntest eq 0 then begin
              ;; Try reading un-commented file name (now # is embedded)
              prs = strsplit(instr_fil[qq],'#',/extract,count=nprs)
              if nprs gt 1 then test = prs[0] + prs[1] $
              else test = instr_fil[qq]
              
              test2 = file_search(test,count=ntest2)
              if ntest2 eq 0 then stop,'civ_calcewn: no spectra ',instr_fil[qq] $
              else instr_fil[qq] = test2
           endif
           print, 'civ_calcewn: Reading ', instr_fil[qq]

           ;; Open data
           if qq LE 6 then $
              fx = x_readspec(instr_fil[qq],SIG=sig,wav=wave,NPIX=npix,inflg=3)$
           else begin           ; STIS
              spos = strpos(instr_fil[qq], 'f.fits')
              sig_fil = strmid(instr_fil[qq], 0, spos)+'e.fits'
              fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, $
                              fil_sig=sig_fil, inflg=0)
           endelse
           ;; Shift wave
           wave = wave + inst_dw[qq]*inst_w0[qq]/wave
        endelse                 ; normal instrument file
     endelse  

     ;; Sort?
     srt = sort(wave)
     wave = wave[srt]
     fx = fx[srt]
     sig = sig[srt]
     ;; dwv
     dwv = wave - shift(wave,1)
     dwv[0] = dwv[1]

     ;; Trim flux (reduce effect of spurious pixels)
     ;; Un-physical, indicates redux problem or echelle gap
     if keyword_set(clean) then begin
        gd = lindgen(npix)
        ngd = npix
        nbd = 0
     endif else $
        bd = where(sig le 0.,nbd,complement=gd,ncomplement=ngd) 
     if nbd ne 0 then begin
        sig[bd] = 0.
        fx[bd] = -1.
     endif 
     if ngd eq 0 then begin
        for ii=0L,nlin-1 do begin
           nn = array_indices(strct.instr,lin[ii])
           strct[nn[1]].zabs[nn[0]] = 0.5*(strct[nn[1]].wv_lim[nn[0],1]+$
                                           strct[nn[1]].wv_lim[nn[0],0])/$
                                      strct[nn[1]].wrest[nn[0]]-1.
           zsig = -1.
           strct[nn[1]].ew[nn[0]] = 0.
           strct[nn[1]].sigew[nn[0]] = 0.
           strct[nn[1]].ncolm[nn[0]] = 0.
           strct[nn[1]].signcolm[nn[0]] = 0.
        endfor 
        goto,skip_instr
     endif 
     bd = where(fx[gd] lt -abs(sig[gd]),nbd)
     if nbd ne 0 then fx[gd[bd]] = -abs(sig[gd[bd]]) ;floor
     bd = where(fx[gd] gt 1+sig[gd],nbd)
     if nbd ne 0 then fx[gd[bd]] = 1+sig[gd[bd]] ;ceil

     for ii=0L,nlin-1 do begin
        nn = array_indices(strct.instr,lin[ii])
        if nstrct eq 1 then nn = [nn,0]
        ;; Find pixmnx
        mn = min(abs(strct[nn[1]].wv_lim[nn[0],0] - wave), pmin)
        mx = min(abs(strct[nn[1]].wv_lim[nn[0],1] - wave), pmax)

        ;; Calculate EW
        rng = lindgen(pmax-pmin+1) + pmin
        if keyword_set(clean) then begin
           ngd = pmax-pmin+1
           gd = lindgen(ngd) 
           sat = -1
           nsat = 0
        endif else begin
           sat = where(fx[rng] le 0.2*sig[rng] or fx[rng] lt 0.05,nsat)
           gd = where(sig[rng] gt 0.,ngd)
        endelse 

        if ngd eq 0 then begin
           zabs = 0.5*(strct[nn[1]].wv_lim[nn[0],1]+$
                       strct[nn[1]].wv_lim[nn[0],0])/strct[nn[1]].wrest[nn[0]]-1.
           zsig = -1.
           ew = total( (1-fx[rng])*dwv[rng] )
           sigew = -sqrt(total( (sig[rng] * dwv[rng])^2 ))
           ncolm = 0.
           signcolm = 0.
        endif else begin        ;no good data
           rng = rng[gd]
           if keyword_set(clean) then begin
              sat = -1
              nsat = 0
           endif else $
              sat = where(fx[rng] le 0.2*sig[rng] or fx[rng] lt 0.05,nsat)
           ew = total( (1-fx[rng])*dwv[rng] )
           sigew = sqrt(total( (sig[rng] * dwv[rng])^2 ))
           
           ;; Bad data
           if ngd lt 2 then begin
              zabs = 0.5*(strct[nn[1]].wv_lim[nn[0],1]+$
                          strct[nn[1]].wv_lim[nn[0],0])/$
                     strct[nn[1]].wrest[nn[0]]-1.
              zsig = -1.
              ncolm = 0.
              signcolm = 0.
           endif else begin
              ;; Redshift (tau-weight must agree with civ_aodm)
              nwfx = fx[rng]
              if nsat ne 0 then $
                 nwfx[sat] = 0.05 > 0.2*sig[rng[sat]] ;floor
              gd = where(nwfx lt 1.,ngd)       
              if ngd eq 0 then begin
                 gd = lindgen(n_elements(nwfx))
                 nwfx[gd] = 1.          ;still measure error
                 wgt = 1/sig[rng[gd]]^2 ;variance-weight
                 wgtwv = total(wgt*wave[rng[gd]])/total(wgt)
                 zabs = wgtwv / strct[nn[1]].wrest[nn[0]] - 1.
                 zsig = abs(zabs)*1./sqrt(total(wgt))
              endif else begin
                 wgt = alog(1./nwfx[gd])          ;tau-weight
                 wgtwv = total(wgt*wave[rng[gd]])/total(wgt) ;centroid
                 zabs = wgtwv / strct[nn[1]].wrest[nn[0]] - 1.
                 wgtwvvar = [total(((1./nwfx[gd])*wave[rng[gd]]*sig[rng[gd]])^2)/$
                             (total(wgt*wave[rng[gd]]))^2,$
                             total(((1./nwfx[gd])*sig[rng[gd]])^2)/(total(wgt))^2]
                 zsig = abs(1.+zabs)*sqrt(total(wgtwvvar))
              endelse 

              ;; Measure ncolm and error
              velo = 2.998e5*(wave[rng]-wgtwv)/wgtwv
              x_aodm,wave[rng],nwfx,sig[rng],strct[nn[1]].wrest[nn[0]],$
                     ncolm,signcolm,/log,velo=velo ;should be consistent

              ;; Check
              if ncolm lt 0. then begin
                 ;; Especially if due to nwfx = 1. (no fx < 1)
                 ;; x_aodm sets ncolm = -9.99 and signcolm not log
                 ncolm = 0.             ; alog10(signcolm)
                 signcolm = 1./alog(10) ;1-sigma measurement (so flg_colm=5)
              endif 
           endelse              ; lt 2 pixels for measurement (gaps)
        endelse                 ; EW, z, ncolm measured

        if not finite(zabs) then stop,'civ_calcewn: non-finite redshift ',$
                                      strct[nn[1]].wrest[nn[0]],2^qq

        ;; Convert to rest EW (preserving highest significance measurement)
        if strct[nn[1]].zsig[nn[0]] eq 0. then begin ;nothing measured
           strct[nn[1]].zabs[nn[0]] = zabs
           strct[nn[1]].zsig[nn[0]] = zsig
           strct[nn[1]].EW[nn[0]] = ew /(1.+(strct[nn[1]].zabs[nn[0]] > 0.))
           strct[nn[1]].sigEW[nn[0]] = $
              sigew /(1.+(strct[nn[1]].zabs[nn[0]] > 0.))
           strct[nn[1]].ncolm[nn[0]] = ncolm 
           strct[nn[1]].signcolm[nn[0]] = signcolm
        endif else begin
           if ew/sigew gt strct[nn[1]].EW[nn[0]]/strct[nn[1]].sigEW[nn[0]] $
              and sigew gt 0. then begin
              ;; Higher significance
              strct[nn[1]].zabs[nn[0]] = zabs
              strct[nn[1]].zsig[nn[0]] = zsig
              strct[nn[1]].EW[nn[0]] = ew /(1.+(strct[nn[1]].zabs[nn[0]] > 0.))
              strct[nn[1]].sigEW[nn[0]] = $
                 sigew /(1.+(strct[nn[1]].zabs[nn[0]] > 0.))
              strct[nn[1]].instr[nn[0]] = 2^qq
              strct[nn[1]].ncolm[nn[0]] = ncolm
              strct[nn[1]].signcolm[nn[0]] = signcolm
           endif
        endelse                 ; new measurement to evaluate
        
        ;; Analyze flag
        strct[nn[1]].flg_colm[nn[0]] = strct[nn[1]].flg_colm[nn[0]] $
                                       or 1 

        ;; Test 3 sigma measurement and not lower limit (flg=2)
        if signcolm gt 1./(3*alog(10)) $
        then strct[nn[1]].flg_colm[nn[0]] = $
           strct[nn[1]].flg_colm[nn[0]] or 4 $ ;upper limit flag 
        else strct[nn[1]].flg_colm[nn[0]] = $
           strct[nn[1]].flg_colm[nn[0]] and not 4 ; remove
        if nsat gt 0 $
        then strct[nn[1]].flg_colm[nn[0]] = $
           strct[nn[1]].flg_colm[nn[0]] or 2 $ ;lower limit flag 
        else strct[nn[1]].flg_colm[nn[0]] = $
           strct[nn[1]].flg_colm[nn[0]] and not 2 ; remove

        ;; Handle case where flagged lower _and_ upper limit
        ncolm_cog = alog10(ew_to_colm(strct[nn[1]].wrest[nn[0]],$
                                      strct[nn[1]].ew[nn[0]]*1000.,/silent)) ; mA
        ncolm_cog = ncolm_cog[0] ; just to be nice
        ;; Estimate the error
        signcolm_cog = (ew_to_colm(strct[nn[1]].wrest[nn[0]],$
                                   (strct[nn[1]].ew[nn[0]]+strct[nn[1]].sigew[nn[0]])*1000.,/silent) - $
                        ew_to_colm(strct[nn[1]].wrest[nn[0]],$
                                   (strct[nn[1]].ew[nn[0]]-strct[nn[1]].sigew[nn[0]])*1000.,/silent))/$
                       (alog(10.) * 10.^ncolm_cog)
        ncolm_uplim = ncolm     ;+ alog10(1. + alog(10.^2.)*signcolm)
        ;; Try setting upper limit with assuming linear
        ;; portion of COG if ncolm = 0. (from x_aodm)
        if strct[nn[1]].ncolm[nn[0]] eq 0. and $
           strct[nn[1]].signcolm[nn[0]] eq 1./alog(10.) then begin
           strct[nn[1]].flg_colm[nn[0]] = strct[nn[1]].flg_colm[nn[0]] or 12 ; upper limit 4, COG 8
           if not finite(ncolm_cog) then begin
              ;; Using EW as is isn't enough, use 1-sigma EW error
              ncolm_cog = alog10(ew_to_colm(strct[nn[1]].wrest[nn[0]],$
                                            strct[nn[1]].sigew[nn[0]]*1000.,/silent))
              ncolm_cog = ncolm_cog[0] ; just to be nice
              signcolm_cog = ew_to_colm(strct[nn[1]].wrest[nn[0]],$
                                        2*strct[nn[1]].sigew[nn[0]]*1000.,/silent)/$
                             (alog(10.) * 10.^ncolm_cog)          
              if not finite(ncolm_cog) then stop
           endif
           strct[nn[1]].ncolm[nn[0]] = ncolm_cog
           strct[nn[1]].signcolm[nn[0]] = signcolm_cog
        endif 

        if (strct[nn[1]].flg_colm[nn[0]] and 6) ge 6 then begin
           ;; Ambiguous upper- and lower- limit flags
           if strct[nn[1]].ew[nn[0]]/strct[nn[1]].sigew[nn[0]] ge 3. then begin
              ;; >= 3-sigma in EW, so not weak; lower limit b/c saturated
              if ncolm_cog lt ncolm then begin
                 strct[nn[1]].ncolm[nn[0]] = ncolm_cog
                 strct[nn[1]].signcolm[nn[0]] = signcolm_cog ; lower limit
                 strct[nn[1]].flg_colm[nn[0]] = $
                    strct[nn[1]].flg_colm[nn[0]] or 8 ; COG flag
              endif                                   ; else keep AODM
              strct[nn[1]].flg_colm[nn[0]] = strct[nn[1]].flg_colm[nn[0]] and not 4
           endif else begin  ; EW >= 3*sigEW
              ;; < 3-sigma in EW; so weak (make upper limit)
              ;; Pick linear COG column density or 2-sigma upper
              ;; limit AOD
              if ncolm_cog gt ncolm_uplim then begin
                 strct[nn[1]].ncolm[nn[0]] = ncolm_cog
                 strct[nn[1]].signcolm[nn[0]] = signcolm_cog ; upper limit
                 strct[nn[1]].flg_colm[nn[0]] = $
                    strct[nn[1]].flg_colm[nn[0]] or 8 ; COG flag
              endif                                   ; else keep AODM
              strct[nn[1]].flg_colm[nn[0]] = strct[nn[1]].flg_colm[nn[0]] and not 2
           endelse              ; EW < 3*sigEW
        endif else begin        ; upper and lower limit
           ;; Non-conflcting errors
           if  (strct[nn[1]].flg_colm[nn[0]] and 4) ge 4 then begin
              ;; Upper limit
              ;; Pick linear COG column density or 2-sigma upper
              ;; limit AOD
              if ncolm_cog gt ncolm_uplim then begin
                 strct[nn[1]].ncolm[nn[0]] = ncolm_cog
                 strct[nn[1]].signcolm[nn[0]] = signcolm_cog ; upper limit
                 strct[nn[1]].flg_colm[nn[0]] = $
                    strct[nn[1]].flg_colm[nn[0]] or 8 ; COG flag
              endif                                   ; else keep AODM
           endif 
        endelse                 ; limit

        if strct[nn[1]].flg_colm[nn[0]] eq 0 then $
           stop,'civ_calcewn: curiously the line has no flg_colm'

     endfor                     ; loop nlin (in instr)
     skip_instr:
  endfor                        ; loop nlist (instr)


  ;; mA
  strct.EW = strct.EW * 1000
  strct.sigEW = strct.sigEW * 1000

  ;; Overwrite
  if size(strct_fil,/type) eq 7 then begin
     spawn,'cp '+strct_fil+' '+strtrim(strct_fil,2)+'.bkp'
     if keyword_set(savfil) then begin
        civcand = strct
        save,civcand,filename=strct_fil
     endif else mwrfits, strct, strct_fil, /create,/silent
     print,'civ_calcewn: overwrote ',strct_fil
  endif else strct_fil = strct

  return
end
