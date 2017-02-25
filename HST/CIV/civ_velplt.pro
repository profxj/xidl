;+ 
; NAME:
; civ_velplt
;  V1.1
;
; PURPOSE:
;    Velcoity plot of FUSE transitions for a given absorption 
;    system
; CALLING SEQUENCE:
;   
;   civ_velplt, strct_fil, instr_list, NTOT=, CSIZE=,
;   LSIZE=, PSROOT=, XTINT=, ENCAPSULATED=, OUTLINE=, LABEL=
;
; INPUTS:
;
; RETURNS:
;  strct_fil -- FITS file for the FUSE abs lin structure
;  instr_list -- List of instrument files
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NTOT -- Number of plots per page [default: 164]
;  LSIZE -- Label size [default: 1.8]
;  CSIZE -- Numbering character size [default: 1.8]
;  LTHICK -- Line thickness (Default = 1.)
;  XTINT -- xtick interval
;  ENCAPSULATED -- create encapsulated postscript
;  OUTLINE -- indicate region used to measure EW with darker line
;  LABEL -- print zabs upper right corner and wobs 
;  BIN -- bin spectra (simple sum, not S/N-weighted)
;  inspec -- single spectrum structure {wave, flux, error}
;
; OPTIONAL OUTPUTS:
;  PSROOT -- Postscript filename root (append z##.#####) or 
;            if set to number, use QSO name
;  CONCAT -- make only one plot (just append)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   civ_prntcivcand_srtwrest
;   civ_aodm_mtch
;
; REVISION HISTORY:
;   14-Feb-2008  adopted from XIDL fuse_velplt, KLC
;   21-Feb-2008  complete overhaul to use civcandstrct
;    4-Mar-2008  use consistent method for computing AOD agreement
;-
;------------------------------------------------------------------------------
@civ_prntcivcand                ;resolve civ_prntcivcand_srtwrest
@civ_aodm                       ;resolve civ_aodm_mtch

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_velplt, strct_fil, instr_list, NTOT=NTOT, CSIZE=csize, $
                LSIZE=lsize, LTHICK=lthick, ENCAPSULATED=ENCAPSULATED,$
                PSROOT=psroot, XTINT=xtint, OUTLINE=OUTLINE, LABEL=LABEL, $
                HDF5_FIL=hdf5_fil, PLTTITLE=plttitle, $
                BIN=BIN,YMNX=YMNX, CONCAT=CONCAT, FIXVMNX=fixvmnx, $
                PLTERR=plterr, PLT_AOD=plt_aod,$
                VMNX=vmnx, LAPTOP=laptop, inspec=inspec, _EXTRA=extra

  if (N_params() LT 1) then begin 
     print,'Syntax - ' + $
           'civ_velplt, strct_fil, NTOT=, CSIZE=, ' + $
           'LSIZE=, PSROOT=, XTINT= (v1.1)' 
     return
  endif 
  
  if not keyword_set(plt_aod) then plt_aod = 0 ; no

  if not keyword_set( LSIZE ) then lsize = 1.9
  if not keyword_set( CSIZE ) then csize = 3.6
  sv_csize = csize
  if not keyword_set( LTHICK ) then lthick = 2.5
  if not keyword_set( YMNX ) then ymnx = [-0.1,1.48]
  if keyword_set(ntot) then svntot = ntot
  if keyword_set(vmnx) then sv_vmnx = vmnx
  if not keyword_set(xtint) then xtint = 0
  if not keyword_set(xtint) then ytint = 0
  spaces = replicate(' ',30) 

  ;; Determie root directory
  mlss_dir = getenv('MLSS_DIR')
  if strtrim(mlss_dir,2) eq '' then $
     mlss_dir = strtrim(getenv('HOME'),2)+'/MLSS'
  mlss_dir = mlss_dir + '/'

  ;;  Read instrument file list (use only one)
  if keyword_set(instr_list) then begin
     readcol, mlss_dir+instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'
     instr_fil = mlss_dir + instr_fil
     nlist = n_elements(instr_fil)
     if nlist lt 7 then stop,'civ_velplt: too few spectra'
  endif else begin
     if keyword_set(inspec) then nlist = 1 ; only one spectrum
  endelse 

  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_velplt: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent)
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  ;; Select first non-blank ion
  gd = where(strtrim(strct[0].ion,2) ne '')
  prs = strsplit(strct[0].ion[gd[0]],/extract,count=nprs) 
  dblt_name = prs[0]

  ;;;;;;;;;;;;;;
  ;; Begin looping structure
  ;;;;;;;;;;;;;;
  if keyword_set(psroot) and keyword_set(concat) then $
     stop,'civ_velplt: cannot make individual and appended plot'
  if keyword_set(concat) then $
     x_psopen, concat, /portrait,encapsulated=encapsulated $
  else  psarr = strarr(nstrct)

  svinstr_list = ''
  for ss=0,nstrct-1 do begin
     lin = civ_prntcivcand_srtwrest(strct[ss],nlin)

     mx = max(2.998e5*(strct[ss].wv_lim[lin,0]/$
                       (strct[ss].wrest[lin]*$
                        (1+strct[ss].zabs[lin]))-1.),/absolute)
     if keyword_set(sv_vmnx) then vmnx = sv_vmnx else begin
        vmnx = [-abs(mx)-25,abs(mx)+25]
        mx = max(2.998e5*(strct[ss].wv_lim[lin,1]/$
                          (strct[ss].wrest[lin]*$
                           (1+strct[ss].zabs[lin]))-1.),/absolute)
        if abs(mx) gt vmnx[1] then vmnx = [-abs(mx)-25,abs(mx)+25]
        if vmnx[1] gt 549 then vmnx = [-549,549] ;cap
     endelse 

     ;; Read instrument file list (when different)
     
     if not keyword_set(instr_list) and $
        strct[ss].instr_fil ne svinstr_list and $
        not keyword_set(hdf5_fil) and $
        not keyword_set(inspec) then begin 
        readcol, mlss_dir+strtrim(strct[ss].instr_fil,2), $
                 instr_fil, inst_dw, inst_w0, format='a,f,f'
        instr_fil = mlss_dir + instr_fil
        nlist = n_elements(instr_fil)
        if nlist LT 7 then stop,'civ_velplt: too few spectra'
        svinstr_list = strct[ss].instr_fil
     endif 
     
     ;; If only two plot windows, goofy scaling, hard code fix
     if nlin EQ 2 and plt_aod eq 0 then csize = 0.5*sv_csize $
     else csize = sv_csize

     ;; PLOT
     if not keyword_set( svntot ) then begin
        if nlin ge 14L-plt_aod then ntot = 16L-plt_aod else ntot = 14L-plt_aod
        if nlin le ntot/2 then begin
           ;; Effectiveness of limits and one column
           if vmnx[1] le 199. and not keyword_set(vmnx) then $
              vmnx = [-199,199]
        endif else begin
           ;; Effectiveness of limits and 2 columns
           if vmnx[1] le 149. and not keyword_set(vmnx) then $
              vmnx = [-149,149] $
           else begin
              ntot = 2*nlin     ;force one column
              if vmnx[1] le 199. and not keyword_set(vmnx) then $
                 vmnx = [-199,199]
           endelse 
        endelse 
     endif                      ; setting NTOT
     if keyword_set(fixvmnx) then vmnx = fixvmnx

     ;; Save space for AODM profile plot of CIV (why extra 1's floating)
     if nlin+plt_aod LE ntot then begin
        if nlin GT ntot/2 then begin
           npx = 2 
           npy = (nlin+plt_aod)/2 + ((nlin+plt_aod) MOD 2)
        endif else begin
           npx = 1
           npy = nlin+plt_aod
        endelse
     endif else begin
        npx = 2
        npy = ntot/2
     endelse
     xmrg = [6,0.5]
     ymrg = [0,0]

     if keyword_set( LAPTOP ) then begin                     
        csize = 3.9                                          
        xmrg = [6, 0.5]                                      
        lsize=2.9                                            
        lthick = 2.                                          
     endif                                                   

     ;; PSFILE
     if keyword_set( PSROOT ) and not keyword_set(concat) then begin
        ;; psfile = root+'z##.#####'
        if size(psroot,/type) eq 7 then $
           psfile = strtrim(psroot,2)+'z'+$ 
                    strtrim(string(strct[ss].zabs[lin[0]],format='(f8.5)'),2)+$
                    '.ps' $
        else psfile = strtrim(strct[ss].qso,2)+strlowcase(dblt_name)+'z'+$
                      strtrim(string(strct[ss].zabs[lin[0]],format='(f8.5)'),2)+$
                      '.ps'                       
        psarr[ss] = psfile
        x_psopen, psfile, /portrait,encapsulated=encapsulated
     endif 
     if keyword_set( LAPTOP ) then x_tiffset


     ;; starting at p.multi[0], plot p.multi[1] columns and p.multi[2] rows
     ;; with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
     !p.multi=[0,npx,npy,0,1]
     !y.omargin=[3,1.5]
     clr = getcolor(/load)

     ;; LOOP lines
     wrest = 0.d
     svinst = -1     
     for ii=0L,nlin-1 do begin
        if ii EQ npy then !p.multi=[npy,npx,npy,0,1] ;new column

        if keyword_set(inspec) then begin
           wave = inspec.wave
           fx = inspec.flux
           sig = inspec.error

           if keyword_set(bin) then begin
              nwwv = smooth(wave,bin,/edge_truncate)
              nwfx = smooth(fx,bin,/edge_truncate)
              nwsig = smooth(sig,bin,/edge_truncate)
              wave=nwwv
              fx = nwfx
              sig = nwsig
              endif 

        endif else begin           

           ;; Get instrument value
           if keyword_set(hdf5_fil) then instr = strct[ss].instr[lin[ii]] $
           else begin
              lgv = alog(double(strct[ss].instr[lin[ii]])) / alog(2)
              instr = fix(lgv+0.00001)        
           endelse
           
           ;; Read data
           if instr NE svinst then begin
              svinst = instr
              if keyword_set(hdf5_fil) then begin
                 fx = hdf5_readspec(getenv('MLSS_DIR')+'/'+strct[ss].instr_fil,$
                                    instr,wav=wave,sig=sig,npix=npix,$
                                    _extra=extra) ; /clean, snr=,/ion_only, etc
                 test = where(sig eq 0.,ntest)
                 if ntest eq npix then clean = 1 ; flag for later
              endif else begin
                 ;;  Note that STIS starts at instrument 7
                 test = file_search(instr_fil[instr],count=ntest)
                 if ntest eq 0 then begin
                    ;; Try reading un-commented file name (now # is embedded)
                    prs = strsplit(instr_fil[instr],'#',/extract,count=nprs)
                    if nprs gt 1 then test = prs[0] + prs[1] $
                    else test = instr_fil[instr]
                    
                    test2 = file_search(test,count=ntest2)
                    if ntest2 eq 0 then $
                       stop,'civ_velplt: no spectra ',instr_fil[instr] $
                    else instr_fil[instr] = test2
                 endif
                 
                 if instr LE 6 then $
                    fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                                    NPIX=npix, inflg=3)$
                 else begin     ; HST
                    spos = strpos(instr_fil[instr], 'f.fits')
                    sig_fil = strmid(instr_fil[instr], 0, spos)+'e.fits'
                    fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                                    NPIX=npix,fil_sig=sig_fil, inflg=0)
                 endelse
                 ;; Shift wave
                 wave = wave + inst_dw[instr]*inst_w0[instr]/wave
              endelse 
              
              if keyword_set(bin) then begin
                 nwwv = smooth(wave,bin,/edge_truncate)
                 nwfx = smooth(fx,bin,/edge_truncate)
                 nwsig = smooth(sig,bin,/edge_truncate)
                 wave=nwwv
                 fx = nwfx
                 sig = nwsig
              endif 
              
              ;; Sort?
              srt = sort(wave)
              wave = wave[srt]
              fx = fx[srt]
              sig = sig[srt]
           endif 
           
        endelse                 ; inspec=0

        ;; Set vel array (to CIV 1548 zabs)
        x_pixminmax, wave, strct[ss].wrest[lin[ii]], strct[ss].zabs[lin[0]], $
                     vmnx[0], vmnx[1], PIXMIN=pmn, PIXMAX=pmx, VELO=velo

        ;; Tick interval
;        ytint = (ymnx[1]-ymnx[0])/3
;        ytint = round(10*ytint)/10.


        ;; Plot (save last x-axis label for AODM)
        ;; ii ne nlin-1+plt_aod b/c AODM plot last (if set)
        ;; ii ne npy-1 and (ii+plt_aod) ne npx is for first column
        ;; ii ne npy-1 and (ii+plt_aod) mod (npx*npy) ne 0 next will be AODM
        if (ii ne nlin-1+plt_aod) and $
           (ii ne npy-1+plt_aod) then begin 
           plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor, /yminor, $
                 yrange=ymnx, xtickn=spaces, xmargin=xmrg, xtickinterval=xtint,$
                 ymargin=ymrg, /NODATA, charsize=csize, psym=10, $
                 background=clr.white, color=clr.black, $
                 yticks=1,ytickv=[0.5,1.2],$
                 /xstyle, /ystyle, thick=lthick, ytickinterval=ytint
        endif else begin
           plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor, /yminor, $
                 yrange=ymnx, xtitle='Relative Velocity (km s!E-1!N)',$
                 xmargin=xmrg, xtickinterval=xtint,$
                 ymargin=ymrg, /NODATA, charsize=csize, psym=10, $
                 background=clr.white, color=clr.black, $
                 yticks=1,ytickv=[0.5,1.2],$
                 /xstyle, /ystyle, thick=lthick, ytickinterval=ytint
        endelse 
        
        oplot, velo[pmn:pmx], fx[pmn:pmx], psym=10, thick=lthick, $
               color=clr.black
        if keyword_set(plterr) then begin
           ;; Plot error
           oplot,velo[pmn:pmx],sig[pmn:pmx],psym=10,color=clr.red,thick=lthick
        endif                   ; /plter
        
        if keyword_set(outline) then begin
           ;; Indicate range used to measure EW (in 1548 frame)
           vmnbnd = 2.998e5*(strct[ss].wv_lim[lin[ii],0]/$
                             (strct[ss].wrest[lin[ii]]*$
                              (1.+strct[ss].zabs[lin[0]])) - 1.)
           vmxbnd = 2.998e5*(strct[ss].wv_lim[lin[ii],1]/$
                             (strct[ss].wrest[lin[ii]]*$
                              (1.+strct[ss].zabs[lin[0]])) - 1.)
           dum = min(abs(velo-vmnbnd),vmn)
           dum = min(abs(velo-vmxbnd),vmx)
           oplot,velo[vmn:vmx],fx[vmn:vmx],psym=10,color=clr.black,$
                 thick=3*lthick
        endif     
        
        ;; Labels 
        xyouts, 0.02*(vmnx[1]-vmnx[0])+vmnx[0], $ ;ION
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                strtrim(strct[ss].ion[lin[ii]],2), $
                color=clr.black, charsize=LSIZE
        
        ;; zabs in upper right corner; observed wavelength
        if keyword_set(label) then begin
           ;; Indicate saturated
           if (strct[ss].flg_colm[lin[ii]] and 2) ge 2 then $
              xyouts,'!ES!N',color=clr.red,charsize=LSIZE
           ;; Indicated lower limit
           if (strct[ss].flg_colm[lin[ii]] and 4) ge 4 then $
              xyouts,'!EW!N',color=clr.red,charsize=LSIZE
        endif 

        ;; Lines
        oplot, [0., 0.], ymnx, color=clr.black, linestyle=2, $
               thick = 2.5*lthick
        oplot, [-10000., 10000.], [0.,0.], color=clr.blue, linestyle=3, $
               thick = 2.5*lthick
        oplot, [-10000., 10000.], [1.,1.], color=clr.limegreen, linestyle=3, $
               thick = 2.5*lthick
     endfor                     ;loop nlin
     
     ;; Axis label 
     if plt_aod eq 1 then $
        xyouts, 0.045, 0.67, 'Normalized Flux', $
                alignment=0.5, ORIENTATION=90., /normal, charsize=lsize $
     else $
        xyouts, 0.045, 0.5, 'Normalized Flux', $
                alignment=0.5, ORIENTATION=90., /normal, charsize=lsize

     lscl = 0.7
     if keyword_set(plttitle) then begin
        xyouts,0.56,0.97,plttitle[ss],charsize=lscl*LSIZE,alignment=0.5,/normal
     endif else begin
        ;; Lya Forest tag
        if (strct[ss].flg_sys[0] and 8) eq 0 then $
;           xyouts,0.5*(0.88+0.25),0.97,'Ly!9a!X Forest',$
           xyouts,0.6,0.97,'Ly!9a!X Forest',$
                  color=clr.red,charsize=lscl*LSIZE, alignment=0.5,/normal
        
        ;; Print QSO (RA) in upper left corner
        x_radec,ra_str,dec_str,strct[ss].ra,strct[ss].dec,/flip
        dum = strtrim(strct[ss].qso,2)+' ('+strtrim(ra_str,2)+')'
;        xyouts,0.12,0.97,dum,charsize=lscl*LSIZE,alignment=0,/normal
        xyouts,0.01,0.97,dum,charsize=lscl*LSIZE,alignment=0,/normal
        
        ;; Reshift label upper right corner
        wr = strtrim(strmid(strct[ss].wrest[0],0,strpos(strct[ss].wrest[0],'.',/reverse_search)),2)
        xyouts,0.99,0.97,'!8z!X!D'+wr+'!N = '+$
               strtrim(string(strct[ss].zabs[lin[0]],format='(f8.5)'),2),$
               charsize=lscl*LSIZE,alignment=1.0,/normal
     endelse                     ; /label

    
     ;;;;;;;;;;;;;;
     ;; AODM PLOTS
     ;;;;;;;;;;;;;;
     if plt_aod eq 1 then begin
        ii = nlin
        if ii EQ npy then !p.multi=[npy,npx,npy,0,1] ;next column

        gd = where(strct[ss].aodm_civ[*,1] gt 0.,ngd) ;CIV 1550
        if ngd eq 0 then begin
           mn = 0.
           mx = 1.
        endif else begin
           mx = max(strct[ss].aodm_civ[gd,1])
           mn = min(strct[ss].aodm_civ[gd,1])
        endelse 
        gd = where(strct[ss].aodm_civ[*,0] gt 0.,ngd) ;CIV 1548
        if ngd ne 0 then begin
           mn = min(strct[ss].aodm_civ[gd,0]) < mn
           mx = max(strct[ss].aodm_civ[gd,0]) > mx
        endif 
        aodcolmnx = [mn > 11.1,mx+0.3]

        ;; Tick interval
        ytint = (aodcolmnx[1]-aodcolmnx[0])/2
        ytint = round(ytint)

        ;; Plot
;        ytickv = lindgen(10)*2  ; only want to show a few ticks
;        gd = where(ytickv ge aodcolmnx[0] and ytickv le aodcolmnx[1],ngd)
        getion, strct[ss].wrest[lin[0]], ion, elm
        dum = 'log !8N!X('+elm+'!E+'+strtrim(ion-1,2)+'!N)'
        plot, vmnx, aodcolmnx, xrange=vmnx, $
              yrange=aodcolmnx, xmargin=xmrg, ymargin=ymrg, /NODATA, $
              charsize=csize, psym=10, background=clr.white, $
              color=clr.black, ytickinterval=ytint, $
              xstyle=1, /ystyle, thick=lthick, xtickinterval=xtint, $
;              yticks = ngd-1, ytickv = ytickv[gd], $
              yticks = 1, ytickv = [12,13], /xminor, $
              xtitle='Relative Velocity (km s!E-1!N)', $ 
              linestyle=lsty,ytitle=dum
        oplot,strct[ss].aodm_vel[*,0],strct[ss].aodm_civ[*,0],$
              psym=10,color=clr.black,thick=lthick
        ;; Offset of vcent(1550) from vcent(1548)
        dv = 2.998e5*(strct[ss].zabs[lin[1]]-strct[ss].zabs[lin[0]])/$
             (1.+strct[ss].zabs[lin[0]])
        oplot,strct[ss].aodm_vel[*,1]+dv,strct[ss].aodm_civ[*,1],$
              psym=10,color=clr.red,thick=lthick,linestyle=1 ; dotted

        ;; CIV 1548
        vmnxciv = 2.998e5*(strct[ss].wv_lim[lin[0],*]/$
                           (strct[ss].wrest[lin[0]]*$
                            (1+strct[ss].zabs[lin[0]]))-1.)
        mn = min(strct[ss].aodm_vel[*,0]-vmnxciv[0],imn,/absolute)
        mx = min(strct[ss].aodm_vel[*,0]-vmnxciv[1],imx,/absolute)
        if keyword_set(outline) then $
           oplot,strct[ss].aodm_vel[imn:imx,0],$
                 strct[ss].aodm_civ[imn:imx,0],psym=10,color=clr.black,$
                 thick=3*lthick

        ;; CIV 1550
        vlim = 2.998e5*(strct[ss].wv_lim[lin[1],*]/$
                        (strct[ss].wrest[lin[1]]*$
                         (1+strct[ss].zabs[lin[1]]))-1.)
        mn = min(strct[ss].aodm_vel[*,1]-vlim[0],imn,/absolute)
        mx = min(strct[ss].aodm_vel[*,1]-vlim[1],imx,/absolute)
        if keyword_set(outline) then $
           oplot,strct[ss].aodm_vel[imn:imx,1]+dv,$
                 strct[ss].aodm_civ[imn:imx,1],psym=10,color=clr.red,$
                 thick=3*lthick,linestyle=1 ; dotted
        vmnxciv = [vmnxciv[0] < vlim[0], vmnxciv[1] > vlim[1]]
        
        ;; Labels
        if keyword_set(label) then begin
           flg_mtch = civ_aodm_mtch(strct[ss])
           xyouts,0.05*(vmnx[1]-vmnx[0])+vmnx[0],$
                  0.85*(aodcolmnx[1]-aodcolmnx[0])+aodcolmnx[0],$
                  strtrim(round(flg_mtch*100),2)+'% AOD agree',$
                  color=clr.black,charsize=0.5*lsize
        endif

        ;; Lines
        oplot, [0., 0.], [0.,20.], color=clr.black, linestyle=2 ,$ ;vertical
               thick=2.5*lthick
     endif 
     if keyword_set( PSROOT ) and not keyword_set(concat) then begin
        x_psclose
        !p.multi=[0,1,1]
        !y.omargin=[0,0]
     endif
     if keyword_set( LAPTOP ) then x_tiffset, /unset
  endfor                        ;loop nstrct

  if keyword_set(concat) then begin
     x_psclose
     !p.multi=[0,1,1]
     !y.omargin=[0,0]
  endif

  if keyword_set(psroot) and not keyword_set(concat) then begin
     if n_elements(psarr) gt 1 then begin
        if size(psroot,/type) eq 7 then $
           psfile = strtrim(psroot,2)+'_all.ps' $
        else psfile = 'civ_velplt_all.ps'
        psconcat,psfile,psarr
        print,'civ_velplt: created ',psfile
     endif 
  endif

  return
end
