;+ 
; NAME:
; fuse_velplt
;  V1.1
;
; PURPOSE:
;    Velcoity plot of FUSE transitions for a given absorption 
;    system
; CALLING SEQUENCE:
;   
;   fuse_velplt, strct_fil, instr_list, vel_fil, NTOT=, CSIZE=,
;   LSIZE=, PSFILE=, XTINT=, ENCAPSULATED=, MULTI=, OUTLINE=,
;   LABEL=, MSGFIL=, VOIGT=
;
; INPUTS:
;
; RETURNS:
;  strct_fil -- FITS file for the FUSE abs lin structure
;  instr_list -- List of instrument files
;  vel_fil -- Input file for velocity plot
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NTOT -- Number of plots per page [default: 16]
;  LSIZE -- Label size [default: 1.8]
;  CSIZE -- Numbering character size [default: 1.8]
;  LTHICK -- Line thickness (Default = 1.)
;  XTINT -- xtick interval
;  ENCAPSULATED -- create encapsulated postscript
;  MULTI -- overplot multiple detections
;  OUTLINE -- indicate region used to measure EW with darker line
;  LABEL -- print zabs upper right corner and wobs 
;  VOIGT -- [ncolm, dopb] or [ion, ncolm, dopb] structure to overplot
;           Voigt profile 
;  VZABS -- redshift to center voigt profiles (default is input z + float)
;  BIN -- bin spectra (simple sum, not S/N-weighted)
;  PLTXTR -- array of wrest, instr to plot extra (doesn't work for
;            multi)
;  PLTERR -- plot error array (changes colors)
;  ROOTDIR - pre-pend string to instr_list file names
;  SUFFIX - for STIS-like files, suffix to search and replace
;           (default: ['f.fits','e.fits'])
;;
; OPTIONAL OUTPUTS:
;  PSFILE -- Postscript filename
;  MSGFIL -- filename to receive messages from MULTI function
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2003 Written by JXP
;    1-Apr-2005 MULTI capability added by KLC
;   28-Apr-2005 Added zabs in upper right-hand corner of MULTI plots
;   13-Feb-2006 Added /outline option, KLC
;    7-Mar-2006 Handle H2 features, KLC
;    7-Jun-2006 Improve handling of combined instrument flags, add /label
;   21-Sep-2006 Update /outline,/multi combo
;   16-Oct-2006 added VOIGT, KLC
;    8-Dec-2006 added BIN, KLC
;-
;------------------------------------------------------------------------------
pro fuse_velplt, strct_fil, instr_list, vel_fil, NTOT=NTOT, CSIZE=csize, $
                 LSIZE=lsize, LTHICK=lthick, ENCAPSULATED=ENCAPSULATED,$
                 PSFILE=psfile, XTINT=xtint,LAPTOP=laptop, $
                 MULTI=multi,MSGFIL=msgfil,OUTLINE=OUTLINE, LABEL=LABEL, $
                 VOIGT=voigt, VZABS=vzabs, BIN=BIN, PLTXTR=PLTXTR, PLTERR=PLTERR, $
                 ROOTDIR=ROOTDIR, SUFFIX=SUFFIX, _EXTRA=extra

  ;; lowzovi_prsdat -- Reads in DLA data to a structure
  
  if (N_params() LT 3) then begin 
     print,'Syntax - ' + $
           'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
           'LSIZE=, PSFILE=, XTINT= (v1.1)' 
     return
  endif 
  
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.6
  lscl = 0.7                    ; label scale
  if not keyword_set( CSIZE ) then csize = 3.2
  sv_csize = csize
  if not keyword_set( LTHICK ) then lthick = 2.5
;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'
  if not keyword_set( SUFFIX ) then suffix = ['f.fits','e.fits']
  if not keyword_set(xtint) then xtint = 0
  if not keyword_set(xtint) then ytint = 0

  ;;  Read instrument file list
  readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'
  if keyword_set(rootdir) then instr_fil=rootdir+instr_fil

  nlist = n_elements(instr_fil)
  if nlist LT 7 then stop

  ;;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)

  ;; Open vel_fil
  close, /all
  openr, 11, vel_fil
  readf, 11, zabs, FORMAT='(f12.7)'
  readf, 11, vmin,vmax
  vmnx = [vmin,vmax]
  nlin = 0
  readf, 11, nlin

  if keyword_set(pltxtr) then begin
     npltxtr = n_elements(pltxtr[0,*])
     nlin = nlin+npltxtr
  endif else npltxtr = 0

  ;; PLOT
  xmrg = [6,0.5]
  ymrg = [0,0]
  if nlin LE ntot then begin
     if nlin GT 8 then begin
        npx = 2 
        npy = nlin/2 + (nlin MOD 2)
     endif else begin
        npx = 1
        npy = nlin
;        LSIZE = 2.2
     endelse
  endif else begin
     npx = 2
     npy = ntot/2
  endelse

  ;; Funny busines s with font size when only two
  if nlin eq 2 then csize = 0.5*sv_csize $
  else csize = sv_csize 


  if keyword_set( LAPTOP ) then begin
     csize = 3.9
     xmrg = [6, 0.5]
     lsize=2.9
     lthick = 2.
  endif

  ;; PSFILE
  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait,$
                                          encapsulated=encapsulated
  if keyword_set( LAPTOP ) then x_tiffset


  ;;starting at p.multi[0], plot p.multi[1] columns and p.multi[2] rows
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[0,npx,npy,0,1]
  !y.omargin=[3,1.5]
  clr = getcolor(/load)
  ;;color of spectra in order as listed in instrument file
                                ;SiC 1B, SiC 2A, LiF 2B, LiF 1A, SiC
                                ;1A, LiF 1B, LiF 2A, STIS 
  clrspec = [clr.brown,clr.steelblue,clr.maroon,clr.yellow,clr.purple,$
             clr.cyan,clr.magenta,clr.navy]

  ;; MSGFIL
  if keyword_set( MSGFIL ) then begin 
     openw,17,msgfil
     printf,17,zabs
     printf,17,'Color Legend'
     printf,17,'SiC 1B   Brown'
     printf,17,'SiC 2A   Steel Blue'
     printf,17,'LiF 2B   Maroon'
     printf,17,'LiF 1A   Yellow'
     printf,17,'SiC 1A   Purple'
     printf,17,'LiF 1B   Cyan'
     printf,17,'LiF 2A   Magenta'
     printf,17,'STIS     Navy'
  endif 


  ;; Vel strct


  wrest = 0.d
  svinst = -1         ;because SiC 1B is binary 1, which makes instr list # = 0
  nblnd = 0


  ;; LOOP
  if not keyword_set(MULTI) then begin
     for ii=0L,nlin-npltxtr-1 do begin
        if ii EQ npy then !p.multi=[npy,npx,npy,0,1]
        readf, 11, wrest, ymin, ymax, nblnd
        ymnx = [ymin,ymax]

        ;; Get instrument value
        a = where(abs(strct.zabs-zabs) LT 0.002 AND $
                  abs(strct.wrest-wrest) LT 0.003, na)
        if na EQ 0 then begin
           stop, 'fuse_velplt: No wave! ', wrest
        endif else begin
           if na eq 1 then a = a[0] $
           else begin
              ;;Assume closest in redshift is best
              dum = min(abs(strct[a].zabs-zabs),dd)
              a = a[dd]
           endelse 
        endelse
        lgv = alog(double(strct[a].instr)) / alog(2)
        instr = fix(lgv+0.00001)        
;      if lgv - jj LT 0.00001 then instr = jj else

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Added by KLC, 7 Jun 06
        ;; LiF 1A (8) better than SiC 1A (16)
        ;; LIF 2A (64) or LiF 1B (32) better than STIS (128)
        nn = strct[a].instr - 2^instr
        case 2^instr of
           16: if nn ne 0 then instr = fix(alog(double(nn)+0.00001)/alog(2))
           128: if nn ne 0 then $
              instr = fix(alog(double(nn)+0.00001)/alog(2))
           else:                ;should be good
        endcase
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        ;; Read data
        if instr NE svinst then begin
           svinst = instr
           ;;  Note that STIS starts at instrument 7
           if instr LE 6 then $
              fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                              NPIX=npix, inflg=3)$
           else begin           ; STIS
              spos = strpos(instr_fil[instr], suffix[0])
              sig_fil = strmid(instr_fil[instr], 0, spos)+suffix[1]
              fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                              NPIX=npix,fil_sig=sig_fil, inflg=0)
           endelse
           ;; Shift wave
           wave = wave + inst_dw[instr]*inst_w0[instr]/wave

           if keyword_set(bin) then begin
              nwpx = npix/bin + 1
              nwwv = fltarr(nwpx)
              nwfx = fltarr(nwpx)
              nwsig = fltarr(nwpx)
              rng = lonarr(2)
              for mm=1L,nwpx-1 do begin
                 rng = [(mm-1)*bin,((mm*bin-1) < (npix-1))]
                 nwwv[mm-1] = total(wave[rng[0]:rng[1]])/bin
                 nwfx[mm-1] = total(fx[rng[0]:rng[1]])/bin
                 nwsig[mm-1] = total(sig[rng[0]:rng[1]])/bin
                                ;if mm eq nwpx - 2 then stop
              endfor 

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

        ;; Set vel array
        x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                     PIXMAX=pmx, VELO=velo

        ;; Y-axis
;      if npx EQ 2 AND (ii+1) MOD ntot GT npy then flg_yax = 1 else flg_yax=0

;        if ymnx[1]-ymnx[0] GT 0.6 then ytint = 0.4 else ytint = 0

        ;; Plot
;      if ii GT npy-1 then ysty=5 else ysty=1
;        ysty=1
        if keyword_set(outline) then begin
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           ;; Indicate range used to measure EW, added by KLC
           if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
              ((ii+1) MOD (npx*npy) NE 0) AND $
              (((ii+1) MOD npy) NE 0) then begin
              spaces = replicate(' ',30)
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xtickn=spaces, xmargin=xmrg,xtickinterval=xtint,$
                    ymargin=ymrg, /NODATA, charsize=csize, psym=10, $
                    background=clr.white, color=clr.black, $
                    yticks=1, ytickv=[0.5,1.2], $
                    /xstyle, /ystyle, thick=lthick, ytickinterval=ytint
           endif else begin
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xmargin=xmrg, ymargin=ymrg, /NODATA, $
                    charsize=csize, psym=10, background=clr.white, $
                    color=clr.black, $
                    yticks=1, ytickv=[0.5,1.2], $
                    /xstyle, /ystyle, thick=lthick, xtickinterval=xtint, $
                    xtitle='Relative Velocity (km s!E-1!N)', ytickinterval=ytint
           endelse
           oplot, velo[pmn:pmx], fx[pmn:pmx], psym=10, thick=lthick, $
                  color=clr.black
           
           vmnbnd = 2.998e5*(strct[a].wv_lim[0]/$
                             (strct[a].wrest*(1.+zabs)) - 1.)
           vmxbnd = 2.998e5*(strct[a].wv_lim[1]/$
                             (strct[a].wrest*(1.+zabs)) - 1.)
           dum = min(abs(velo-vmnbnd),vmn)
           dum = min(abs(velo-vmxbnd),vmx)
           oplot,velo[vmn:vmx],fx[vmn:vmx],psym=10,color=clr.black,$
                 thick=3*lthick
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        endif else begin 

           if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
              ((ii+1) MOD (npx*npy) NE 0) AND $
              (((ii+1) MOD npy) NE 0) then begin
              spaces = replicate(' ',30)
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xtickn=spaces, xmargin=xmrg, $
                    ymargin=ymrg, NODATA=nblnd, charsize=csize, psym=10, $
                    background=clr.white, color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2], xtickinterval=xtint,$
                    /xstyle, /ystyle, thick=lthick, ytickinterval=ytint
           endif else begin
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xmargin=xmrg, ymargin=ymrg, NODATA=nblnd, $
                    charsize=csize, psym=10, background=clr.white, $
                    color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2],$
                    /xstyle, /ystyle, thick=lthick, xtickinterval=xtint, $
                    xtitle='Relative Velocity (km s!E-1!N)', ytickinterval=ytint
           endelse
        endelse 
        if keyword_set(plterr) then $
           oplot,velo[pmn:pmx],sig[pmn:pmx],psym=10,color=clr.red,thick=lthick

        
;      if ii GT npy-1 then begin
;          spaces = replicate(' ',30)
;          axis, yrange=ymnx, ystyle=1, yaxis=1, charsize=csize, $
;            ytickinterval=ytint
;          axis, yrange=ymnx, ystyle=1, yaxis=0, ytickn=spaces, $
;            ytickinterval=ytint
;      endif

        ;; BLENDS
        nlow = pmn
        for jj=0L,nblnd-1 do begin
           readf, 11, vmin, vmax, FORMAT='(f,f)'
           mn = min(abs(velo-vmin),pixmin)
           mn = min(abs(velo-vmax),pixmax)

           ;; Plot good
           nlow = nlow < pixmin
           oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
                  color=clr.black, psym=10, thick=lthick
           ;; Plot blend
           oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
                  psym=10, linestyle=2, thick=1.
           ;; End
           if jj EQ nblnd-1 then begin
              if pixmax LT pmx then $
                 oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                        psym=10, thick=lthick
           endif
           ;; Set nlow
           nlow = pixmax
        endfor
        

        ;; Labels
        getfnam, wrest, fv, nam, _extra=extra ;eg fil=getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'

        ;;Label H2 features
        h2 = where(strcmp(strtrim(nam,2),''),nh2)
        if nh2 ne 0 then $
           nam[h2] = 'H!D2!N '+string(wrest[h2],format='(f6.1)')

        xyouts, 0.02*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                strtrim(nam,2), color=clr.black, charsize=LSIZE
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; zabs in upper right corner; observed wavelength
        ;; added by KLC
        if keyword_set(label) then begin
           wobs = wrest*(1.+zabs)
           wobsnam = ' ('+strmid(strtrim(wobs,2),0,7)+',' +$
                     strtrim(2^instr,2)+')'
           ;; append to most recent positions
           xyouts, wobsnam, color=clr.black,charsize=lscl*LSIZE   
           ;; Reshift label upper right corner
           xyouts,0.99,0.97,'!8z!X!Dabs!N = '+strtrim(zabs,2),$
                  charsize=lscl*LSIZE,alignment=1.0,/normal
           ;; Print input file name in upper left corner
           dum = [strpos(vel_fil,'/',/reverse_search)+1,$
                  strpos(vel_fil,'.',/reverse_search)]
           xyouts,0.01,0.97,strmid(vel_fil,dum[0],dum[1]-dum[0]),$
                  charsize=lscl*LSIZE,alignment=0,/normal
        endif                   ; label
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        if (ii+1) MOD ntot EQ 1 then begin
;           if npx EQ 1 then $
;              xyouts, 0.045, 0.5, 'Normalized Flux', $
;                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize $
;           else $
              xyouts, 0.045, 0.5, 'Normalized Flux', $
                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize
        endif
        
        ;; Lines
        if keyword_set(plterr) then begin
           oplot, [0., 0.], ymnx, color=clr.black, linestyle=2
           oplot, [-10000., 10000.], [0.,0.], color=clr.blue, linestyle=3
        endif else begin
           oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2
           oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3
        endelse 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Overplot Voigt profile with given N and b
        ;; added by KLC 16 Oct 2006
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if keyword_set(voigt) then begin
           vplin = x_setline(wrest)
           vplin.zabs = strct[a].zabs
           if keyword_set(vzabs) then vplin.zabs = vzabs else $
              ;; 30 km s!E-1!N wiggle room
              if abs(strct[a].zabs-zabs) gt 0.0001 then vplin.zabs = zabs
           
           sz = size(voigt,/dimensions)
           if size(voigt,/n_dimensions) eq 1 then begin ;one element to overlay VP
              case sz[0] of
                 2: begin       ;default: only HI lines
                    if stregex(nam,'HI',/fold_case,/boolean) then begin
                       vplin.n = float(voigt[0])
                       vplin.b = float(voigt[1])
                    endif else vplin.n = -1
                 end
                 3: begin       ;only specified element
                    if stregex(nam,voigt[0],/fold_case,/boolean) then begin
                       vplin.n = float(voigt[1])
                       vplin.b = float(voigt[2])
                    endif else vplin.n = -1
                 end
                 else: stop,'fuse_velplt: voigt wrong size'
              endcase 
           endif else begin     ;multiple elements to overlay VP
              for vv=0,sz[1]-1 do begin
                 if stregex(nam,voigt[0,vv],/fold_case,/boolean) then begin
                    vplin.n = float(voigt[1,vv])
                    vplin.b = float(voigt[2,vv])
                    break
                 endif else vplin.n = -1
              endfor 
           endelse 

           if instr le 6 then fwhm = 5 $ ;FUSE
           else fwhm = 2                 ;STIS
           
           if vplin.n ne -1 then $ ;check current line correct
              vptau = x_voigt(wave[pmn:pmx], vplin, fwhm=fwhm) $
           else vptau = replicate(1.,pmx-pmn+1)
           
           oplot, velo[pmn:pmx],vptau,color=clr.limegreen, linestyle=3
        endif else begin
           oplot, [-10000., 10000.], [1.,1.], color=clr.limegreen, linestyle=3
        endelse 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     endfor 

     ;; EXTRA PLOTS
     for ii=nlin-npltxtr,nlin-1 do begin
        if ii EQ npy then !p.multi=[npy,npx,npy,0,1]
        wrest = pltxtr[0,ii-(nlin-npltxtr)]
        
        a = long(pltxtr[1,ii-(nlin-npltxtr)])
        lgv = alog(a) / alog(2)
        instr = fix(lgv+0.00001)        
        
        ;; Read data
        if instr NE svinst then begin
           svinst = instr
           ;;  Note that STIS starts at instrument 7
           if instr LE 6 then $
              fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                              NPIX=npix, inflg=3)$
           else begin           ; STIS
              spos = strpos(instr_fil[instr], suffix[0])
              sig_fil = strmid(instr_fil[instr], 0, spos)+suffix[1]
              fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                              NPIX=npix,fil_sig=sig_fil, inflg=0)
           endelse
           ;; Shift wave
           wave = wave + inst_dw[instr]*inst_w0[instr]/wave

           if keyword_set(bin) then begin
              nwpx = npix/bin + 1
              nwwv = fltarr(nwpx)
              nwfx = fltarr(nwpx)
              nwsig = fltarr(nwpx)
              rng = lonarr(2)
              for mm=1L,nwpx-1 do begin
                 rng = [(mm-1)*bin,((mm*bin-1) < (npix-1))]
                 nwwv[mm-1] = total(wave[rng[0]:rng[1]])/bin
                 nwfx[mm-1] = total(fx[rng[0]:rng[1]])/bin
                 nwsig[mm-1] = total(sig[rng[0]:rng[1]])/bin
              endfor 
              
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
        
        ;; Set vel array
        x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                     PIXMAX=pmx, VELO=velo
        
        if ymnx[1]-ymnx[0] GT 0.6 then ytint = 0.4 else ytint = 0
        
        ;; Plot
        ysty=1
        if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
           ((ii+1) MOD (npx*npy) NE 0) $
           (((ii+1) MOD npy) NE 0) then begin
           spaces = replicate(' ',30)
           plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                 yrange=ymnx, xtickn=spaces, xmargin=xmrg, $
                 ymargin=ymrg, NODATA=nblnd, charsize=csize, psym=10, $
                 yticks=1,ytickv=[0.5,1.2],xtickinterval=xtint,$
                 background=clr.white, color=clr.black, $
                 /xstyle, /ystyle, thick=lthick, ytickinterval=ytint, $
                 linestyle=lsty
        endif else begin
           plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                 yrange=ymnx, xmargin=xmrg, ymargin=ymrg, NODATA=nblnd, $
                 charsize=csize, psym=10, background=clr.white, $
                 color=clr.black, $
                 yticks=1,ytickv=[0.5,1.2],$
                 /xstyle, /ystyle, thick=lthick, xtickinterval=xtint, $
                 xtitle='Relative Velocity (km s!E-1!N)', ytickinterval=ytint, $
                 linestyle=lsty
        endelse
        if keyword_set(plterr) then $
           oplot,velo[pmn:pmx],sig[pmn:pmx],color=clr.red,thick=lthick
        
        ;; Labels
        getfnam, wrest, fv, nam, _extra=extra ;eg fil=getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
        
        ;;Label H2 features
        h2 = where(strcmp(strtrim(nam,2),''),nh2)
        if nh2 ne 0 then $
           nam[h2] = 'H!D2!N '+string(wrest[h2],format='(f6.1)')
        
        xyouts, 0.02*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                strtrim(nam,2), color=clr.black, charsize=LSIZE
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; zabs in upper right corner; observed wavelength
        ;; added by KLC
        if keyword_set(label) then begin
           wobs = wrest*(1.+zabs)
           wobsnam = ' ('+strmid(strtrim(wobs,2),0,7)+','+$
                     strtrim(2^instr,2)+')'
           ;; append to most recent positions
           xyouts, wobsnam, color=clr.black,charsize=0.5*LSIZE        
;            ;; Reshift label upper right corner
;            xyouts,0.9,0.98,'z!Dabs!N = '+strtrim(zabs,2),$
;              charsize=0.5*LSIZE,alignment=0.5,/normal
;            ;; Print input file name in upper left corner
;            dum = [strpos(vel_fil,'/',/reverse_search)+1,$
;                   strpos(vel_fil,'.',/reverse_search)]
;            xyouts,0.25,0.98,strmid(vel_fil,dum[0],dum[1]-dum[0]),$
;              charsize=0.5*LSIZE,alignment=0.5,/normal
        endif 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        
        if (ii+1) MOD ntot EQ 1 then begin
;           if npx EQ 1 then $
;              xyouts, 0.045, 0.5, 'Normalized Flux', $
;                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize $
;           else $
              xyouts, 0.045, 0.5, 'Normalized Flux', $
                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize
        endif
        
        ;; Lines
        oplot, [-10000., 10000.], [1.,1.], color=clr.limegreen, linestyle=3,$
               lthick=2.5*lthick
        if keyword_set(plterr) then begin
           oplot, [0., 0.], ymnx, color=clr.black, linestyle=2, $
                  lthick=2.5*lthick
           oplot, [-10000., 10000.], [0.,0.], color=clr.blue, linestyle=3, $
                  lthick=2.5*lthick
        endif else begin
           oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, $
                  lthick=2.5*lthick
           oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3, $
                  lthick=2.5*lthick
        endelse 
     endfor  ;;end extra



  endif else begin 
     ;;MULTI keyword to overplot features from different channels
     ;;added by KLC
     for ii=0L,nlin-1 do begin
        if ii EQ npy then !p.multi=[npy,npx,npy,0,1]
        readf, 11, wrest, ymin, ymax, nblnd
        ymnx = [ymin,ymax]

        ;; Get instrument value
        a = where(abs(strct.zabs-zabs) LT 0.001 AND $
                  abs(strct.wrest-wrest) LT 0.003, na)
        if na EQ 0 then begin
           stop, 'fuse_velplt: No wave! ', wrest
        endif else begin
           if na eq 1 then a = a[0] $
           else begin
              ;;Assume closest in redshift is best
              dum = min(abs(strct[a].zabs-zabs),dd)
              a = a[dd]
           endelse 
        endelse

        ;;converts binary instrument number to
        ;;correspond with number in instrument list
        lgv = alog(double(strct[a].instr)) / alog(2)
        instr = fix(lgv+0.00001)

        ;;Dissect instrument flag to constituent flags
        tmp = strct[a].instr - 2^instr
        if tmp ne 0 then done = 0 else done = 1
        while not done do begin
           nn = alog(double(tmp))/alog(2)
           instr = [instr,fix(nn+0.00001)]
           tmp = tmp - 2^instr[n_elements(instr)-1]
           if tmp eq 0 then done = 1
        endwhile 
        na = n_elements(instr)

        ;; Plot first spectra
        ;; Read data only if need new instrument spectrum
        if instr[0] NE svinst then begin
           svinst = instr[0]
           ;;  Note that STIS starts at instrument 7
           if instr[0] LE 6 then $
              fx = x_readspec(instr_fil[instr[0]], SIG=sig, wav=wave, $
                              NPIX=npix, inflg=3)$
           else begin           ; STIS
              spos = strpos(instr_fil[instr[0]], suffix[0])
              sig_fil = strmid(instr_fil[instr[0]], 0, spos)+suffix[1]
              fx = x_readspec(instr_fil[instr[0]], SIG=sig, wav=wave, $
                              NPIX=npix,fil_sig=sig_fil, inflg=0)
           endelse
           ;; Shift wave
           wave = wave + inst_dw[instr[0]]*inst_w0[instr[0]]/wave

           if keyword_set(bin) then begin
              nwpx = npix/bin + 1
              nwwv = fltarr(nwpx)
              nwfx = fltarr(nwpx)
              nwsig = fltarr(nwpx)
              rng = lonarr(2)
              for mm=1L,nwpx-1 do begin
                 rng = [(mm-1)*bin,((mm*bin-1) < (npix-1))]
                 nwwv[mm-1] = total(wave[rng[0]:rng[1]])/bin
                 nwfx[mm-1] = total(fx[rng[0]:rng[1]])/bin
                 nwsig[mm-1] = total(sig[rng[0]:rng[1]])/bin
              endfor 

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

        ;; Set vel array
        x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                     PIXMAX=pmx, VELO=velo

        ;; Y-axis

        if ymnx[1]-ymnx[0] GT 0.6 then ytint = 0.4 else ytint = 0

        ;; Plot
        ysty=1
        ;;if not last line and not last column  and not at end of lines
        if keyword_set(outline) then begin
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           ;; Indicate range used to measure EW, added by KLC

           if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
              ((ii+1) MOD (npx*npy) NE 0) $
              (((ii+1) MOD npy) NE 0) then begin
              spaces = replicate(' ',30)
              ;;xtickn is xtickname which is annotation for x-axis
              ;;NODATA will be activated if nblnd NE 0 xtickinterval
              ;;only activated if xtint set in this program 
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xtickn=spaces, xmargin=xmrg, $
                    ymargin=ymrg, /NODATA, charsize=csize, psym=10, $
                    background=clr.white, color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2],xtickinterval=xtint,$
                    /xstyle, /ystyle, thick=lthick, ytickinterval=ytint
           endif else begin
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xmargin=xmrg, ymargin=ymrg, /NODATA, $
                    charsize=csize, psym=10, background=clr.white, $
                    color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2],$
                    /xstyle, /ystyle, thick=lthick, xtickinterval=xtint, $
                    xtitle='Relative Velocity (km s!E-1!N)', ytickinterval=ytint
           endelse
           oplot,velo[pmn:pmx],fx[pmn:pmx],psym=10,color=clr.black,$
                 thick=lthick

           vmnbnd = 2.998e5*(strct[a[0]].wv_lim[0]/$
                             (strct[a[0]].wrest*(1.+zabs)) - 1.)
           vmxbnd = 2.998e5*(strct[a[0]].wv_lim[1]/$
                             (strct[a[0]].wrest*(1.+zabs)) - 1.)
           dum = min(abs(velo-vmnbnd),vmn)
           dum = min(abs(velo-vmxbnd),vmx)
           oplot,velo[vmn:vmx],fx[vmn:vmx],psym=10,color=clr.black,$
                 lthick=3.*lthick
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        endif else begin
           if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
              ((ii+1) MOD (npx*npy) NE 0) $
              (((ii+1) MOD npy) NE 0) then begin
              spaces = replicate(' ',30)
              ;;xtickn is xtickname which is annotation for x-axis
              ;;NODATA will be activated if nblnd NE 0 xtickinterval
              ;;only activated if xtint set in this program 
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xtickn=spaces, xmargin=xmrg, $
                    ymargin=ymrg, NODATA=nblnd, charsize=csize, psym=10, $
                    background=clr.white, color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2],xtintinterval=xtint,$
                    /xstyle, /ystyle, thick=lthick, ytickinterval=ytint,$
                    linestyle=lsty
           endif else begin
              plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, /xminor,/yminor,$
                    yrange=ymnx, xmargin=xmrg, ymargin=ymrg, NODATA=nblnd, $
                    charsize=csize, psym=10, background=clr.white, $
                    color=clr.black, $
                    yticks=1,ytickv=[0.5,1.2],$
                    /xstyle, /ystyle, thick=lthick, xtickinterval=xtint, $
                    xtitle='Relative Velocity (km s!E-1!N)', ytickinterval=ytint, $
                    linestyle=lsty
           endelse
        endelse 
        if keyword_set(plterr) then $
           oplot,velo[pmn:pmx],sig[pmn:pmx],color=clr.red,thick=lthick

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Overplot Voigt profile with given N and b
        ;; Must preserve values for primary spectrum
        ;; added by KLC 16 Oct 2006
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if keyword_set(voigt) then begin
           vplin = x_setline(wrest)
           vplin.zabs = strct[a].zabs
           if keyword_set(vzabs) then vplin.zabs = vzabs else $
              ;; 30 km/s wiggle room
              if abs(strct[a].zabs-zabs) gt 0.0001 then vplin.zabs = zabs 

           sz = size(voigt,/dimensions)
           if size(voigt,/n_dimensions) eq 1 then begin ;one element to overlay VP
              case sz[0] of
                 2: begin       ;default: only HI lines
                    if stregex(nam,'HI',/fold_case,/boolean) then begin
                       vplin.n = float(voigt[0])
                       vplin.b = float(voigt[1])
                    endif else vplin.n = -1
                 end
                 3: begin       ;only specified element
                    if stregex(nam,voigt[0],/fold_case,/boolean) then begin
                       vplin.n = float(voigt[1])
                       vplin.b = float(voigt[2])
                    endif else vplin.n = -1
                 end
                 else: stop,'fuse_velplt: voigt wrong size'
              endcase 
           endif else begin     ;multiple elements to overlay VP
              for vv=0,sz[1]-1 do begin
                 if stregex(nam,voigt[0,vv],/fold_case,/boolean) then begin
                    vplin.n = float(voigt[1,vv])
                    vplin.b = float(voigt[2,vv])
                    break
                 endif else vplin.n = -1
              endfor 
           endelse 

           if instr le 6 then fwhm = 5 $ ;FUSE
           else fwhm = 2                 ;STIS
           if vplin.n ne -1 then $       ;check current line correct
              vptau = x_voigt(wave[pmn:pmx], vplin, fwhm=fwhm) $
           else vptau = replicate(1.,pmx-pmn+1)

                                ;if n_elements(voigt) ne 2 then $
                                ;  stop,'fuse_velplt: voigt wrong size'
                                ;
                                ;vplin = x_setline(wrest)
                                ;vplin.zabs = zabs
                                ;vplin.n = voigt[0]
                                ;vplin.b = voigt[1]
                                ;if instr le 6 then fwhm = 5 $ ;FUSE
                                ;else fwhm = 2       ;STIS
                                ;vptau = x_voigt(wave[pmn:pmx], vplin, fwhm=fwhm)
           vpvelo = velo[pmn:pmx]
        endif
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        ;; Overplot other detections
        if na gt 0 then begin
           for jj=1,na-1 do begin
              if instr[jj] NE svinst then begin
                 svinst = instr[jj]
                 ;;  Note that STIS starts at instrument 8 (7 wrt index)
                 if instr[jj] LE 6 then $
                    fx = x_readspec(instr_fil[instr[jj]], SIG=sig, $
                                    wav=wave,NPIX=npix, inflg=3)$
                 else begin     ; STIS
                    spos = strpos(instr_fil[instr[jj]], suffix[0])
                    sig_fil = strmid(instr_fil[instr[jj]], 0, spos)+$
                              suffix[1]
                    fx = x_readspec(instr_fil[instr[jj]], SIG=sig, $
                                    wav=wave, $
                                    NPIX=npix,fil_sig=sig_fil, inflg=0)
                 endelse
                 ;; Shift wave
                 wave = wave + inst_dw[instr[jj]]*inst_w0[instr[jj]]/wave

                 if keyword_set(bin) then begin
                    nwpx = npix/bin + 1
                    nwwv = fltarr(nwpx)
                    nwfx = fltarr(nwpx)
                    nwsig = fltarr(nwpx)
                    rng = lonarr(2)
                    for mm=1L,nwpx-1 do begin
                       rng = [(mm-1)*bin,((mm*bin-1) < (npix-1))]
                       nwwv[mm-1] = total(wave[rng[0]:rng[1]])/bin
                       nwfx[mm-1] = total(fx[rng[0]:rng[1]])/bin
                       nwsig[mm-1] = total(sig[rng[0]:rng[1]])/bin
                    endfor 

                    wave=nwwv
                    fx = nwfx
                    sig = nwsig
                 endif 

                 ;; Sort?
                 srt = sort(wave)
                 wave = wave[srt]
                 fx = fx[srt]
                 sig = sig[srt]
              endif else begin
                 if keyword_set(MSGFIL) then begin
                    printf,17,'# Multiple lines within tolerance: '+$
                           'wrest, wobs, zabs?, wrest?, wobs?, instr?'
                    printf,17,wrest,wrest*(1.+zabs),strct[a[jj]].zabs,$
                           strct[a[jj]].wrest,$
                           strct[a[jj]].wrest*(1.+strct[a[jj]].zabs),$
                           strct[a[jj]].instr,$
                           format='(2(f9.4,tr2),tr1,f8.6,tr2,2(f9.4,tr2),i3)'
                 endif else begin
                    print,'fuse_velplt: Confused detections! ',zabs
                    print,'Skipping ',strct[a[jj]].wrest,' at ',$
                          strct[a[jj]].zabs
                 endelse 
                 goto,next_line
                                ;stop
              endelse 

              ;; Set vel array
              x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                           PIXMAX=pmx, VELO=velo
              
              if keyword_set(outline) then begin
                 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                 ;; Indicate range used to measure EW, added by KLC
                 vmnbnd = 2.998e5*(strct[a].wv_lim[0]/$
                                   (strct[a].wrest*(1.+zabs)) - 1.)
                 vmxbnd = 2.998e5*(strct[a].wv_lim[1]/$
                                   (strct[a].wrest*(1.+zabs)) - 1.)
                 dum = min(abs(velo-vmnbnd),vmn)
                 dum = min(abs(velo-vmxbnd),vmx)
                 oplot,velo[vmn:vmx],fx[vmn:vmx],psym=10,thick=3*lthick,$
                       color=clrspec[instr[jj]]
                 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              endif else $
                 oplot,velo[pmn:pmx],fx[pmn:pmx],color=clrspec[instr[jj]],$
                       psym=10,thick=lthick,linestyle=lsty

           endfor 
           next_line:
                                ;print,'Finished plot ',ii
           !p.color = clr.black
        endif


        goto,skip_blnds
        ;; BLENDS
        nlow = pmn
        for jj=0L,nblnd-1 do begin
           readf, 11, vmin, vmax, FORMAT='(f,f)'
           mn = min(abs(velo-vmin),pixmin)
           mn = min(abs(velo-vmax),pixmax)

           ;; Plot good
           nlow = nlow < pixmin
           oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
                  color=clr.black, psym=10, thick=lthick
           ;; Plot blend
           oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
                  psym=10, linestyle=2, thick=1.
           ;; End
           if jj EQ nblnd-1 then begin
              if pixmax LT pmx then $
                 oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                        psym=10, thick=lthick
           endif
           ;; Set nlow
           nlow = pixmax
        endfor
        skip_blnds:


        ;; Labels
        getfnam, wrest, fv, nam, _extra=extra ;eg fil=getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
        ;;Label H2 features
        h2 = where(strcmp(strtrim(nam,2),''),nh2)
        if nh2 ne 0 then $
           nam[h2] = 'H!D2!N '+string(wrest[h2],format='(f6.1)')

        xpos = 0.02*(vmnx[1]-vmnx[0])+vmnx[0]
        ypos = ymnx[0]+ (ymnx[1]-ymnx[0])*0.10
        xyouts, xpos, ypos,strtrim(nam,2), color=clr.black, charsize=LSIZE
        
        wobs = wrest*(1.+zabs)
        wobsnam = ' ('+strmid(strtrim(wobs,2),0,7) + $ ;xxxx.xx
                  ', '+strtrim(2^instr[0],2)
        ;; Cross correlate each line with all other systems and
        ;; indicate overlap
        cross = where(abs(strct.wrest*(1.+strct.zabs)-wobs) LT 0.05 AND $
                      strct.wrest NE wrest,ncross)
        if ncross NE 0 then begin
           wobsnam += '; '+strtrim(ncross,2)+')' 
           if keyword_set(MSGFIL) then begin
              printf,17,'# Possible identifications for ',strtrim(nam,2),$
                     ' at ',strtrim(zabs,2),' (',strtrim(wobs,2),'):'
              for kk=0,ncross-1 do $
                 printf,17, strct[cross[kk]].zabs,strct[cross[kk]].wrest,$
                        (1.+strct[cross[kk]].zabs)*strct[cross[kk]].wrest,$
                        strct[cross[kk]].instr,format='(f8.6,2(tr2,f9.4),tr2,i3)'
           endif else begin
              print,'Line overlap for ',zabs,wobs,':'
              for kk=0,ncross-1 do $
                 print,strct[cross[kk]].zabs,strct[cross[kk]].wrest
           endelse 
        endif else wobsnam += ')'                   
        xyouts, wobsnam, color=clr.black,charsize=0.5*LSIZE
        ;; Reshift label upper right corner
        xyouts,0.9,0.98,'z!Dabs!N = '+strtrim(zabs,2),$
               charsize=0.5*LSIZE,alignment=0.5,/normal
        ;; Print input file name in upper left corner
        dum = [strpos(vel_fil,'/',/reverse_search)+1,$
               strpos(vel_fil,'.',/reverse_search)]
        xyouts,0.25,0.98,strmid(vel_fil,dum[0],dum[1]-dum[0]),$
               charsize=0.5*LSIZE,alignment=0.5,/normal

        if (ii+1) MOD ntot EQ 1 then begin
;           if npx EQ 1 then $
;              xyouts, 0.045, 0.5, 'Normalized Flux', $
;                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize $
;           else $
              xyouts, 0.045, 0.5, 'Normalized Flux', $
                      alignment=0.5, ORIENTATION=90., /normal, charsize=lsize 
        endif
        
        ;; Lines
        if keyword_set(plterr) then begin
           oplot, [0., 0.], ymnx, color=clr.black, linestyle=2, $
                  thick=2.5*lthick
           oplot, [-10000., 10000.], [0.,0.], color=clr.blue, linestyle=3, $
                  thick=2.5*lthick
        endif else begin
           oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, $
                  thick=2.5*lthick
           oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3, $
                  thick=2.5*lthick
        endelse

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Overplot Voigt profile with given N and b
        ;; added by KLC 16 Oct 2006
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if keyword_set(voigt) then oplot,vpvelo,vptau,color=clr.limegreen $
        else oplot, [-10000., 10000.], [1.,1.], color=clr.limegreen, linestyle=3
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     endfor
  endelse

  close,11
  if keyword_set( PSFILE ) then x_psclose
  if keyword_set( LAPTOP ) then x_tiffset, /unset
  if keyword_set( MSGFIL ) then begin
     close,17
     print,'fuse_velplt: ',msgfil,' created'
  endif
  !p.multi=[0,1,1]
  !y.omargin=[0,0]

  return
end
