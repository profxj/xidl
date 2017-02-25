;+
; NAME:
; g_mkvoigt
; Version 1.0
;
; PURPOSE:
;   create voigt profile flux arrays for output files from vpfit
;
; CALLING SEQUENCE:
;   g_mkvoigt, vpfil, fluxfil, TLIST=, FITSFIL=, PSFIL=, /DAT
;
; INPUTS:
;  vpfil       - output fort.26 vpfil
;  fluxfil     - wavelength, flux, sig of qso to be fit
;  DAT         - set this flag if the input fluxfil is ASCII 
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  FITSFIL     - output fits file
;
; OPTIONAL OUTPUTS:
;  PSFIL       - output postscript plots of fitting regions
;
; COMMENTS:
;   this procedure (because of call to x_voigt) is VERY expensive!
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   x_voigt
;
; REVISION HISTORY:
;   21-Feb-2005  Created by GEP
;  
;-
;------------------------------------------------------------------------------


pro g_mkvoigt, vpfil, fluxfil, TLIST=tlist, DEBUG=debug, FITSFIL=fitsfil, $
               PSFIL=psfil, DAT=dat, VEL=vel, ZABS=ZABS

    if (N_params() LT 2) then begin
        print, 'Syntax - ' + $
               'g_mkvoigt, vpfil, fluxfil, TLIST=tlist, DEBUG=debug, FITSFIL=fitsfil,' + $
               ' PSFIL=psfil, DAT=dat, VEL=vel, ZABS=ZABS'
        return
    endif


    all_lin = xmrdfits(getenv('XIDL_DIR')+'/Spec/Lines/all_lin.fits',$
                        1,structyp='abslinstrct',/silent)

    g_vpparse, vpfil, VPSTR=vpstr

    bd = where(strmid(vpstr.ion,1,1) EQ ' ',nbd)
    if nbd NE 0 then begin              ;  Need to fix the names to match 'all_lin'
        vpstr.ion[bd] = strmid(vpstr.ion[bd],0,1)+strmid(vpstr.ion[bd],2)
    endif
    name = vpstr.ion[0:vpstr.nion-1]

    lls = where(name EQ 'HI',nlls)
    metals = where(name NE 'HI',nmetals)
    nmet = 0
    for i=0,nmetals-1 do begin
        use = where(strmatch(strtrim(all_lin.ion,2),'*'+name[metals[i]]+' *') EQ 1,nuse)
        nmet = nmet+nuse
    endfor

    nlines = nlls*19+nmet
    nset = n_elements(name)
    
    tmp = {abslinstrct}
    lines = replicate(tmp,nlines)
    
    ;; SET LINES

    if not keyword_set(tlist) then tlist = getenv('XIDL_DIR')+'/LLS/Lines/LLS_std.lst'
    
    readcol, tlist, wr, skipline=1

    nlin=0
    for i=0,nset-1 do begin
        if name[i] EQ 'HI' then begin
            for j=0,n_elements(wr)-1 do begin
                tmp = x_setline(wr[j])
                lines[nlin+j] = tmp
                lines[nlin+j].n = vpstr.n[i]
                lines[nlin+j].nsig = vpstr.nerr[i]
                lines[nlin+j].b = vpstr.b[i]
                lines[nlin+j].bsig = vpstr.berr[i]
                lines[nlin+j].zabs = vpstr.z[i]
                lines[nlin+j].zsig = vpstr.zerr[i]
                lines[nlin+j].set = i+1
            endfor
            nlin = nlin+n_elements(wr)
        endif else begin
            use = where(strmatch(strtrim(all_lin.ion,2),'*'+name[i]+' *') EQ 1,nuse)
            for j=0, nuse-1 do begin
                lines[nlin+j] = all_lin[use[j]]
                lines[nlin+j].n = vpstr.n[i]
                lines[nlin+j].nsig = vpstr.nerr[i]
                lines[nlin+j].b = vpstr.b[i]
                lines[nlin+j].bsig = vpstr.berr[i]
                lines[nlin+j].zabs = vpstr.z[i]
                lines[nlin+j].zsig = vpstr.zerr[i]
                lines[nlin+j].set = i+1
            endfor
            nlin = nlin+nuse
        endelse
    endfor

    if keyword_set(DAT) then $
        readcol, fluxfil, wave, flux, sig, FORMAT='F,F,F' $
    else begin
        wave = xmrdfits(fluxfil)
        flux = xmrdfits(fluxfil,1)
        sig = xmrdfits(fluxfil,2)
    endelse
    npix = n_elements(wave)
    fit = replicate(1.,npix)
    ifit = replicate(1.,nset,npix)

    nwave = vpstr.nreg

    if not keyword_set(VEL) then begin
        for i=0,nset-1 do begin
            for j=0,nwave-1 do begin
                use = where(lines.set EQ i+1, nuse)
                gd = where(wave GE vpstr.reg_beg[j] AND wave LE vpstr.reg_end[j],ngd)
                if max(wave[gd]) LT 4780. then FWHM = 3.384 else FWHM = 2.881
                ifit[i,gd] = x_voigt(wave[gd],lines[use],FWHM=FWHM)
            endfor
        endfor
        for j=0,nwave-1 do begin
            use = where(lines.set NE 0)
            gd = where(wave GE vpstr.reg_beg[j] AND wave LE vpstr.reg_end[j],ngd)
            if max(wave[gd]) LT 4780. then FWHM = 3.384 else FWHM = 2.881
            fit[gd] = x_voigt(wave[gd],lines[use],FWHM=FWHM)
        endfor
    endif
    if keyword_set(DEBUG) then stop

    if keyword_set(PSFIL) then begin
        clr = getcolor(/load)
        x_psopen, psfil, /portrait
        !p.multi=[0,2,6,0,1]
    endif
    if not keyword_set(VEL) then begin
        for j=0,nwave-1 do begin
            spaces = replicate('!17 ',30)
            gd = where(wave GE vpstr.reg_beg[j] AND wave LE vpstr.reg_end[j],ngd)
            plot, wave[gd], flux[gd], xrange=[vpstr.reg_beg[j],vpstr.reg_end[j]], $
                  yrange=[-0.11,1.09], $
                  xtickn=spaces, xmargin=[9,3], ymargin=[0,0], charsize=1.8, psym=10, $
                  background=clr.white, color=clr.black, xstyle=1, ystyle=1, thick=2.5
            oplot, wave[gd], fit[gd], psym=10, color=clr.green, thick=2.5
;            for i=0,nset-1 do begin
;                use = where(lines.set EQ i+1, nuse)
;                oplot, wave[gd], ifit[i,gd], psym=10, color=clr.red
;            endfor
        endfor
    endif 
    if keyword_set(VEL) then begin
        if not keyword_set(ZABS) then begin
            print, 'ZABS not set -- do so now then continue...'
            stop
        endif
        vlin = numlines(vel)
        rwave = fltarr(vlin)
        tname = strarr(vlin)
        vmnx = fltarr(2,vlin)
        numwidths=0
        sysnum = intarr(vlin)
        ovmnx = [0.,0.]

        openr,1,vel

        dumc = ''
        for i=0,vlin-1 do begin
            readf, 1, dumc
            prs = strsplit(dumc,' ',/extract)
            rwave[i] = float(prs[0])
            tname[i] = prs[1]+' '+prs[2]
            vmnx[0,i] = float(prs[3])
            vmnx[1,i] = float(prs[4])
            if ovmnx[0] NE vmnx[0,i] OR ovmnx[1] NE vmnx[1,i] then numwidths = numwidths+1
            ovmnx[0] = vmnx[0,i]
            ovmnx[1] = vmnx[1,i]
            sysnum[i] = numwidths
        endfor
        close, 1
        
        unum = ceil(float(vlin)/2.)
        !p.multi=[0,2,unum,0,1]
        for i=1,numwidths do begin
            use = where(sysnum EQ i,nuse)
            if nuse EQ 0 then stop
            uvmnx = [vmnx[0,use[0]],vmnx[1,use[0]]]
            all_velo = x_allvelo(wave, zabs, rwave[use], vmnx, ALL_PMNX=all_pmnx)
            for j=0, nuse-1 do begin
                luse = where(lines.set NE 0)
                gd = where(wave GE wave[all_pmnx[0,j]] AND wave LE wave[all_pmnx[1,j]],ngd)
                if wave[ngd-1] LE 4780.0 then FWHM=3.384 else FWHM=2.881
                fit[gd] = x_voigt(wave[gd],lines[luse],FWHM=FWHM)
                plot, all_velo[0:all_pmnx[2,j],j], flux[all_pmnx[0,j]:all_pmnx[1,j]],$
                      xrange=[all_velo[0],all_velo[all_pmnx[2,0]]], $
                      yrange=[-0.11,1.09], $
                      xtickn=spaces, xmargin=[9,3], ymargin=[0,0], charsize=1.8, psym=10, $
                      background=clr.white, color=clr.black, xstyle=1, ystyle=1, thick=2.5
                oplot, all_velo[0:all_pmnx[2,j],j], fit[all_pmnx[0,j]:all_pmnx[1,j]],$
                      psym=10, color=clr.green, thick=2.5
                xyouts, 0.07*(vmnx[1,use[j]]-vmnx[0,use[j]])+vmnx[0,use[j]], $
                        1.2*0.05, tname[use[j]], color=clr.black, charsize=1.5
                oplot, [0.,0.], [-0.11,1.09], color=clr.orange, linestyle=3
                oplot, [-10000.,10000.], [0.,0.], color=clr.blue, linestyle=3
                oplot, [-10000.,10000.], [1.,1.], color=clr.blue, linestyle=3
            endfor
        endfor
    endif
    if keyword_set(PSFIL) then x_psclose

    if keyword_set(FITSFIL) then begin
        mwrfits, wave, fitsfil, /create
        mwrfits, flux, fitsfil
        mwrfits, sig, fitsfil
        mwrfits, fit, fitsfil
        for i=0,nset-1 do begin
          mwrfits, ifit[i], fitsfil
        endfor
    endif

end
