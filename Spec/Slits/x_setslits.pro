;+ 
; NAME:
;   x_setslits
; PURPOSE:
;    Finds slit edges on the 'flattened' flat in the presence of rotation
; CALLING SEQUENCE:
;   slackers_setslits, flatimg, slitstr
; INPUTS:
;   flat        - Flat image or fits file
;   slitstr     - Slit structure
; COMMENTS:
;   This piece of code allows for a slight rotation (-2 to 2 deg) of
;   the mask, which can cause overlapping slits and other annoying
;   things.
; REVISION HISTORY:
;  14-Feb-2002 Written by JXPj
;   6-Apr-2004 Revised by MRB 
;-
;------------------------------------------------------------------------------
pro match_slits, peak, guess, imatch

ipeak=sort(peak)
iguess=sort(guess)

j=0L
imatch=lonarr(n_elements(guess))-1L
for i=0L, n_elements(iguess)-1L do begin
    while(peak[ipeak[j]] lt guess[iguess[i]] and $
          j lt n_elements(ipeak)-1L) do j=j+1L
    jlo=j-1L
    jhi=j
    if(j eq 0L) then begin
        jlo=0L
        jhi=1L
    endif
    if(j eq n_elements(peak)-1L) then begin
        jlo=n_elements(peak)-2L
        jhi=n_elements(peak)-1L
    endif
    dm0=abs(peak[ipeak[jlo]]-guess[iguess[i]])
    dm1=abs(peak[ipeak[jhi]]-guess[iguess[i]])
    if(dm0 le dm1) then imatch[iguess[i]]=jlo else imatch[iguess[i]]=jhi
endfor

end
;
pro x_setslits, flat, slitstr, YSTRT=ystrt, DEBUG=debug, $
                       MINMAXTHETA=minmaxtheta, NTHETA=ntheta, $
                       MINMAXSHIFT=minmaxshift, NSHIFT=nshift, $
                       NOFIT=nofit, THETA=theta, SHIFT=shift

; Error catching
if  N_params() LT 2 then begin 
    print,'Syntax - ' + $
      'slackers_setslits, flat, slitstr'
    return
endif 

; Optional Keywords
if not keyword_set( YSTRT ) then ystrt = 400L
if not keyword_set( MINMAXTHETA ) then minmaxtheta = [-2.,2.] ; in deg
if not keyword_set( NTHETA ) then ntheta=200L
if not keyword_set( MINMAXSHIFT ) then minmaxshift = [-50.,50.] ; in pix
if not keyword_set( NSHIFT ) then nshift=200L
if not keyword_set( WIDTH ) then width=80L
if not keyword_set( NSIG ) then nsig=50.
if not keyword_set( CLOSESLIT ) then closeslit=5.

if(keyword_set(NOFIT)) then begin
    nx=(size(flat,/dim))[0]
    ny=(size(flat,/dim))[1]
    nslit=n_elements(slitstr)
    slitstr.yedg_flt=slitstr.yedg-shift-theta*!DPI* $
      (replicate(1.,2)#slitstr.xpos-0.5*float(nx))/180.
endif else begin

; Allow flat to be fits file
    dat = x_readimg(flat, /fscale)
    sz_img = size(dat, /dimensions)

; Now find the initial guess for slit positions by trying a whole
; bunch of offsets and rotation angles. 
    nslit = n_elements(slitstr)

    mintheta=minmaxtheta[0]*!DPI/180.
    maxtheta=minmaxtheta[1]*!DPI/180.
    thetas=mintheta+(maxtheta-mintheta)*(findgen(ntheta)+0.5)/float(ntheta)
    minshift=minmaxshift[0]
    maxshift=minmaxshift[1]
    shifts=minshift+(maxshift-minshift)*(findgen(nshift)+0.5)/float(nshift)

    smsh = djs_median(dat[ystrt-width/2L:ystrt+width/2L,*],1)
    saw = shift(smsh,1) - shift(smsh,-1)
    x_fndpeaks, abs(saw), center, NSIG=nsig, PEAK=peak, /thin, $
      /force, pkwdth=1, /silent

    nx=(size(flat,/dim))[0]
    ny=(size(flat,/dim))[1]
    chi2=fltarr(ntheta,nshift)
    for i=0L, ntheta-1L do begin
        for j=0L, nshift-1L do begin
            guess=slitstr.yedg-shifts[j]-thetas[i]* $
              (replicate(1.,2)#slitstr.xpos-0.5*float(nx))
            match_slits, peak, guess, imatch
            chi2[i,j]=total((peak[imatch]-guess)^2,/double)
        endfor
    endfor
    minchi2=min(chi2, imin)
    theta=thetas[imin mod ntheta]
    shift=shifts[imin/ntheta]
    guess=slitstr.yedg-shift-theta* $
      (replicate(1.,2)#slitstr.xpos-0.5*float(nx))

; Proceed with this guess --- find the nearest bottom to each bottom,
;                             the nearest top to each top, within the
;                             distance CLOSESLIT; if none, revert to
;                             the guess
    tops=peak[where(saw[peak] GT 0, ntop)]
    if(ntop eq 0) then message, 'no tops in the flat, must not be a flat image'
    match_slits, tops, guess[1,*], imatch
    for i=0L, nslit-1L do begin
        dm=abs(guess[1,i]-tops[imatch[i]])
        if(dm lt closeslit) then begin
            slitstr[i].yedg_flt[1]=tops[imatch[i]]
        endif else begin
            slitstr[i].yedg_flt[1]=guess[1,i]
            splog,'using guess for top of slit '+string(i)+' at '+ $
              string(guess[1,i])
        endelse
    endfor

    bots=peak[where(saw[peak] LT 0, nbot)]
    if(nbot eq 0) then message, 'no bots in the flat, must not be a flat image'
    match_slits, bots, guess[0,*], imatch
    for i=0L, nslit-1L do begin
        dm=abs(guess[0,i]-bots[imatch[i]])
        if(dm lt closeslit) then begin
            slitstr[i].yedg_flt[0]=bots[imatch[i]]
        endif else begin
            slitstr[i].yedg_flt[0]=guess[0,i]
            splog,'using guess for bot of slit '+string(i)+' at '+ $
              string(guess[0,i])
        endelse
    endfor
    slitstr.yedg_flt=(slitstr.yedg_flt > 1.) < (float(ny)-2.)

; check for zero-width slits
    for i=0L, nslit-1L do begin
        if(slitstr[i].yedg_flt[0] ge slitstr[i].yedg_flt[1]) then $
          message, 'zero-width slit found, that is not supposed to happen'
    endfor

endelse

; Now that we have the full slits determined, make them
; non-overlapping 
;   a- count number of slits covering each row
ninrow=lonarr(ny)
for i=0L, nslit-1L do begin
    rlo=long(djs_floor(slitstr[i].yedg_flt[0]))
    rhi=long(djs_ceil(slitstr[i].yedg_flt[1]))
    ninrow[rlo:rhi]=ninrow[rlo:rhi]+1L
endfor
;   b- within each slit, find largest contiguous set with exactly one coverage 
for i=0L, nslit-1L do begin
    rlo=long(djs_floor(slitstr[i].yedg_flt[0]))
    rhi=long(djs_ceil(slitstr[i].yedg_flt[1]))
    ione=where(ninrow[rlo:rhi] eq 1, none)
    if(none eq 0) then begin 
        splog, 'lost entire slit '+string(i)+ $
          ' (priority= '+string(slitstr[i].priority)+'), setting flg_anly=0'
        slitstr[i].flg_anly=0
    endif else begin 
        inot=where(ninrow[rlo:rhi] ne 1, nnot)
        if(nnot gt 0) then begin
            ncontig=0L
            for r=rlo, rhi do begin
                for rp=r, rhi do begin
                    if(rp-r+1L gt ncontig) then begin
                        inot=where(ninrow[r:rp] ne 1, nnot)
                        if(nnot eq 0) then begin
                            contighi=rp
                            contiglo=r
                            ncontig=rp-r+1L
                        endif
                    endif
                endfor
            endfor
            slitstr[i].yedg_flt[0]=contiglo-0.5
            slitstr[i].yedg_flt[1]=contighi+0.5
        endif
    endelse
endfor

if(NOT keyword_set(NOFIT)) then begin
; Finally, determine which are vignetted and set them to flg_anly==0
    inone=where(ninrow eq 0, none)
    background=median(smsh[inone])
    sigback=djsig(smsh[inone],sigrej=3)
    slitflux=fltarr(nslit)
    for i=0L, nslit-1L do begin
        rlo=long(djs_floor(slitstr[i].yedg_flt[0]))
        rhi=long(djs_ceil(slitstr[i].yedg_flt[1]))
        slitflux[i]=median(smsh[rlo:rhi])
        if(slitflux[i] lt background+sigback) then begin
            splog, 'looks like vignetted slit '+string(i)+ $
              ' (priority= '+string(slitstr[i].priority)+ $
              '), setting flg_anly=0'
            slitstr[i].flg_anly=0
        endif
    endfor
endif
    
return
end

