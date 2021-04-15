;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_velplt.pro
; Author: Kathy Cooksey                      Date: 9 Apr 2013
; Project: Precious Metals in SDSS Quasar Spectra
; Description: Plot a velocity plot given spectrum and redshift
; Input:
;   spec_fil -- file name
;   zabs -- redshift of system
; Optional:
;   dir= -- directory path to prepend to spec_fil (default: ./)
;   linlst= -- standard XIDL line list
;              (default: $XIDL_DIR/SDSS/CIV/sdss_civ.lst)
;   nx= -- number of columns (default: 2)
;   ny= -- number of rows (default: 6)
;   csize= -- character/font size (default: 3)
;   lthick= -- thickness of lines (default: 5)
;   xrng= -- 2-element array of velocity limits for plot (km/s)
;            (default: [-950,950])
;   yrng= -- 2-element array of (normalized) flux for plot
;            (default: [-0.1,1,3])
;   title= -- to print to top of postscript
;   label= -- flag (number) or input (string) to show in panels
;   nbin= -- number of pixels to smooth over for display purposes
;   _extra includes:
;      inflg= for x_readspec() [5 for SDSS]
;      fil_sig= for x_readspec()
; Output:
;   psfil= -- name of postscript to plot to (default: sdss_velplt.ps).
;             If number of lines greater than nx*ny, then use root 
;             to create multiple outputs (e.g., sdss_velplta.ps)
; Example:
; History:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_velplt,spec_fil,zabs,norm=norm,dir=dir,linlst=linlst,$
                nx=nx,ny=ny, csize=csize,lthick=lthick,$
                xrng=xrng,yrng=yrng,psfil=psfil,title=title, $
                label=label,nbin=nbin,spec2=spec2,_extra=extra
  
  c = 2.998e5                   ; km/s

  ;; Params
  angstrom = STRING("305B)   
  gesign = ' !9'+string("263B)+'!X '
  lesign = ' !9'+string("243B)+'!X '
  times = '!9'+string("264B)+'!X'
  divid = '!9'+string("270B)+'!X'
  langle = '!9'+string("341B)+'!X'
  rangle = '!9'+string("361B)+'!X'
  sqrt = "!9"+string("326B)+"!X"

  ;; Linestyle
  ;; 0: solid, 1; dotted, 2: dashed, 3: dash-dot, 4: dash-dot-dot-dot,
  ;; 5: long dash
  ;; Psym
  ;; 1: +, 2: *, 3: ., 4: Diamond, 5: Triangle, 6: Square, 7: X, 8: user, 9: undef, 10: histogram

  if not keyword_set(dir) then dir = './'
  if not keyword_set(nx) then nx = 2 ; columns
  if not keyword_set(ny) then ny = 6 ; rows
  if not keyword_set(linlst) then $
     linlst = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civ.lst'
  if size(linlst,/type) eq 7 then begin
     if size(linlst,/n_dim) gt 1 then begin
        wrest = float(linlst[*,0])
        ionnam = linlst[*,1]
     endif else begin
        readcol,linlst,wrest,ion,skip=1,format='(f9.4,a5)',count=nlin
        ionnam = strtrim(ion,2)+' '+strtrim(floor(wrest),2)
     endelse
  endif else begin
     wrest = linlst.wrest
     ionnam = linlst.ion
  endelse 
  nion = n_elements(wrest)
  if not keyword_set(psfil) then psfil = 'sdss_velplt.ps'

  ;; Test and call recursively
  if nion gt nx*ny then begin
     ;; calculate plots
     nplot = ceil(nion/float(nx*ny))
     psroot = strmid(psfil,0,strpos(psfil,'.',/reverse_search))
     istart = 0
     istop = nx*ny - 1
     for pp=0,nplot-1 do begin
        ltr = string(97B+byte(pp)) ; letter
        sub_linlst = strarr(istop-istart+1,2)
        sub_linlst[*,0] = wrest[istart:istop]
        sub_linlst[*,1] = ionnam[istart:istop]
        sdss_velplt,spec_fil,zabs,norm=norm,dir=dir,linlst=sub_linlst,nx=nx,ny=ny,$
                    csize=csize,lthick=lthick,$
                    xrng=xrng,yrng=yrng,psfil=psroot+ltr+'.ps',title=title, $
                    label=label,_extra=extra
        istart = istop + 1
        istop = (istart + nx*ny - 1) < (nion-1)
     endfor ; loop pp=nplot
     return ; exit
  endif

  if not keyword_set(csize) then csize = 3
  if not keyword_set(lthick) then lthick = 5.

  if not keyword_set(yrng) then yrng = [-0.1,1.3]   ; flux
  if not keyword_set(xrng) then xrng = [-950.,950.] ; km/s
  dx = xrng[1] - xrng[0]
  dy = yrng[1] - yrng[0]
  if not keyword_set(title) then title = ['Relative Velocity (km s!E-1!N)','Normalized Flux']

  ;; Plot
  x_psopen,psfil,/portrait
  ;; [0] plots left on page; [1] no. col. per page; [2] no. row per
  ;; page; [3] no. plots in z dir.; [4] 0 for L->R, T->B; 1 for T->B,
  ;; L->R 
  !p.multi = [0,nx,ny,0,1]
  !x.omargin = [4,1]         ; left and right border
  !y.omargin = [4,0.5]         ; bottom and top border

  clr = getcolor(/load)

  ;; read spectrum
  flux = x_readspec(spec_fil,wav=wave,sig=sigma,_extra=extra) ; inflg=, fil_sig=
  if keyword_set(norm) then begin
     ;; Normalize but need wavelength from above
     sdss_normspec, spec_fil, norm, norm_spec, /ret_spec, _extra=extra
     flux = norm_spec[*,0]
     sigma = norm_spec[*,2]
  endif 

  if keyword_set(spec2) then begin
;     tags = tag_names(extra)
;     chk = where(stregex(tags,'INFLG',/boolean))
;     if chk[
     fx2 = x_readspec(spec2[0],wav=wv2,sig=er2,fil_sig=spec2[1])
     flux = [flux,fx2]
     wave = [wave,wv2]
     sigma = [sigma,er2]
  endif
;  parse_sdss,dir+spec_fil,flux,wave,sig=sigma


  ;; Smooth
  if keyword_set(nbin) then begin
     wave = smooth(wave,nbin,/edge_truncate)
     flux = smooth(flux,nbin,/edge_truncate)
     sigma = smooth(sigma,nbin,/edge_truncate)
  endif
     
  ;; Mask
  bd = where(finite(flux) eq 0)
  if bd[0] ne -1 then begin
     flux[bd] = !values.f_nan
     sigma[bd] = !values.f_nan
  endif 


  spaces = replicate(' ',30) 
  ymrg = [0,0]
  xmrg = [4,0]
  ymnr = 1
  xmnr = 1
  iion = 0
  for xx=0,nx-1 do begin
     for yy=0,ny-1 do begin
        
        if (xx eq 0) or (yy eq ny-1) or (iion eq nion-1) then begin
           if (xx eq 0) and (yy eq ny-1) then $
              ;; All axis
              plot,xrng,yrng,/nodata,ystyle=1,xstyle=1,background=clr.white,$
              color=clr.black,charsize=csize,xminor=xmnr,yminor=ymnr,$
              xmargin=xmrg,ymargin=ymrg,$
              yticks=3,ytickv=[0.,0.5,1.] $
           else begin
              if xx eq 0 then $
                 ;; Vertical only
                 plot,xrng,yrng,/nodata,ystyle=1,xstyle=1,background=clr.white,$
                 color=clr.black,charsize=csize,xminor=xmnr,yminor=ymnr,xmargin=xmrg,ymargin=ymrg,$
                 yticks=3,ytickv=[0.,0.5,1.],xtickn=spaces $
              else $
                 ;; Horizontal only
                 plot,xrng,yrng,/nodata,ystyle=1,xstyle=1,background=clr.white,$
                 color=clr.black,charsize=csize,xminor=xmnr,yminor=ymnr,xmargin=xmrg,ymargin=ymrg,$
                 yticks=3,ytickv=[0.,0.5,1.],ytickn=spaces
           endelse
        endif else begin
           ;; No axis 
           plot,xrng,yrng,/nodata,ystyle=1,xstyle=1,background=clr.white,$
                color=clr.black,charsize=csize,xminor=xmnr,yminor=ymnr,xmargin=xmrg,ymargin=ymrg,$
                yticks=3,ytickv=[0.,0.5,1.],xtickn=spaces,ytickn=spaces
        endelse

        ;; Set scale
        wobs = wrest[iion]*(1+zabs)
        velo = (wave/wobs - 1) * c
        
        ;; Verify in range
        gd = where(velo ge xrng[0] and velo lt xrng[1],ngd)
        while gd[0] eq -1 do begin
           iion++
           if iion ge nion then begin
              gd = !values.f_nan
              continue
           endif
              
           wobs = wrest[iion]*(1+zabs)
           velo = (wave/wobs - 1) * c
           gd = where(velo ge xrng[0] and velo lt xrng[1],ngd)
        endwhile

        if iion ge nion then break
        
        oplot,velo,sigma,psym=10,color=clr.red,thick=lthick
        oplot,velo,flux,psym=10,color=clr.black,thick=lthick

        ;; Label
        xyouts,0.5,0.85*dy+yrng[0],ionnam[iion],charsize=0.5*csize,$
               color=clr.black,alignment=0.5

        ;; Extra lines
        oplot,xrng,[1,1],linestyle=3,thick=lthick,color=clr.limegreen ; conti
        oplot,xrng,[0,0],linestyle=3,thick=lthick,color=clr.blue ; 0
        oplot,[0,0],yrng,linestyle=1,thick=lthick,color=clr.slategray; vertical

        iion++
        if iion ge nion then break
     endfor                     ; loop yy=ny
     if iion ge nion then break
  endfor                        ; loop xx=nx
  
  ;; Add labels
  if nion le nx then $          ; one column
     xyouts,0.55,0.0,title[0],alignment=0.5,charsize=0.8*csize,$
            color=clr.black, /normal $ ; xtitle
  else $
     xyouts,0.55,0.02,title[0],alignment=0.5,charsize=0.8*csize,$
            color=clr.black, /normal ; xtitle
  xyouts,0.05,0.55,title[1],orientation=90,alignment=0.5,charsize=0.8*csize,$
         color=clr.black,/normal ; ytitle

  if keyword_set(label) then begin
     if size(label,/type) ne 7 then $
        lbl = '!8z!X!Dabs!N = '+string(zabs,format='(f7.5)') $
     else lbl = label
     xyouts,0.55,1.02,lbl,alignment=0.5,charsize=0.8*csize,color=clr.black,$
            /normal
  endif 
  
  ;; Close
  x_psclose
  print,'sdss_velplt: created ',psfil

end                             ; sdss_velplt
