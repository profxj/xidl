;+ 
; NAME:
; x_specobjcr   
;   Version 1.1
;
; PURPOSE:
;    Flags CRs given 2 or more object images.  
;
; CALLING SEQUENCE:
;   
;  x_specobjcr, fil, EXPT=, RTIO=, LTHRESH=, GROW=, /NOWRT
;
; INPUTS:
;   fil --  List of filenames
;
; RETURNS:
;
; OUTPUTS:
;   CR removed images.  
;   If /nowrt is set, images are named the same as the original file names, 
;   with 'crrej_' appended to the front.
;   If /nowrt is not set, original images are overwritten.
;
; OPTIONAL KEYWORDS:
;   GROW= -  Value to grow an identified CR by [default: 0]
;   RTIO  -  Ratio CR must exceed to be flagged (default: 9)
;   LTHRESH-
;   NOWRT - Don't overwrite the input files with the CR rejected files. (default: overwrite files)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Algorithm for 3 or more images is not well tested.
;   It would be nice to add the GROW option.
;
; EXAMPLES:
;      x_specobjcr, hires[exp].img_final, EXPT=hires[exp].exp, $
;        RTIO=rtio, LTHRESH=lthresh;, GROW=grow
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;   24-Feb-2005 Ported to XIDL by JXP
;   28-Feb-2007 Minor revisions by LKP
;   10-May-2007 Fixed 'grow' bug JXP
;   17-May-2012 Monkeyed with the headers
;
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_specobjcr_grow, crmap, grow

  ;; Announce
  print, 'x_specobjcr:  Growing CR hits by ', grow

  ;;
  sz = size(crmap, /dimens)
  cridx = where(crmap LT 0., ncr)
  if ncr EQ 0 then return
  xcr = cridx mod sz[0]
  ycr = cridx / sz[0]

  ximg = lindgen(sz[0]) # replicate(1L,sz[1])
  yimg = replicate(1L,sz[1]) # lindgen(sz[1])
  x0 = (xcr - grow) > 0
  x1 = (xcr + grow) < (sz[0]-1)
  y0 = (ycr - grow) > 0
  y1 = (ycr + grow) < (sz[1]-1)

  ;; Loop
  for qq=0L,ncr-1 do begin
      crmap[x0[qq]:x1[qq],y0[qq]:y1[qq]] = -99
  endfor
return
end
  

pro x_specobjcr, fil, EXPT=expt, CHK=chk, RTIO=rtio, LTHRESH=lthresh, $
                 GROW=grow, NOWRT=nowrt, EXT=ext

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'x_specobjcr, fil, EXPT=, LTHRESH=, /CHK, RTIO=, /NOWRT [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 9.
  if not keyword_set(LTHRESH) then lthresh = 80.
  if n_elements(GROW) EQ 0 then grow = 1L  
  if n_elements(EXT) EQ 0 then ext = 0L 
  
  nfil = n_elements(fil)
  if nfil LT 2 then begin
      print, 'x_specobjcr:  Need at least 2 images!'
      return
  endif

  if not keyword_set(EXPT) then expt = replicate(1., n_elements(fil))
  ;; make sure that exposure time array isn't integer type, or else math will be messed up later on.
  expt = float(expt)

  ;;  Read in first image
  if not keyword_set( SILENT ) then print, 'x_specobjcr: Reading images..'

  tmp  = xmrdfits(fil[0], ext, head0, /silent)
  sz   = size(tmp, /dimensions)

  flux = fltarr(sz[0],sz[1],nfil)
  ivar = fltarr(sz[0],sz[1],nfil)
  svexp = fltarr(nfil)

  ;; Header mumbo-jumbo
  all_head = strarr(100000L)
  nlin_head = lonarr(100L)

  nlin_head[0] = n_elements(head0)
  all_head[0:nlin_head[0]-1] = head0

  flux[*,*,0] = temporary(tmp) / expt[0]
  ivar[*,*,0]  = xmrdfits(fil[0], 1, /silent) * (expt[0]^2)
  svexp[0] = expt[0]

  ;; Read in the rest
  for q=1L,nfil-1 do begin
      flux[*,*,q]  = xmrdfits(fil[q], 0, head, /silent) / expt[q]
      ivar[*,*,q]  = xmrdfits(fil[q], 1, /silent) * (expt[q]^2)
      svexp[q] = expt[q]
      ;; Header
      i0 = total(nlin_head)
      nlin_head[q] = n_elements(head)
      all_head[i0:i0+nlin_head[q]-1] = head
  endfor

  if nfil EQ 2 then begin
      diff = (flux[*,*,0] - flux[*,*,1])
      ;; Deal with zero ivar values
      lovar = (ivar[*,*,0]>0.) < (ivar[*,*,1]>0.)
      gd = where(lovar GT 0.)
      hivar = fltarr(sz[0],sz[1])
      snr = hivar
      hivar[gd] = 1./lovar[gd]
;      hivar = (1./(ivar[*,*,0]>0.)) > (1./(ivar[*,*,1]>0.))
      ;; First image
      snr[gd] = diff[gd] / sqrt(hivar[gd])
      cr = where(SNR GT rtio AND flux[*,*,0] GT lthresh/svexp[0], na)
      ;; CR
      crmap = ivar[*,*,0] / (svexp[0]^2)
      if na NE 0 then begin
          crmap[cr] = -99.
          ;; Grow
          if GROW GT 0L then x_specobjcr_grow, crmap, grow
          crmap[where(crmap LT 0.)] = -1.
      endif
      print, 'x_specobjcr:  Identified ', strtrim(na,2), ' cosmic rays in image'+$
        fil[1]
      ivar[*,*,0] = crmap  ;was the intent to change ivar to ivar/(exp^2), because that's what happens here!  (Of course, CR pixels are set to -1.)
      ;; Second Image
      crtot = cr
      cr = where(-1.*SNR GT rtio AND flux[*,*,1] GT lthresh/svexp[1], na)
      crtot = [crtot,cr]
      ;; CR
      crmap = ivar[*,*,1] / (svexp[1]^2)
      if na NE 0 then begin
          crmap[cr] = -99.
          ;; Grow
          if GROW GT 0L then x_specobjcr_grow, crmap, grow
          crmap[where(crmap LT 0.)] = -1.
      endif
      print, 'x_specobjcr:  Identified ', strtrim(na,2), $
        ' cosmic rays in image'+$
        fil[1]
      ivar[*,*,1] = crmap
      if keyword_set( CHK ) then begin
          snr[crtot] = -99.
          xatv, snr, /block, min=-99., max=20.
          print, 'x_specobjcr:  Continue or first set NOWRT=1 to ' + $
            'NOT write the CR subtracted files.'
          stop
      endif
      
  endif else begin
      ;; Median several images
      medimg = djs_median(flux,3)
      print, 'x_specobjcr: Looping!'
      ;; Loop on Images
      crtot = [0L]
      for q=0L,nfil-1 do begin
          tmp = flux[*,*,q]
          rto = tmp/medimg
;              if keyword_set( CHK ) then xatv, rto, /block
          ;; CR
          cr = where(rto GT rtio AND tmp GT lthresh/svexp[q], ncr)
          crtot = [crtot, cr]
          crmap = ivar[*,*,q] / (svexp[q]^2)
          if ncr NE 0 then begin
              crmap[cr] = -99.
              ;; Grow
              if GROW GT 0L then x_specobjcr_grow, crmap, grow
              crmap[where(crmap LT 0.)] = -1.
          endif
          print, 'x_specobjcr:  Identified ', ncr, ' cosmic rays..'
          ;; Set VAR
          ivar[*,*,q] = crmap
      endfor
      if keyword_set( CHK ) then begin
          print, 'x_specobjcr: Check for issues!'
          mdd = median(medimg)
          medimg[crtot] = -99.
          xatv, medimg, /block, min=-3*mdd, max=3*mdd
          print, 'x_specobjcr:  Continue as you wish.. Or set NOWRT=1 '
          stop
      endif
  endelse
  
  ;; OUTPUT
  if not keyword_set( NOWRT ) then begin
     if not keyword_set( SILENT ) then print, 'x_specobjcr: Writing images...'
     for q=0L,nfil-1 do begin
        ;; Note that flux has been altered slightly through the floating point errors associated with multiplying and dividing by exposure time.  
        ;; May want to fix this.
        ;; Header
        if q GT 0 then i0 = total(nlin_head[0:q-1]) else i0 = 0
        head = all_head[i0+lindgen(nlin_head[q])]
        mwrfits, float(flux[*,*,q]*expt[q]), fil[q], head, /create
        mwrfits, float(ivar[*,*,q]), fil[q]
        ;; COMPRESS
        spawn, 'gzip -f '+fil[q]
     endfor
  endif else begin
     if not keyword_set( SILENT ) then print, 'x_specobjcr: Writing images...'
     for q=0L,nfil-1 do begin
        mwrfits, float(flux[*,*,q]*expt[q]), 'crrej_'+fil[q], /create
        mwrfits, float(ivar[*,*,q]), 'crrej_'+fil[q]
        ;; COMPRESS
        spawn, 'gzip -f '+fil[q]
     endfor
  endelse
  if not keyword_set( SILENT ) then print, 'x_specobjcr: All done!'
  return
end
