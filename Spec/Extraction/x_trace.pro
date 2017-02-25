;+ 
; NAME:
; x_trace   
;    Version 1.2
;
; PURPOSE:
;    Trace a spectrum. Ok for quick reductions
;
; CALLING SEQUENCE:
;   
;   fit = x_trace(img, aper, cline, PHSIG=, PLSIG=, TFFIT=,
;      NTRC=, PFUNC=, PNORD=, PHSIG=, PLSIG=, NOINTER=, YBUFF=, /PINTER)
;
; INPUTS:
;   img   - 2D spectral image 
;   aper  - Aperture defining the object
;   cline - Column to begin trace at
;
; RETURNS:
;   FIT    - Trace at each column
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NTRC -  Number of columns to median for the trace
;   TFITSTR - Structure for the fitting the trace
;   PFUNC- Profile fitting function
;   PNORD- Profile fitting order
;   PHSIG
;   PLSIG
;   INTER - Interactive fitting for trace
;   
; OPTIONAL OUTPUTS:
;   TFFIT -  Output fit info for the trace
;
; COMMENTS:
;
; EXAMPLES:
;   x_trace, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Nov-2001 Written by JXP
;   02-Feb-2002 JXP -- Added trace_crude option via CRUDE
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_trace, img, guess, cline, VAR=var, TFITSTR=tfitstr, $
                  NTRC=ntrc, PFUNC=pfunc, PNORD=pnord, PHSIG=phsig, $
                  PLSIG=plsig, INTER=inter, YBUFF=ybuff, PINTER=pinter, $
                  CRUDE=crude, ROT=rot, TNORD=tnord, MXSHFT=mxshft
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'fit = x_trace(img, guess, cline, TFITSTR=, /INTER, VAR='
    print, '      NTRC=, PFUNC=, PNORD=, PHSIG=, PLSIG=, YBUFF=, /PINTER '
    print, '      /ROT)  [v1.2]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( NTRC ) then ntrc = 50
  if not keyword_set( YBUFF ) then ybuff = 2L
  if not keyword_set( TNORD ) then tnord = 4

  ; Trace fit structure
  if not keyword_set( TFITSTR ) then begin
      tfitstr = { fitstrct }
      tfitstr.func = 'POLY'
      tfitstr.nord = tnord
      tfitstr.hsig = 3.
      tfitstr.lsig = 3.
      tfitstr.niter = 2
      tfitstr.flg_rej = 1
      tfitstr.maxrej = 10L
      tfitstr.minpt = 10L
  endif

; Variance
  sz = size(img, /dimensions)
  if not keyword_set(VAR) then var = replicate(1., sz[0], sz[1])


;;;;;;;;;;;;;;; My Way (probably not as good as CRUDE) ;;;;;;;

  if not keyword_set( CRUDE ) then begin
      stop
      hfwdth = sz[0]/(2*ntrc)   ; Half width
      totwdth = hfwdth*2 + 1    ; Total width (keep it odd)

; Aperture

      halfaper = (max(aper) - min(aper) + 1)/2.
      cenaper = long(float(aper[0]+aper[1])/2.)


; Setup Trace arrays

      nreg = 1 + $              ; center region
        (cline - hfwdth)/totwdth + $ ; Left
        long( (cline-hfwdth) MOD totwdth NE 0) + $ ; Excess left
        (sz[0] - 1 - cline - hfwdth)/totwdth + $ ; Right
        long( (sz[0]-1-cline-hfwdth) MOD totwdth NE 0) ; Excess right
      
      xpnt = dblarr(nreg)
      ypnt = dblarr(nreg)

; Start at center
      
      ymax = (long(max(aper)) + ybuff) < (sz[1]-1)
      ymin = (long(min(aper)) - ybuff) > 0
      
      xmin = (cline - hfwdth) > 0
      xmax = (cline + hfwdth) < sz[0]-1
      

      xpnt[0] = double(cline)

  ; median to eliminate CRs
      yfit = djs_median(img[xmin:xmax,ymin:ymax],1)
      xfit = findgen(n_elements(yfit)) + float(ymin)
  ; variance :: Assume median improves S/N by sqrt(N)
      vfit = djs_median(var[xmin:xmax,ymin:ymax],1)/float(xmax-xmin+1)

; FIT
      if not keyword_set( PNORD ) then $
        pnord = (long( n_elements(xfit) * 0.75) < 5L)
  ; Make PNORD odd
      if PNORD mod 2 NE 1 then PNORD = PNORD - 1
      if not keyword_set( PFUNC ) then pfunc = 'BSPLIN'

      if keyword_set( PINTER ) then $
        ypnt[0] = x_centerpk( xfit, yfit, FUNC=pfunc, NORD=pnord, HSIG=5., /inter) $
      else $
        ypnt[0] = x_centerpk( xfit, yfit, FUNC=pfunc, NORD=pnord, HSIG=5.)

; Head Left

      qq = 1
      while( xmin GT 0 ) do begin
      ; Set Region limits
          xmin = ( cline - totwdth*qq - hfwdth) > 0
          xmax = ( cline - totwdth*qq + hfwdth) < (sz[0] - 1)
          ymax = long((ypnt[qq-1] + halfaper + ybuff)) < (sz[1]-1)
          ymin = long((ypnt[qq-1] - halfaper - ybuff)) > 0

      ; Set x
          xpnt[qq] = double(xmin+xmax)/2.

      ; Set profile to fit
          yfit = djs_median(img[xmin:xmax,ymin:ymax],1)
          xfit = dindgen(n_elements(yfit)) + float(ymin)
      ; variance :: Assume median improves S/N by sqrt(N)
          vfit = djs_median(var[xmin:xmax,ymin:ymax],1)/float(xmax-xmin+1)

      ; Center on profile
          center = x_centerpk( xfit, yfit, VAR=vfit, FUNC=pfunc, $
                               NORD=pnord, HSIG=5.)
          if center NE -1 then ypnt[qq] = center else begin
              print, 'x_trace: Lost the trace, taking last point'
              ypnt[qq]=ypnt[qq-1]
          endelse
          
      ; Next
          qq = qq + 1
      endwhile

; And now to the right!
      uu = 1
      while( xmax LT sz[0]-1 ) do begin
                                ; Set Region limits
          xmin = ( cline + totwdth*uu - hfwdth) > 0
          xmax = ( cline + totwdth*uu + hfwdth) < (sz[0] - 1)
          if uu NE 1 then begin
              ymax = long((ypnt[qq-1] + halfaper + ybuff)) < (sz[1]-1)
              ymin = long((ypnt[qq-1] - halfaper - ybuff)) > 0
          endif else begin      ; Reset to the center line for the first one
              ymax = long((ypnt[0] + halfaper + ybuff)) < (sz[1]-1)
              ymin = long((ypnt[0] - halfaper - ybuff)) > 0
              flg_first = 1
          endelse

      ; Set x
          xpnt[qq] = double(xmin+xmax)/2.
          
      ; Set profile to fit
          yfit = djs_median(img[xmin:xmax,ymin:ymax],1)
          xfit = dindgen(n_elements(yfit)) + float(ymin)
      ; variance :: Assume median improves S/N by sqrt(N)
          vfit = djs_median(var[xmin:xmax,ymin:ymax],1)/float(xmax-xmin+1)

      ; Center on profile
          center = x_centerpk( xfit, yfit, VAR=vfit, FUNC=pfunc, NORD=pnord, $
                               HSIG=5.)
          if center NE -1 then ypnt[qq] = center else begin
              print, 'x_trace: Lost the trace, taking last point'
              ypnt[qq]=ypnt[qq-1]
          endelse
          
      ; Next
          qq = qq + 1
          uu = uu + 1
      endwhile


  endif else begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CRUDE ;;;;;;;;;;;;;;
;;; CRUDE ;;;;;;;;;;;;;;
;      guess = (aper[0] + aper[1])/2.
      
      if keyword_set( ROT ) then begin
          nwimg = transpose(img)
          nwvar = transpose(1./var)
      endif else begin
          nwimg = img
          nwvar = 1./var
      endelse

      ; Refine the guess (trace crude kludge)
      nwgss = trace_fweight(nwimg, guess, cline, invvar=nwvar)
      nwgss = trace_fweight(nwimg, nwgss, cline, invvar=nwvar)
      nwgss = trace_fweight(nwimg, nwgss, cline, invvar=nwvar)
      
      ; trace crude
      ypnt = trace_crude( nwimg, nwvar, xstart=nwgss, ystart=cline, $
                          xerr=xerr, yset=xpnt, maxshifte=mxshft ) 
  endelse
                              
; Sort
  spnt = sort(xpnt)

; Fit
  fit = x1dfit(xpnt[spnt], ypnt[spnt], INTER=inter, FITSTR=tfitstr)

; Calculate at each column
  if keyword_set( ROT ) then xfit = dindgen(sz[0]) $
  else xfit = dindgen(sz[1])
  
  return, x_calcfit(xfit, FITSTR=tfitstr)

end


