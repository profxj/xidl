;+ 
; NAME:
; esi_lwdtrcarc   
;     Version 1.0
;
; PURPOSE:
;    Trace an Arc Image
;
; CALLING SEQUENCE:
;   
;  esi_lwdtrcarc, esi, slit
;
; INPUTS:
;   esi   -  ESI structure
;   slit  -  Slit size
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdtrcarc, esi, slit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_lwdtrcarc, esi, slit, REFROW=refrow, LINLIST=linlist, $
                   TRCSTR=trcstr

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdtrcarc, esi, slit, LINLIST=, [v1.0]'
      return
  endif 

;  Optional Keywords
  
  if not keyword_set(LINLIST) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/XeHgNe.lst'
  if not keyword_set(REFROW) then refrow = 820L

; Slit
  c_s = esi_slitnm(slit)

; Grab Arc
  arc_fil = 'Arcs/ArcLWD_'+c_s+'.fits'
  print, 'esi_lwdtrcarc: Reading Arc file ', arc_fil
  img_arc = mrdfits(arc_fil, /silent)
  sz_arc = size(img_arc, /dimensions)

  ;; Inverse variance
  ivar = 1./ (img_arc > 1)

  ;; Get Arc Spec + Fit
  arc_spec = djs_median( img_arc[*,refrow-4:refrow+4], 2)
  fitfil = 'Arcs/AFIT_LWD'+c_s+'.fits'
  a = findfile(fitfil, count=na)
  if na EQ 0 then begin
      print, 'esi_lwdtrcarc: No Arc Fit found!'
      stop
  endif
  x_fitstrtofits, fitstr, fitfil, /reverse


  ;; Get starting points
  if not keyword_set( LINTRC ) then begin
      lintrc = [ 4624.28D, $    ; Xe
                 4671.23D, $
                 4734.15D, $
                 7642.02D, $
                 7887.40D, $
                 8231.64D, $
                 8280.12D, $
                 8346.82, $
                 8409.19D, $
                 8819.41D, $
                 8952.25D, $
                 9045.45D, $
                 9162.65D, $
                 4046.563D, $   ; HgNe
                 4077.831D, $
                 4358.327D, $
                 4916.068D, $  
                 5460.735D, $   ; SAT?
                 5769.598D, $  
                 5790.663D, $  
                 5852.488D, $  
                 5881.895D, $ 
                 5944.834D, $    
                 6029.997D, $    
                 6096.163D, $ 
                 6143.062D, $   
                 6217.281D, $
                 6266.495D, $    
                 6304.789D, $    
                 6334.428D, $    
                 6402.246D, $   
                 6506.528D, $   
                 6598.953D, $    
                 6678.276D, $    
                 6717.043D, $    
                 7032.413D, $  
                 7173.938D, $  
                 7245.166D, $   
                 7438.898D $   
               ]
;             8988.58, $    
;             5975.534D, $    
;             8495.359D, $   
;             9148.68D, $    
;             9201.76D $
;             8300.324D, $   
;             8377.606D  $
  endif
  ntrc = n_elements(lintrc)
      
  ;; Get Pix values
  xstart = dblarr(ntrc)
  xval = 475. + findgen(sz_arc[0]-475L+1)
  dum1=findgen(sz_arc[0])
  for ii=0L,ntrc-1 do begin
      fitpix = x_fndfitval(lintrc[ii], fitstr, xval, $
                           TOLER=10.d-6)
      xstart[ii] = x_centspln(dum1[round(fitpix)-3L:round(fitpix)+3L], $
                              arc_spec[round(fitpix)-3L:round(fitpix)+3], $
                              /SILENT)
      if xstart[ii] EQ -1 then xstart[ii] = fitpix
      ;; Center with fweight
      for j=0L,9 do $
        xstart[ii] = trace_fweight(img_arc, xstart[ii], REFROW, radius=1.5, $
                                   invvar=ivar)
  endfor

  ;; Tracing Lines

  print, 'esi_lwdtrcarc: Tracing with trace_crude'
  xcen_pos = trace_crude(img_arc, ivar, yset=ycen_pos, XSTART=xstart, $
                         radius=2., ystart=REFROW, xerr=xerr_pos, $
                         MAXSHIFTE=0.5, NMED=7, NAVE=5)
  ;; Parse out bad lines!
  print, 'esi_lwdtrcarc: Parsing bad trace!'
  gdends = lonarr(ntrc,2)
  ;; Tops
  for q=0L,ntrc-1 do begin
      dx = shift(xcen_pos[1000L:sz_arc[1]-1,q],-20) - xcen_pos[1000L:sz_arc[1]-1,q] 
      a = where(dx GT 0, na)
      if (a[0]+1000L) GT 3070L then gdends[q,1] = 3070L $
      else gdends[q,1] = a[0] + 700L
  endfor
  ;; Bottom
  duml = lindgen(501)
  for q=0L,ntrc-1 do begin
      dx = shift(xcen_pos[0L:500L,q],-20) - xcen_pos[0L:500L,q] 
      a = where(dx LT 0 AND duml LT 475L, na)
      if na EQ 0 then gdends[q,0] = 230L else begin
          if a[na-1] LT 230L then gdends[q,0] = 230L $
          else gdends[q,0] = a[na-1] + 300L
      endelse
  endfor

  ;; Output
  outfil = 'Arcs/ATRC_LWD'+esi_slitnm(slit)+'.fits'
  trcstr = { $
          lintrc: lintrc, $
          ends: gdends, $
          xerr: xerr_pos, $
          xcen: xcen_pos $
  }
  mwrfits, trcstr, outfil, /create, /silent

  ;; DONE
  print, 'esi_lwdtrcarc: All done! Output in ', outfil 

  return
end
