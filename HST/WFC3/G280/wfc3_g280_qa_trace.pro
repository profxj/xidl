;+ 
; NAME:
;  wfc3_g280_qa_trace
;
; PURPOSE:
;   This code generates QA for the tracing algorithm (DEPRICATED).
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_qa_trace, wfc3_g280_strct, ii, specim, BEAM=, QADIR=
;
; INPUTS:
;  wfc3_g280_strct -- the wfc3_g280 structure
;  ii -- the index of the object in the structure
;
; RETURNS:
;
; OUTPUTS:
;  psfile -- QA file of trace
;
; OPTIONAL KEYWORDS:
;  BEAM= which beam to plot
;  QADIR= -- Directory for the QA file
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_qa_trace, wfc3_g280_strct, ii, specim, BEAM=beam, QADIR=qadir
;
; PROCEDURES CALLED:
;  plotting stuff
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Updated to handle structure + multibeam by MN
;------------------------------------------------------------------------------
pro wfc3_g280_qa_trace, wfc3_g280_strct, ii, specim, BEAM=beam, QADIR=qadir

  if (N_params() LT 2) then begin 
     print, 'Syntax - ' + $
            'wfc3_g280_qa_trace, psfile, trace, specim, NAME= [v1.0]'
    return
  endif 

  case beam of
     0: begin
        cnt=wfc3_g280_strct(ii).cnta
        trace_x=wfc3_g280_strct(ii).trace_xa
        trace_y_orig=wfc3_g280_strct(ii).trace_ya_orig
        trace_yfit=wfc3_g280_strct(ii).trace_ya_fit
        trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fita
        bm='A'
     end
     1: begin
        cnt=wfc3_g280_strct(ii).cntc
        trace_x=wfc3_g280_strct(ii).trace_xc
        trace_y_orig=wfc3_g280_strct(ii).trace_yc_orig
        trace_yfit=wfc3_g280_strct(ii).trace_yc_fit
        trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fitc
        bm='C'       
     end
     else: stop
  endcase
  
  ;; PSFILE
  psfile=qadir+wfc3_g280_strct(ii).name+'_trace'+bm+'.ps'
  x_psopen, psfile, /maxs

  ;; Image
  dim=size(specim,/dim)
  tmpgd=where(trace_x gt 0 and trace_x lt dim(0))
  gd=[min(tmpgd):max(tmpgd)]
  xmn = min(trace_x(gd), max=xmx)
  xmnx = round([(dim(0)-1)<xmn>0L,(dim(0)-1)<xmx>0L])
  ymn = min(trace_y_orig(gd), max=ymx)
  YBUFF = 10L
  ymn = (dim(1)-1)<(ymn-YBUFF*3)>0L 
  ymx = (dim(1)-1)<(ymx+YBUFF)>0L
  ymnx = round([ymn,ymx])
  img_cut = specim[xmnx[0]:xmnx[1], ymnx[0]:ymnx[1]]
  szc = size(img_cut,/dimen)

  devicefactor=2540. ;; cm to inch
  imsize = 9.0       ;; inch
  irange = [-10, 300.]
  
  xpos1 = 1.0
  ypos1 = 2.
  
  ;; Linear
  if not keyword_set(XNCOLORS) then xncolors=200L
  scaled = bytscl(img_cut, min=irange[0], max=irange[1], $
                  top=(xncolors - 1))
     
  ;; Plot
  ctload, 0, ncolors=xncolor, /rever
  tv, scaled, xpos1, ypos1, xsize=imsize, /inches

  ;; Axes
  dims = size(scaled,/dim)
  xlabel = findgen(dims[0]+1)
  ylabel = findgen(dims[1]+1)
  thisPosition = devicefactor*[xpos1, ypos1, $
                               xpos1+imsize, $
                               ypos1+(imsize*dims[1]/dims[0])]
  ctload, 0, ncolors=xncolor
  plot, xlabel, ylabel, /nodata, /device, /noerase, $
        position=thisPosition, $
        xrange=[min(xlabel),max(xlabel)], $
        yrange=[min(ylabel),max(ylabel)], $
        xstyle=5, ystyle=5
  plot, [0], [0], /device, /noerase, xrange=xmnx, $
        yrange=ymnx, xtitle='Column', charsiz=csz2, ytickn=yspaces, $
        ytitle='Row', /nodata, xsty=1, ysty=1, ytickint=10, $
        position=thisPosition ;, ytickname=['-2','-1','0','+1','+2']
  
  clr = getcolor(/load)
  oplot, trace_x(gd), trace_y_orig(gd), color=clr.red, thick=1, linest=1
  oplot, trace_x(gd), trace_yfit(gd), color=clr.cyan, thick=1, linest=2
  oplot, trace_x(gd), trace_yfit(gd)+trace_sigma_fit(gd), color=clr.yellow, $
         thick=1, linest=1
  oplot, trace_x(gd), trace_yfit(gd)-trace_sigma_fit(gd), color=clr.yellow, $
         thick=1, linest=1
  oplot, trace_x(gd), trace_yfit(gd)+6*trace_sigma_fit(gd), color=clr.blue, $
         thick=1, linest=1
  oplot, trace_x(gd), trace_yfit(gd)-6*trace_sigma_fit(gd), color=clr.blue, $
         thick=1, linest=1

  ;; Label
  xyouts, 0.5, 0.6, wfc3_g280_strct(ii).name+' Trace of Beam '+bm, color=clr.black, $
          charsiz=2., /norma, align=0.5 
  x_psclose

  return
end
