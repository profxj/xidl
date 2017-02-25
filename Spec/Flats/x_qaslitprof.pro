;+ 
; NAME:
; x_qaslitprof
;     Version 1.1
;
; PURPOSE:
;   Creates the QA files for the x_slitflat routine
;    
; CALLING SEQUENCE:
;  qa_profile, qastr, qafil
;
; INPUTS:
;  qastr -- The QA structure
;  qafil -- The name of the QA file
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2004 Written by SB
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_qaslitprof, qastr, qafil
      

   norder = n_elements(qastr)

   x_psopen, qafil, /maxs
   clr = getcolor(/load)
   !p.multi=[0,3,1]

   for ii = 0,norder-1 do begin
     e = [min(qastr[ii].cen), max(qastr[ii].cen)]
     plot, e, [-0.25,1.15], /nodata, /xs,/ys, xmar=[0,0], ychars=2., xchars=2.
     xyouts, [-0.1], [0.8], 'Order #'+string(qastr[ii].order,format='(i4)')
     xyouts, [-0.2], [0.75], 'Npixels fit'+string(qastr[ii].npix, format='(i6)')
     xyouts, [-0.2], [0.70], 'Npixels rej'+string(qastr[ii].nrej, format='(i6)')
     xyouts, [-0.1], [0.65], 'Chi^2 '+string(qastr[ii].chi2, format='(f7.2)')

;     oplot, qastr[ii].cen, qastr[ii].median, ps=4
     err = (qastr[ii].p95 - qastr[ii].p05)/2.
     y = (qastr[ii].p95 + qastr[ii].p05)/2.
     djs_oploterr, qastr[ii].cen, y, yerr=err, xlen=0., thick=0.5
     djs_oploterr, qastr[ii].cen, (y-qastr[ii].mid)-0.1, yerr=err, $
             xlen=0., thick=0.5
     oplot, qastr[ii].cen, qastr[ii].mid, ps=10, thick=0.5
     oplot, qastr[ii].cen, qastr[ii].edgel, color=clr.blue, thick=0.3
     oplot, qastr[ii].cen, qastr[ii].edger, color=clr.red, thick=0.3
     oplot, e, poly(e,qastr[ii].ab), thick=2
   endfor


   x_psclose
   !p.multi=[0,1,1]

   replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
   ps_replacetitle, replace_title, qafil

   spawn, 'gzip -f '+qafil

   
return
end


