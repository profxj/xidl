;+ 
; NAME:
;  wfc3_g280_qa_position
;
; PURPOSE:
;   Creates a QA file of the position (DEPRICATED)
;
; CALLING SEQUENCE:
;   
;   wfc3_g280_qa_position, wfc3_g280_strct, SRCH=srch, QADIR=qadir
;
; INPUTS:
;
;   wfc3_g280_strct -- the wfc3_g280 structure
; 
; RETURNS:
;
; OUTPUTS:
;   Creates a QA file for the requested object of the position
;
; OPTIONAL KEYWORDS:
;  SRCH= -- Size of the search box
;  QADIR= -- Directory for the QA file
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; wfc3_g280_qa_position, wfc3_g280_strct, SRCH=srch, QADIR=qadir
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_qa_position, wfc3_g280_strct, SRCH=srch, QADIR=qadir


  for ii=0L,n_elements(wfc3_g280_strct)-1 do begin

     ;; naming the file
     psfile = QADIR+wfc3_g280_strct(ii).name+'_position.eps'

     ;; open the file
     ps_start, psfile, /nomatch, /encapsulated, xs=4, ys=4
     
     ;; readfile
     ext=7-3*wfc3_g280_strct(ii).chip
     img = xmrdfits(wfc3_g280_strct[0].img_fil,ext,head)
     
     ;; trim image
     buff=10L
     offset=srch+buff

     x0=round(wfc3_g280_strct(ii).x0)
     y0=round(wfc3_g280_strct(ii).y0)
     img_trim=img[x0-offset:x0+offset,y0-offset:y0+offset]

     ;; the plotting
     pos=[0.12,0.12,0.9,0.9]
     cgimage, img_trim, ctindex=0, pos=pos, $
              minval=-10, maxval=200
     cgplot, [0], [0], xr=[-1*offset-0.5,offset+0.5], yr=[-1*offset-0.5,offset+0.5], $
             /nodata, /noerase, pos=pos, xs=1, ys=1, $
             xtit='X (pixels)', ytit='Y (pixels)', xcharsize=0.5, $
             ycharsize=0.5
     cgplot, [wfc3_g280_strct(ii).x0-x0],[wfc3_g280_strct(ii).y0-y0], $
             psym=16, color='blue', /ov
     cgplot, [wfc3_g280_strct(ii).x0-x0],[wfc3_g280_strct(ii).y0-y0], $
             psym=9, color='blue', /ov, symsize=3
     cgplot, [wfc3_g280_strct(ii).xguess-wfc3_g280_strct(ii).x0], $
             [wfc3_g280_strct(ii).yguess-wfc3_g280_strct(ii).y0], $
             psym=16, color='red', /ov
     
     cgtext, -1.*offset, 1.1*offset, 'Pixel centroid = ('+ $
             strcompress(string(round(wfc3_g280_strct(ii).x0)),/r)+','+$
             strcompress(string(round(wfc3_g280_strct(ii).y0)),/r)+')'
     
     ps_end

  endfor   
end
