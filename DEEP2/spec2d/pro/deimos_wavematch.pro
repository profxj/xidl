;+
; NAME:
;   deimos_wavematch
;
; PURPOSE:
;   Find inital match between found and catalogue lines
;
; CALLING SEQUENCE:
;   deimos_wavematch, model_lambda, spec, lamps,  arcline_x
;
; INPUTS:
;   model_lambda -- optical model structure, lambda(y),y(lambda) (if
;                   0, then brute force solution is sought)
;   spec         -- input spectrum, 1-d
;   lamps        -- structure of lamp wavelengths, intensity
;   chipno       -- which chip number? (1-8)
;   grating      -- grating [lines/mm]
;   lambda_c     -- central wavelength [Ang]
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   slitwid      -- slit width in pixels
;   plot         -- set if you want plots
;
; OUTPUTS:
;   arcline_x    -- pixel positions for found lines 
;   wset - the set of coefficients after cross-correlation
;   
; COMMENTS:
;   We use Drew's optical model, or failing that, a brute force
;   search.
;
;
; REVISION HISTORY:
;
;       Sat Feb 23 11:19:50 2002, Douglas Finkbeiner (dfink)
;		split from 2dtest.pro
;       comments by MD 12Apr02
;
;----------------------------------------------------------------------
     
pro deimos_wavematch, model_lambda, spec, lamps, chipno, grating, lambda_c, $
          arcline_x, slitwid=slitwid, plot=plot, wset=wset1
  
  xtemp = findgen(4096)
  print


; -------- see if we have an optical model
  if keyword_set(model_lambda) then begin
      nomodel = (model_lambda.lambda_y[0] eq 0 $
                 OR model_lambda.lambda_y[1] lt 0) 
      if nomodel then begin
           if (model_lambda.lambda_y_bottom[0] NE 0 AND $
               model_lambda.lambda_y_top[0] NE 0) then $
              model_lambda.lambda_y=( model_lambda.lambda_y_bottom + $
                                      model_lambda.lambda_y_top)/2. $
            else if (model_lambda.lambda_y_bottom[0] NE 0) then $
              model_lambda.lambda_y=model_lambda.lambda_y_bottom $
           else if (model_lambda.lambda_y_top[0] NE 0) then $
              model_lambda.lambda_y=model_lambda.lambda_y_top
      endif    
      nomodel = (model_lambda.lambda_y[0] eq 0) 
  endif else nomodel = 1B

 
;    Use brute force if no model
;     print, 'no optical model available for this slitlet '
;     print, 'USING grating: ', grating
     t1 = systime(1)
     
     color = 'blue'
     if(chipno gt 4) then color = 'red' ;select red or blue!

     if keyword_set(slitwid) then begin 
        smspec = smooth(spec, slitwid)
     endif else smspec = spec

  if nomodel then begin 

     wset1 = discrete_arcfit_guess(smspec, lamps.lambda, grating, lambda_c, $
                                   color=color, plot=plot, func=func )
;                                   arcsat=arcsat)

  endif else begin
      ; use model to get initial lambda, then tweak

	degree=n_elements(model_lambda.lambda_y)
      if degree eq 6 then $
	print,'Optical model guess:',model_lambda.lambda_y,$
	  format='(A21,2f9.3,4f9.5)'

     wset1=discrete_tweak_omodel(smspec, lamps.lambda, $
                 model_lambda.lambda_y, plot=plot, func=func )
  endelse


     if size(wset1, /tname) NE 'STRUCT' then begin 
        message, 'ABORT: No arcfit guess was found!', /info
        arcline_x = 0
        return
     endif 
;     print, 'Tguess: ', systime(1)-t1, '  lambda center', wset1.coeff[0], $
;       format='(A,F6.2,A,F7.1)'

	degree=n_elements(wset1.coeff)
     if degree eq 6 then $
	print,'Discrete-corr result:',wset1.coeff,form='(A21,2f9.3,4f9.5)'

     traceset2xy, wset1, xtemp, lambda

     if keyword_set(plot) then begin
        splot, lambda, spec
        soplot, lamps.lambda, lamps.intensity*4, ps=7
     endif 

     arcline_x = traceset2pix(wset1, lamps.lambda,/silent)



  if keyword_set(pspath) then begin 
;     Plot 4 - wavematch
    
     dfpsplot, pspath+'wavematch.ps', /sq
     plot, lambda, spec, thick=1, xtit='lamgda [Ang]', ytit='counts', $
       title='Wavelength matching', chars=1.5, /xst
     oplot, lamps.lambda, lamps.intensity*4, ps=7, syms=1.5, thick=2
     dfpsclose
  endif 
        

  
  return
end

