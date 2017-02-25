;+ 
; NAME:
; x_suboscan
;     Version 1.1
;
; PURPOSE:
;     Called by x_mkbias  and x_subbias.
;     This routine does the overscan subtraction using the col and row
;     sections at the right and top, respectively.  The code averages
;     the overscan region (with clipping), identifies bad rows
;     (SVBAD), and then uses SAVGOL to create a smoothed
;     representation of the overscan.  This is then subtracted from
;     the columns in the image.  Finally, the bias row (written at the
;     top of the image) is subtracted from each row after a SAVGOL
;     processing.
;
; CALLING SEQUENCE:
;   
;  x_suboscan, raw, head, ovimg, rbin, cbin, [imtype], /NOBIASROW,
;                  /DEBUG, SVBAD=
;
;
; INPUTS:
;   raw -- 2D image
;   head -- Header for the image (only used to update the card)
;   fincol -- Last column of data section
;   imtype= -- Image type:: 'ZRO'  (only used for debugging)
;
; RETURNS:
;
; OUTPUTS:
;   ovimg -- Bias subtracted image
;
; OPTIONAL KEYWORDS:
;   BIASROW - if set, bias row is used. 
;   SKIPOV  - If set, skip overscan subtraction (bias only)
;   DEBUG   -- Turn debug mode on
;
; OPTIONAL OUTPUTS:
;   SVBAD -- Rows with anomolous behavior, most likely related to a
;            transient occurance within the CCD electronics.
;
; COMMENTS:
;
; REVISION HISTORY:
;   18-July-2003   RAB
;   16-Feb-2004    JXP -- added kludge for bad rows
;
;------------------------------------------------------------------------------

pro  x_suboscan, raw, head, ovimg, fincol, IMTYPE=imtype, $
                 BIASROW=biasrow, DEBUG=debug, SKIPOV=skipov, $
                 SVBAD=svbad, SILENT=silent, FINROW=finrow, $
                 RBIN=rbin, CBIN=cbin

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 4  then begin 
      print,'Syntax:  ' + $
        'x_suboscan, raw, head, ovimg, fincol, [imtype], ' + $
        '/NOBIASROW, CBIN=, RBIN=, SVBAD=, /DEBUG, /SILENT [v1.1])'
      return
  endif 

  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( BADLVL ) then badlvl = 10000.
  if keyword_set( SVBAD ) then delvarx, svbad
  if keyword_set( BIASROW ) AND NOT keyword_set(FINROW) then stop
  if not keyword_set(IMTYPE) then imtype='UNK'

  sz = size(raw, /dimensions)

  ;; overscan subtract from right hand side of chip unless 
  if not keyword_set (SKIPOV) then begin
      if not keyword_set(SILENT) then $
        print, 'x_suboscan: Doing oscan on rows.'
      oscan = raw[fincol+1:sz[0]-1,*]
      ofunct = djs_avsigclip(oscan, 1, sigrej=3, maxiter=3)
      ;; Kludge out bad 'single' rows
      bdrow = where( ofunct GT BADLVL, nbad)
      case nbad of
          0: 
          1: svbad = bdrow
          2: begin
              if bdrow[1] NE bdrow[0] + 1 then svbad = bdrow
          end
          else: begin
              s1 = shift(bdrow,1)
              s2 = shift(bdrow,-1)
              sngl = where( bdrow - (s1+1) NE 0 AND $
                            bdrow - (s2-1) NE 0, nsv)
              if nsv NE 0 then svbad = bdrow[sngl]
          end
      endcase
      if keyword_set( SVBAD ) then begin
          ;; Interpolate
          for i=0L,n_elements(svbad)-1 do begin
              case svbad[i] of
                  0: ofunct[0] = ofunct[1]
                  n_elements(ofunct): ofunct[svbad[i]] = ofunct[svbad[i]-1]
                  else: ofunct[svbad[i]] = $
                    (ofunct[svbad[i]-1]+ofunct[svbad[i]+1]) / 2.
              endcase
          endfor
      endif
      
      ;; Convolve
      savfilt = savgol(60/cbin,60/cbin,0,4 < (12/cbin - 1))
      bias_fct = convol(ofunct, savfilt,/EDGE_TRUNCATE)
                                ; bias_fct_fast = ofunct
      savfilt = savgol(2,2,0,2)
      bias_fct_fast = convol(ofunct, savfilt,/EDGE_TRUNCATE)
      dif = (shift(ofunct,1) - ofunct)
      djs_iterstat, dif , sigma=sig
;      print, 'Doing oscan on rows. Jumps set to (7x sig) = ', 7*sig
      ;; Kludge
      jump = where (abs(dif) gt 7*sig)
      if (n_elements(jump) gt 1) then  begin
          ;; extend masked region  50/rbin to right and 10/rbin to left
          jump2 = lonarr(n_elements(jump) * (72/rbin+1 ))
          for j=0,n_elements(jump)-1 do begin
              for i=0,72/rbin do begin
                  jump2[ j*(72/rbin+1) +i ] = jump[j]+ (i-12/rbin)
              endfor
          endfor
          jumpind = uniq(jump2, sort(jump2))
          jump= jump2[jumpind]
          ;; apply mask
          bias_fct(jump) = bias_fct_fast(jump)
          ;; undefine,jump2 & undefine, jumpind
      endif

      if keyword_set(DEBUG) then begin
          x4plt = findgen(n_elements(bias_fct))
          y = minmax(ofunct)
          x_psopen, 'bias_plots.ps',/maxs
          ;  !p.background = colors.white
          a = raw
          !p.multi=[0,1,2]
          yone= y[0] - (y[1]-y[0])*0.25
          ytwo= y[1] + (y[1]-y[0])*0.25
          plot, x4plt, ofunct, color=colors.black  $
                ,yrange=[yone,ytwo],xrange=[-10,20+n_elements(bias_fct)] $
                ,xstyle=1,ystyle=1 $
                ,xtitle='Row number (y)', ytitle='DN/pixel'
          xyouts, 100, y[1] + (y[1]-y[0])*0.15,'Overscan columns'
          oplot, x4plt, bias_fct,color= colors.cyan
          ; oplot, x4plt, bias_fct_fast,color= colors.red
          stop
          if (imtype eq 'ZRO') then begin 
              temp = djs_avsigclip(raw[0:(2048/cbin)-1,*], 1, sigrej=3, maxiter=3)
              oplot, x4plt, temp,color= colors.red
              if (n_elements(jump) gt 1) then $
                oplot, x4plt(jump), temp(jump),color= colors.yellow, thick=2
              stop
              test_dif = temp - bias_fct
              plot, x4plt, test_dif, color= colors.green
              stop

          endif
      endif
      
      bias_im = replicate(1.,sz[0]) # bias_fct
      ovimg = raw - bias_im
      sxaddpar, head, 'OSCANCOL', 'T'

  endif else ovimg = raw

  ;;  Bias Row off the top, too (if desired)

  if keyword_set( BIASROW ) then begin    ; and from top, if set
      if not keyword_set(SILENT) then $
        print, 'x_suboscan: Doing oscan on columns.'
      oscan = ovimg[*,finrow+2:sz[1]-1]
      ofunct = djs_avsigclip(oscan,2, sigrej=3, maxiter=3)
      savfilt = savgol(36/rbin,36/rbin,0,4 < (12/rbin - 1))
      bias_fct = convol(ofunct,savfilt,/EDGE_TRUNCATE)

      if keyword_set( DEBUG ) then begin
          x4plt = findgen(n_elements(bias_fct))
          x = minmax(ofunct)
          xone= x[0] - (x[1]-x[0])*0.25
          xtwo= x[1] + (x[1]-x[0])*0.25
          plot, findgen(n_elements(ofunct)), ofunct, color=colors.black $
                ,yrange=[xone,xtwo],xrange=[-10,20+n_elements(bias_fct)] $
                ,xstyle=1,ystyle=1 $
                ,xtitle='Column number (x)', ytitle='DN/pixel'
          xyouts, 100, x[1] + (x[1]-x[0])*0.15,'Bias rows'
          oplot, findgen(n_elements(bias_fct)), bias_fct,color= colors.cyan
          stop
;          plot, findgen(n_elements(ofunct)), ofunct, color=colors.black $
;            ,yrange=[xone,xtwo], xrange=[0,500]
;          oplot, findgen(n_elements(bias_fct)), bias_fct,color= colors.cyan
          if (imtype eq 'ZRO') then begin 
              temp = djs_avsigclip(ovimg[*,0:(4096/rbin)], 2, $
                                   sigrej=3, maxiter=3)
              test_dif = temp - bias_fct
              plot, x4plt, test_dif, color= colors.green,yrange=[-1,2],/noerase
          endif
      endif

      bias_im = transpose(replicate(1.,sz[1])#bias_fct)
      ovimg = temporary(ovimg) - bias_im
      sxaddpar, head, 'OSCANROW', 'T'

      if keyword_set( DEBUG ) then begin
          x_psclose
          x_psopen, 'bias_imgs.ps'
          b = ovimg
          c = transpose([transpose(  ((b[1:1024,800:1000]+0.)<200)*(-30.))  $
                         ,transpose( ((a[1:1024,800:1000]-200.)<200)*(-30.))    ])
          c[*,0:1]=0
          c[*,200:202]=0
          c[*,400:401]=0
          c[0:1,*]=0
          c[1022:1023,*]=0
          loadct, 0
          tvscl, c
          x_psclose
      endif

      if keyword_set (DEBUG) then xatv, ovimg, /block
      
  endif

  ;;;;; don't trim the bias image. only trim science exposures when
  ;;;;; bias & overscan correction is applied
  ; ovimg = ovimg[strcol:fincol-1,strrow:finrow-1]
  ; sxaddpar, head, 'TRIM', 'T', imsect

  return
end

