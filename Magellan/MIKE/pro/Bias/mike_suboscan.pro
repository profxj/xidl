;+ 
; NAME:
; mike_suboscan
;     Version 1.1
;
; PURPOSE:
;     Called by mike_mkbias  and mike_subbias.
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
;  mike_suboscan, raw, head, ovimg, rbin, cbin, [imtype], /NOBIASROW,
;                  /DEBUG, /REDARC, SVBAD=
;
;
; INPUTS:
;   raw -- 2D image
;   head -- Header for the image (only used to update the card)
;   rbin -- Binning of the rows (1,2,3,4)
;   cbin -- Binning of the columns (1,2,3,4)
;   [imtype] -- Image type:: 'ZRO'  (only used for debugging)
;
; RETURNS:
;
; OUTPUTS:
;   ovimg -- Bias subtracted image
;
; OPTIONAL KEYWORDS:
;   NoBIASROW - if set, bias row is ***NOT*** used. 
;             I.E,  default is to use the bias row.
;   REDARC  -- Because of charge bleeding, the red arcs tend to have
;              improper overscan regions.  Setting this keyword will
;              stop the code from using the overscan for bias
;              subtraction.
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

pro  mike_suboscan, raw, head, ovimg, rbin, cbin, imtype $
                    , NoBIASROW=nobiasrow, DEBUG=debug, REDARC=redarc, $
                    SVBAD=svbad, SILENT=silent

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 5  then begin 
      print,'Syntax:  ' + $
        'mike_suboscan, raw, head, ovimg, rbin, cbin, [imtype], ' + $
        '/NOBIASROW, /REDARC, SVBAD=, /DEBUG, /SILENT [v1.1])'
      return
  endif 

  if not keyword_set( BADLVL ) then badlvl = 10000.
  if keyword_set( SVBAD ) then delvarx, svbad
  if not keyword_set(IMTYPE) then imtype='UNK'

  sz = size(raw, /dimensions)

  strcol = 0L 
  fincol = (2048L / cbin) - 1
  strrow = 0L 
  finrow = (4096L / rbin) - 1 

  if keyword_set( DEBUG ) then begin
      print,'strrow: ', strrow ,' finrow: ', finrow
      print,'strcol: ', strcol ,' fincol: ', fincol
      stop
  endif

  print, '  '

  imsect='['+string(fincol+1)+':'+string(sz[0]-1)+',*/' 

  ; overscan subtract from right hand side of chip unless /REDARC
  if not keyword_set (redarc) then begin
      if not keyword_set(SILENT) then $
        print, 'mike_suboscan: Doing oscan on rows.'
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
          yone= y[0] - (y[1]-y[0])*0.5
          ytwo= y[1] + (y[1]-y[0])*0.5
          plot, x4plt, ofunct, color=colors.white  $
            ,yrange=[yone,ytwo],xrange=[-10,20+n_elements(bias_fct)],xstyle=1,ystyle=1
          oplot, x4plt, bias_fct,color= colors.cyan
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
      ovimg = temporary(raw) - bias_im
      sxaddpar, head, 'OSCANCOL', 'T'

  endif else ovimg = temporary(raw)

  ; if keyword_set(debug) then xatv, ovimg, /block
      

;  now Bias Row off the top, too

  if not keyword_set( NoBIASROW ) then begin    ; and from top, if set
      if not keyword_set(SILENT) then $
        print, 'mike_suboscan: Doing oscan on columns.'
      oscan = ovimg[*,finrow+2:sz[1]-1]
      ofunct = djs_avsigclip(oscan,2, sigrej=3, maxiter=3)
      savfilt = savgol(36/rbin,36/rbin,0,4 < (12/rbin - 1))
      bias_fct = convol(ofunct,savfilt,/EDGE_TRUNCATE)

      if keyword_set( DEBUG ) then begin
          x4plt = findgen(n_elements(bias_fct))
          x = minmax(ofunct)
          xone= x[0] - (x[1]-x[0])*0.5
          xtwo= x[1] + (x[1]-x[0])*0.5
          plot, findgen(n_elements(ofunct)), ofunct, color=colors.white $
            ,yrange=[xone,xtwo]
          plot, findgen(n_elements(ofunct)), ofunct, color=colors.white $
            ,yrange=[xone,xtwo]
          oplot, findgen(n_elements(bias_fct)), bias_fct,color= colors.cyan
          stop
          plot, findgen(n_elements(ofunct)), ofunct, color=colors.white $
            ,yrange=[xone,xtwo], xrange=[0,500]
          oplot, findgen(n_elements(bias_fct)), bias_fct,color= colors.cyan
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

      ;if keyword_set (DEBUG) then xatv, ovimg, /block
      
  endif

  ;;;;; don't trim the bias image. only trim science exposures when
  ;;;;; bias & overscan correction is applied
  ; ovimg = ovimg[strcol:fincol-1,strrow:finrow-1]
  ; sxaddpar, head, 'TRIM', 'T', imsect

  return
end

