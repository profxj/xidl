;+ 
; NAME:
; mike_mkbias   
;     Version 1.1
;
; PURPOSE:
;    Create bias frames from zero exposures.
;    This program is not really necessary as it is VERY unlikely that
;    a bias frame will be used in future reductions.  Overscan on 
;    top and right of frames is sufficient. Combined bias frames
;    created here can be used to verify that this is the case.
;
; CALLING SEQUENCE:
;   
;  mike_mkbias, mike, /NOBIASROW, /CLOBBER, /DEBUG, /SVOV
;
; INPUTS:
;   mike   -  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;   Creates bias frames in 'Bias/Bias[NxN][B,R].fits'
;   where [NxN] indicates binning.
;
; OPTIONAL KEYWORDS:
;    CLOBBER = Overwrite old bias frames if set
;    NOBIASROW - if set, bias row is ***NOT*** used. 
;             I.E,  default is to use the bias row.
;    SVOV -- Save the bias subtracted frames of the individual files
;    DEBUG -- Turn on debug flag
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mkbias, mike, /CLOBBER
;
; PROCEDURES/FUNCTIONS CALLED:
;   mike_oscansub
; 
; REVISION HISTORY:
;   28-Feb-2002 Written by JXP
;   July-2003   Method for bias subtraction changed by RAB. 
;               Now:
;               Overscan subtract both red and blue  0 second exposures.
;               Combine to make a bias frames.
;------------------------------------------------------------------------------

pro mike_mkbias, mike, NoBIASROW= nobiasrow, CLOBBER=clobber, DEBUG=debug, $
                 SVOV=svov

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_mkbias, mike, /NoBIASROW, /CLOBBER, /DEBUG, /SVOV [v1.1]'
;      print,'Recommended - ' + $
;        'mike_mkbias, mike, /NoBIASROW, /CLOBBER ' 
      return
  endif 
  
;  Optional Keywords

; Organize all bias frames

  bintot = mike.colbin + 10*mike.rowbin
  bintyp = bintot[uniq(bintot, sort(bintot))]
  nbin = n_elements(bintyp)

; Loop on Side
  for qq=1L,2 do begin  ; 1 = blue
      if qq EQ 1 then nm = 'B' else nm = 'R' 
      if qq EQ 1 then print, 'mike_mkbias: Creating BLUE bias frames' $
      else print, 'mike_mkbias: Creating RED bias frames'

      ;; Loop on binning
      for i=0L,nbin-1 do begin
          row = bintyp[i]/10
          col = bintyp[i] - row*10
          print, 'mike_mkbias: Creating bias for rowbin = ', row, $
            ' and col bin = ',$
            col
          ;;  Outfil
          outfil = 'Bias/Bias'+strtrim(col,2)+'x'+strtrim(row,2)+nm+'.fits'
          a = findfile(outfil+'*', count=na)
          if na NE 0 and not keyword_set( CLOBBER ) then  begin
              print, 'mike_mkbias: Outfil ', outfil, ' already exists!' 
              continue
          endif
          ;; Grab all the ZRO frames
          zro = where(mike.type EQ 'ZRO' AND mike.rowbin EQ row AND $
                      mike.colbin EQ col and mike.flg_anly NE 0 AND $
                      mike.side EQ qq, nzro)
          if nzro EQ 0 then begin
              print, 'mike_mkbias: No ZRO frames. Just use overscan regions.' 
              continue
          endif

          ovroot = 'OV/ov_'

          for q=0,n_elements(zro)-1 do begin
              print, ' ' 
              print, 'mike_mkbias: Subtracting the overscan for image ', $
                mike[zro[q]].rootpth + mike[zro[q]].img_root, mike[zro[q]].type

              oscanfil = ovroot + mike[zro[q]].img_root
              raw = xmrdfits(mike[zro[q]].rootpth + mike[zro[q]].img_root $
                            , 0, head, /fscale, /silent)
              mike_suboscan, raw, head, ovimg, $
                mike[zro[q]].rowbin, mike[zro[q]].colbin, mike[zro[q]].type $
                ,NoBIASROW= nobiasrow, DEBUG= debug
              mwrfits, ovimg, oscanfil, head, /create, /silent
          endfor

          ;; Combine
          
          if nzro GT 1 then $
            xcombine, ovroot + mike[zro].img_root, $
            img_zro, head, FCOMB=0 $
          else img_zro = xmrdfits(ovroot + mike[zro].img_root, $
                                  /fscale, /silent)
          ;; Output
          mwrfits, img_zro, outfil, /create, /silent
          spawn, 'gzip -f '+outfil
          print, 'mike_mkbias Bias file created ', outfil 

          ;; Delete ov files
          if not keyword_set( SVOV ) then mike_delov, mike, zro, /silent
      endfor
  endfor

  print,' *********************************************** '
  print, ' '
  print,'Bias frames created in Bias/ .'
  print,'Set the flag /usebias in mike_subbias to use them. '

  return
end

  
