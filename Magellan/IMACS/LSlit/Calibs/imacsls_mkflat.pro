;+ 
; NAME:
; imacsls_mkflat   
;     Version 1.1
;
; PURPOSE:
;    Identifies flats, processes them.  Creates one flat file
;      per setup.  This is not setup for non-standard binning.
;
; CALLING SEQUENCE:
;  imacsls_mkflat, imacsls, /IFLAT
;
; INPUTS:
;   imacsls  -  IMACS structure
;   setup    -  Setup ID value
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   /REDOOV  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Only set for 1x1 binning
;
; EXAMPLES:
;   imacsls_mkflat, imacsls, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro imacsls_mkflat, imacsls, setup, REDOOV=redoov, SVOV=svov
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'imacsls_mkflat, imacsls, setup, /SVOV, /REDOOV [v1.1]'
      return
  endif 
  
;  Optional Keywords

; Grab all FLAT files

  if not keyword_set( SIDE ) then side = [1L,2L]
  
  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]

      flt = where(imacsls.mode EQ 1 AND imacsls.flg_anly NE 0 AND $
                  imacsls.side EQ qq AND imacsls.setup EQ setup AND $
                  strtrim(imacsls.type,2) EQ 'QTZ' , nflt)
      if nflt EQ 0 then begin
          print, 'imacsls_mkflat: No Flats of type ', 'QTZ', ' found!' 
          continue
      endif


      ;; Bias Subtract
      if not keyword_set( REDOOV ) then $
        bias = where(imacsls[flt].flg_ov EQ 0, nbias) $
      else begin
          nbias = n_elements(flt)
          bias = lindgen(nbias)
      endelse
      if nbias NE 0 then imacsls_subbias, imacsls, flt[bias]


      ;; Median Combine
      statcol = [cline-5,cline+5]
      statrow = round(3000. / imacsls[flt[0]].rbin) + [-15, 15L]

      STATSEC='['+strtrim(statcol[0],2)+':'+strtrim(statcol[1],2)+','+$
        strtrim(statrow[0],2)+':'+strtrim(statrow[1],2)+']'
      xcombine, 'OV/ov_'+imacsls[flt].img_root, img_flat, head, $
        FCOMB=2, SCALE='MED', GAIN=imacsls[flt[0]].gain, $
        RN=imacsls[flt[0]].readno, STATSEC=statsec
      sz = size(img_flat, /dimensions)
    
      ;; Normalize
      medflt = djs_median( img_flat[cline-5:cline+5,*],1)
      ;; Edges
      medflt[0:3] = medflt[3]
      medflt[round(4090./imacsls[flt[0]].rbin):*] = $
        medflt[round(4090./imacsls[flt[0]].rbin)]
      ;; Fit
      bset = bspline_iterfit(findgen(n_elements(medflt)),medflt, yfit=yfit,$
                             everyn=15) 
      nrmflat = img_flat / ( replicate(1., sz[0]) # yfit)

      ;; Output
      outfil = imacsls_getfil('flat_fil', imacsls[flt[0]].setup, SIDE=qq, /name)
      mwrfits, nrmflat, outfil, head, /create, /silent
      spawn, 'gzip -f '+outfil
      print, 'imacsls_mkflat: Flat created ', outfil+'.gz.'

      ;; Set Obj+Std
      objstd = where(imacsls.mode EQ 1 AND imacsls.flg_anly NE 0 AND $
                     imacsls.setup EQ setup AND imacsls.side EQ qq AND $
                     (imacsls.type EQ 'STD' OR imacsls.type EQ 'OBJ'), nobj)
      if nobj NE 0 then imacsls[objstd].flat_fil = outfil

      ;; DEL OV
      if not keyword_set(SVOV) then imacsls_delov, imacsls, flt

  endfor

  print, 'imacsls_mkflat: All done!'
  return
end
