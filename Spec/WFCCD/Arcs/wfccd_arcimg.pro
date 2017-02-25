;+ 
; NAME:
; wfccd_arcimg
;    Version 1.0
;
; PURPOSE:
;    Creates an arc image given the wfccd structure
;
; CALLING SEQUENCE:
;   
;   wfccd_arcimg, wfccdstr, WFAIMG=
;
; INPUTS:
;   wfccdstr    - WFCCD structure
;   mskid     - Long defining the mask to process
;
; RETURNS:
;
; OUTPUTS:
;   wfaimg      -  WFCCD arc image (fits file)
;
; OPTIONAL KEYWORDS:
;  SLITSTR   - slit structure
;  WFASTR    - WFCCD arc structure
;  NOFITS   - Suppress fits output
;  OUTFIL   - Suppress fits output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_arcimg, wfccdstr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Feb-2002 Written by JXP
;   16-Jul-2002 Modified: Dealt with edge of CCD
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  MAIN DRIVER
pro wfccd_arcimg, wfccd, mask_id, exp_id, WFAIMG=wfaimg, NOFITS=nofits, $
                  OUTFIL=outfil, MAP=map, CLOBBER=clobber, $
                  WFSLIT=wfslit, WFASTR=wfastr, WFARC=wfarc, $
                  NOVAC=novac

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_arcimg, wfccd, mask_id, [exp_id], WFAIMG=, WFARC=, SLITSTR=, '
    print, '          OUTFIL=, /NOFITS, /NOVAC [v1.0]'
    return
  endif 

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

;  Optional Keywords
  if not keyword_set(OUTFIL) then begin
      i = strpos(wfccd[exp].arc_fil, '_')
      outfil = strtrim('Arcs/ArcI'+strmid(wfccd[exp].arc_fil,i),2)
  endif

;;;;;;;;;;;;;;;;;;;;;;;
;;  Check to see if the Image exists

  if not keyword_set( CLOBBER ) then begin
      a = findfile(outfil, count=count)
      if count GT 0 then begin
          print, 'wfccd_arcimg: Arc image (', outfil, ') exists!'
          return
      endif
  endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SET UP THE STRUCTURES AND IMAGES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read in Arc structure
  if not keyword_set( WFASTR ) then begin
      i = strpos(wfccd[exp].arc_fil, '_')
      fil = 'Arcs/ArcS'+strmid(wfccd[exp].arc_fil,i)
      wfccd_readastrct, strtrim(fil,2), wfastr
  endif

; Read in Slit structure
  if not keyword_set( WFSLIT ) then $
    wfslit = xmrdfits(wfccd[exp].slit_fil,1,STRUCTYP='mslitstrct', /silent)
  
; Read in Arc image
  if not keyword_set( WFARC ) then wfarc = xmrdfits(wfccd[exp].arc_fil, /silent)

; Read in map
  if not keyword_set( MAP ) then map = xmrdfits(wfccd[exp].map_fil,0, /silent)


;;;;
;  Read in the LINE LIST
  linelist = $
    getenv('XIDL_DIR')+'/Spec/Arcs/Lists/wfccdB_HeNe.lst'
  x_arclist, linelist, lines

;;;;;;;;;;;;;
; Create output image
  sz_arc = size(wfarc, /dimensions)
  wfaimg = fltarr(sz_arc[0], sz_arc[1])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  LOOP ON SLITS
;;;;;;;;
  ; Find good exp
  gdexp = where(wfastr.flg_anly NE 0, ngd)

  if not keyword_set( SILENT ) then $
    print, 'wfccd_arcimg: Looping on the slits...'

  for i=0L,ngd-1 do begin
      print, 'wfccd_arcimg: Slit '+strtrim(i,2)+' of '+strtrim(ngd-1,2)
      ; Edges
      edgs = round(wfslit[gdexp[i]].yedg_flt)
      ;  Trace Structure
      trcstr = x_tracearc(wfarc, edgs, NSIG=5., RN=5.6, RADIUS=1.5, /SILENT, $
                          YSTRT=wfastr[gdexp[i]].cent) 
      ;  Create Image
      tmp_sz = [sz_arc[0], edgs[1]-edgs[0]+1]
      tmpaimg = $
        x_arcimage(trcstr, wfastr[gdexp[i]].fit, tmp_sz, lines, NSIG=3., $
                  YSTRT=wfastr[gdexp[i]].cent-edgs[0])
      ; VACUUM CORRECTION
      a = where(tmpaimg GT 0.)
      gdaimg = tmpaimg[a]
      if not keyword_set(NOVAC) then begin
          airtovac, gdaimg
          tmpaimg[a] = gdaimg
      endif
      ; SAVE
      wfaimg[*,edgs[0]:edgs[1]] = tmpaimg
  endfor

;;;;;;;
; Reverse Map
;;;;;;;

  wfaimg = x_invertarc(wfaimg, map)

; OUTPUT

  if not keyword_set( NOFITS ) then begin
      mwrfits, wfaimg, outfil, /create, /silent
      spawn, 'gzip -f '+outfil
  endif

  print, 'wfccd_arcimg: All done!'
  return
end
