;+ 
; NAME:
; x_fndslitobj   
;    Version 1.0
;
; PURPOSE:
;    Given the slitstr and the map, find slit positions in the
;    original image and create object structure. Used for the WFCCD
;
; CALLING SEQUENCE:
;   
;   x_fndslitobj, img, wvimg, slitstr
;
; INPUTS:
;   img         - Flux image
;   slitstr     - Slit structure
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndslitobj, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fndslitobj, img, wvimg, slitstr, objstr, WVGUESS=wvguess


;  Error catching
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_fndslitobj, img, wvimg, slitstr, objstr, WVGUESS= [v1.0]'
    return
  endif 


;  Optional Keywords

  if not keyword_set( WVGUESS ) then wvguess = 5720.

; Allow images to be fits file

  dat = x_readimg(img, /fscale)
  wav = x_readimg(wvimg, /fscale)
  sz_img = size(dat, /dimensions)

; Slits
  nslit = n_elements(slitstr)
  gdslit = where(slitstr.flg_anly EQ 1, ngdslit)

; Create objects
  tmp = { specobjstrct }
  objstr = replicate(tmp, 5*ngdslit)
  objstr.slit_fil = ' '
  objstr.img_fil = ' '
  objstr.spec2d_fil = ' '
  objstr.UT = ' '
  objstr.instr_strct = ' '
  objstr.field = slitstr[0].field
  objstr.flg_anly = 1

; Counters
  totobj = 0L
  nmissed = 0L

; Loop on Good slits
  for qq=0L,ngdslit-1 do begin
      splog, 'working on slit '+string(qq)+' ('+string(gdslit[qq])+')'

      ; Grab all of the pixels in the slit
      msk1 = bytarr(sz_img[0],sz_img[1])
      msk2 = bytarr(sz_img[0],sz_img[1])
      for i=0L, sz_img[0]-1 do begin
          yedg = round(slitstr[gdslit[qq]].yedg_orig[i,*])
          msk1[i,yedg[0]+1:yedg[1]-1] = 1
          msk2[i,yedg[0]+1:yedg[1]-1] = 1
      endfor
      slitpix = where(msk1 EQ 1)

      ; Go to wave = wvguess in the image
      mn = min( abs(wav[slitpix]-wvguess), imn)
      clm = slitpix[imn] MOD sz_img[0]

      ; Collapse 41 columns 
      cmin = (clm-20) > 0 
      cmax = (clm+20) < (sz_img[0]-1)
      yedg = round(slitstr[gdslit[qq]].yedg_orig[clm,*])
      sub = where( msk2[cmin:cmax, yedg[0]:yedg[1]] EQ 0)
      subimg = dat[cmin:cmax, yedg[0]:yedg[1]]
      subimg[sub] = median(subimg)
      smsh = djs_median( subimg, 1)

      ; Replace zeroed edges
      md = median(smsh)
      bd = where(smsh LT 0.1*md, nbd)
      if nbd NE 0 then smsh[bd] = md

      ; Find peaks
      x_fndpeaks, smsh, center, NSIG=5., FRACPK=0.2, nordb=2, /FORCE, $
        NEDG=2, AFUNC='POLY', EDGES=edges
      flg_pk = 0

      ; Allow for no significant peak
      if center[0] EQ -1 then begin
          print, 'x_fndslitobj: Warning! Slit: '+strtrim(gdslit[qq],2)+$
            ' No peak found, taking supposed value' 
          nmissed = nmissed + 1
          flg_pk = -1
      endif else begin
          percen = center/float(n_elements(smsh))
          ; Find primary object first -- Key on wavelength
          pobj = slitstr[gdslit[qq]].ypos[0]
          mn = min(abs(percen-pobj), imn)
          if mn*n_elements(smsh) GT 5 then begin
              print, 'x_fndslitobj: Warning! Slit: '+strtrim(gdslit[qq],2)+$
                'Obj off by more than 5 pix; Taking original value'
                flg_pk = -1
          endif else slitstr[gdslit[qq]].ypos[0] = percen[imn]
      end

      ; Fill up Obj structure
      x_slittoobj, objstr, totobj, slitstr[gdslit[qq]], clm

      ; Set ypix
      if center[0] NE -1 AND flg_pk NE -1 then begin
          objstr[totobj].ycen = yedg[0] + center[imn] 
          if edges[0,0] NE -1. AND edges[0,1] NE -1. then begin
              objstr[totobj].aper[0] = min(edges[0,*],max=mx) - center[0] - 1.
              objstr[totobj].aper[1] = mx - center[0] + 1.
          endif else objstr[totobj].aper = 1.3*[-1.,1.] ; Arbitrary
      endif else begin
          objstr[totobj].ycen = yedg[0] + $
            slitstr[gdslit[qq]].ypos[0] * (yedg[1]-yedg[0])
          objstr[totobj].aper = 1.3*[-1.,1.]
      endelse
      objstr[totobj].obj_id = 'a'
      totobj = totobj + 1

      ; Set nobj
      slitstr[gdslit[qq]].nobj = 1
      ; Take all other objects
      if n_elements(center) GT 1 then begin
          for i=0L,n_elements(center)-1 do begin
              if i EQ imn then continue
              slitstr[gdslit[qq]].ypos[slitstr[gdslit[qq]].nobj] = center[i]
              ; Obj structure
              x_slittoobj, objstr, totobj, slitstr[gdslit[qq]], clm
              objstr[totobj].ycen = yedg[0] + center[i]
              if edges[i,0] NE -1. AND edges[i,1] NE -1. then begin
                  objstr[totobj].aper[0] = min(edges[i,*],max=mx) - center[i] - 1.
                  objstr[totobj].aper[1] = mx - center[i] + 1.
              endif else objstr[totobj].aper = 1.3*[-1.,1.] ; Arbitrary
              ; ID
              objstr[totobj].obj_id = x_objnumid(slitstr[gdslit[qq]].nobj)
              ; Increment
              slitstr[gdslit[qq]].nobj = slitstr[gdslit[qq]].nobj + 1
              totobj = totobj + 1
              if slitstr[gdslit[qq]].nobj GT 9L then break
          endfor
      endif
  endfor
      
  ; Missed obj
  if nmissed GT 0L then $
  print, 'x_fndslitobj: Missed ', strtrim(nmissed,2), ' objects.'+$
    ' You should run x_setobjgui'

  ; Trim down the object structure
  objstr = temporary(objstr[0:totobj-1])

  ; Set crude trace (not for extraction but ok for masking)
  for i=0L,totobj-1 do begin
      ; Grab the slit
      slit = where(slitstr.id EQ objstr[i].slit_id)
      ; Take nearest edge 
      mn = min(abs(objstr[i].ycen-slitstr[slit].yedg_orig[objstr[i].xcen,*]),imn)
      ; and Offset
      off = objstr[i].ycen - slitstr[slit].yedg_orig[objstr[i].xcen,imn]
      ; Fill it up
      objstr[i].trace[0:sz_img[0]-1] = slitstr[slit].yedg_orig[0:sz_img[0]-1,imn]$
        + off
  endfor

  ; 
  return
end
