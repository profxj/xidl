;+
; NAME:
;   DEIMOS_FIXGLOW
;
; PURPOSE:
;   To remove the excess counts at the edges of flats (or other frames
;   with lots of counts) coming from the nearby silicon, evidently.
;
; CATEGORY:
;   Spec2d
;
; CALLING SEQUENCE:
;    fixedflat = DEIMOS_FIXGLOW(flatimage [,nfix=nfix,mask=mask,/FIXBOTTOM]
;
; INPUTS:
;    flatimage - a standard, 2kx4k image of one chip's flatfield 
;
; KEYWORD PARAMETERS:
;    nfix=nfix - set this to alter the number of rows/columns
;                corrected (default 5)
;    mask=mask - an array the size of flatimage, 1 on bad pixels.
;                 if mask is set, djs_maskinterp will be applied before medianing,
;                 and no correction will be performed on masked
;                 pixels.
;    /FIXBOTTOM:  - set this to fix the bottom rows of the chip, rather
;                than the top.  This option should be set for chips 5-8
; OUTPUTS:
;    fixedflat - the repaired image.  nfix rows/columns are adjusted
;                (multiplicatively - this was much better than the 
;                subtractive version)
;                so that their local median level
;                matches the nfix rows/columns just interior to the fixed
;                region.  Note that the bottom of the image is not
;                repaired, as it does not exhibit a glow.
;
; MODIFICATION HISTORY:
;     created by JAN 8/2/02
;-

function deimos_fixglow, flatimage, nfix=nfix,mask=mask,fixbottom=fixbottom


  if n_elements(fixbottom) eq 0 then fixbottom=0
  if n_elements(nfix) eq 0 then nfix = 5

  fixedimage = flatimage 
  if n_elements(mask) eq 0 then mask=flatimage*0

  ; Fix top
  if fixbottom eq 0 then begin 
	  tempmask=reform(mask[*, 4096-nfix-1])
	  range = minmax(where(tempmask eq 0))
	  range[0] = range[0] > nfix
	  range[1] = range[1] < 2048-nfix-1
	  useidx = lindgen(range[1]-range[0]+1)+range[0]
	
	  tempmask=reform(mask[useidx, 4096-nfix-1])
	  tempmodel=reform(flatimage[useidx, 4096-nfix-1])
	  tempmodel=djs_maskinterp(tempmodel,tempmask)
	  topmodel = djs_median( tempmodel, width=201, boundary='reflect')
	
	  for i=1, nfix do begin
	     tempmodel=reform(flatimage[useidx, 4096-i])
	     tempmask=reform(mask[useidx, 4096-i])
	     tempmodel=djs_maskinterp(tempmodel,tempmask)
	     rowmodel = djs_median( tempmodel, width=201, boundary='reflect')
	     fixedimage[useidx, 4096-i] = flatimage[useidx, 4096-i]*(topmodel/rowmodel)
	  endfor
	endif else begin
            tempmask=reform(mask[*, nfix])
	  range = minmax(where(tempmask eq 0))
	  range[0] = range[0] > nfix
	  range[1] = range[1] < 2048-nfix-1
	  useidx = lindgen(range[1]-range[0]+1)+range[0]
	
	  tempmask=reform(mask[useidx, nfix])
	  tempmodel=reform(flatimage[useidx, nfix])
	  tempmodel=djs_maskinterp(tempmodel,tempmask)
	  topmodel = djs_median( tempmodel, width=201, boundary='reflect')
	
	  for i=0, nfix-1 do begin
	     tempmodel=reform(flatimage[useidx, i])
	     tempmask=reform(mask[useidx, i])
	     tempmodel=djs_maskinterp(tempmodel,tempmask)
	     rowmodel = djs_median( tempmodel, width=201, boundary='reflect')
	     fixedimage[useidx, i] = flatimage[useidx, i]*(topmodel/rowmodel)
	  endfor
	endelse


  ; Fix right side
  ; need to match gradient as well as normal level due to vignetting

  tempmodarr=fltarr(4096,nfix)
  rawarr=tempmodarr
  tempmodcol=fltarr(nfix)
  rawcol=fltarr(nfix)

  range = [0, 4095]
  fullmask = fltarr(4096)+1

  for i=0,nfix-1 do begin
     tempmask=reform(mask[2048-nfix-1-i, *])
     temprange = minmax(where(tempmask eq 0))
     range[0] = range[0] > temprange[0]
     range[1] = range[1] < temprange[1]
  endfor
     useidx = lindgen(range[1]-range[0]+1)+range[0]
     fullmask[useidx] = 0
     nuse = n_elements(useidx)

 tempmodarr=fltarr(nuse,nfix)
  rawarr=tempmodarr
  tempmodcol=fltarr(nfix)
  rawcol=fltarr(nfix)


  for i=0,nfix-1 do begin
	tempmodel=reform(flatimage[2048-nfix-1-i, useidx])
	tempmodcol[i]=2048-nfix-1-i
	tempmask=reform(mask[2048-nfix-1-i, useidx])
	tempmodel=djs_maskinterp(tempmodel,tempmask)
  	tempmodarr[*,i] = djs_median( tempmodel, width=201, boundary='reflect')
  endfor

	; linear least-squares fit to the model region
	N=nfix
	SS=N
	x=tempmodcol
	sumx=total(tempmodcol)
	t=x-sumx/ss
	tarr=(fltarr(nuse)+1) # t
	sumtsq=total(t^2)
	slopearr=total(tarr*tempmodarr,2)/sumtsq
;	slopearr=djs_median(slopearr,width=101,boundary='reflect')
	constarr=(total(tempmodarr,2)-slopearr*total(x))/ss
;	constarr=djs_median(constarr,width=101,boundary='reflect')

  for i=1, nfix do begin
	  tempmodel=reform(flatimage[2048-i, useidx])
	  rawcol=2048-i
	  tempmask=reform(mask[2048-i, useidx])
	  tempmodel=djs_maskinterp(tempmodel,tempmask)
	  colmodel = djs_median( tempmodel, width=201, boundary='reflect')
	  rightmodel=djs_median(constarr+slopearr*rawcol,width=201,boundary='reflect')
          fixedimage[2048-i, useidx] = flatimage[2048-i,useidx]*(rightmodel/colmodel)
  endfor
	
  ; Fix left side


  range = [0, 4095]
  fullmask = fltarr(4096)+1

  for i=0,nfix-1 do begin
     tempmask=reform(mask[nfix+i, *])
     temprange = minmax(where(tempmask eq 0))
     range[0] = range[0] > temprange[0]
     range[1] = range[1] < temprange[1]
  endfor
     useidx = lindgen(range[1]-range[0]+1)+range[0]
     fullmask[useidx] = 0
     nuse = n_elements(useidx)

 tempmodarr=fltarr(nuse,nfix)
  rawarr=tempmodarr
  tempmodcol=fltarr(nfix)
  rawcol=fltarr(nfix)



  for i=0,nfix-1 do begin
	tempmodel=reform(flatimage[nfix+i, useidx])
	tempmodcol[i]=nfix+i
	tempmask=reform(mask[nfix+i, useidx])
	tempmodel=djs_maskinterp(tempmodel,tempmask)
  	tempmodarr[*,i] = djs_median( tempmodel, width=201, boundary='reflect')
  endfor
 
	; linear least-squares fit to the model region
	N=nfix
	SS=N
	x=tempmodcol
	sumx=total(tempmodcol)
	t=x-sumx/ss
	tarr=(fltarr(nuse)+1) # t
	sumtsq=total(t^2)
	slopearr=total(tarr*tempmodarr,2)/sumtsq
;	slopearr=djs_median(slopearr,width=101,boundary='reflect')
	constarr=(total(tempmodarr,2)-slopearr*total(x))/ss
;	constarr=djs_median(constarr,width=101,boundary='reflect')

  for i=0, nfix-1 do begin
	  tempmodel=reform(flatimage[i, useidx])
	  rawcol=i
	  tempmask=reform(mask[i, useidx])
	  tempmodel=djs_maskinterp(tempmodel,tempmask)
	  colmodel = djs_median( tempmodel, width=201, boundary='reflect')
	 leftmodel=djs_median(constarr+slopearr*rawcol,width=201,boundary='reflect')
         fixedimage[i, useidx] = flatimage[i,useidx]*(leftmodel/colmodel)
  endfor


return, fixedimage
end














