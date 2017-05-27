;+
; NAME:
;    fink_merge
;
; PURPOSE:
;    merges one chip's worth of data into a figure with nice ATV
;    display. returned results can be saved as fits file.
;
; CALLING SEQUENCE:
;    flux=fink_merge(chipno, [yobj, mask])
; 
; INPUTS:
;    chipno  -- which chip (1-8)
;    
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;   flux -- reduced, sky-subtracted data, painted onto the background
;           of the 2k*4k CCD from which it was extracted.  
;   yobj -- position denoting region where source is expected
;   mask -- regions hit by CR during exposure
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   works in current directory.  This is a pixel mapping, not
;   rectified by lambda.  The auxiliary routine
;   write_merged(filename,chipno) will call fink_merge, and the write
;   the output to a plain fits file, followed by a BINTABLE extension
;   providing yobj.
;
; REVISION HISTORY:
;  DPF, Jul02
;  MD, 07Jul02
;  DPF, 13jul02  reverse order of slits so they look like the chip!
;                  simplified code, using new deimos_tables
;----------------------------------------------------------------------
function fink_merge, chipno;, yobj, mask
; inputs: chipno
; outputs: flux, mask, yobj


;  delvarx, yobj
;  planfile = findfile('*.plan', count=nfile)
;  if nfile eq 0 then print, 'no plan files found!'
;  if nfile gt 1 then print, 'too many plan files found! Using the 1st'
;  planfile = planfile[0]
;  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
;             flatnames, arcnames, sciencenames

  objectcat = deimos_tables('*bintabs*', bluslits=slitcoords, slitobjmap=slitobjmap) 

  scale_pix_per_asec = 8.52

  flux = fltarr(4096, 2048)-10
  mask = intarr(4096, 2048)

  flist = findfile('slit*.fits*', count=nfile)
  slitind = lindgen(nfile)
  for i=0, nfile-1 do begin 
     h = headfits(flist[i], ext=1)
     print, i
     if (sxpar(h, 'chipno') EQ chipno) then begin 

        a = mrdfits(flist[i], 1, h, /silent)

        y0 = sxpar(h, 'slitx0')
        y1 = sxpar(h, 'slitx1')
        slitno = sxpar(h, 'slitno')

        slitind = where(slitcoords.slitno eq slitno, ct)
        if ct eq 0 then message, 'No entry in slitcoords for this slit'

        ind = where(slitobjmap.dslitid eq slitcoords[slitind].dslitid, ct)
        if ct eq 0 then stop
        nrow = (size(a.flux, /dim))[1]

        skyslit = slitcoords[slitind].slitwidth/ $
          median(slitcoords.slitwidth) LT 0.9

        if skyslit eq 0 then begin 
           desilen = (slitobjmap[ind].topdist+slitobjmap[ind].botdist)*scale_pix_per_asec
           edgeloss = (desilen-nrow)/2
           yobj0 = slitobjmap[ind].botdist * scale_pix_per_asec + y0 - edgeloss
           yobj = keyword_set(yobj) ? [yobj, yobj0] : [yobj0]
           slits = slitno+yobj0*0
           slitarr = keyword_set(slitarr) ? [slitarr, slits] : [slits]
        endif

        tags = tag_names(a)
        if total(tags eq 'BITMASK') then bitmask = a.bitmask AND 127 else $
              bitmask=a.mask
        print, 'adding  ', flist[i]
        flux[*, round(y0)-nrow:round(y0)-1] = flux[*, round(y0)-nrow:round(y0)-1]+a.flux+10
        mask[*, round(y0)-nrow:round(y0)-1] = mask[*, round(y0)-nrow:round(y0)-1]+bitmask

     endif 
  endfor 

  wid = 12
  atv, flux,  max=30, min=-30
  for i=0, n_elements(yobj)-1 do atvplot, [0, 4096], [1, 1]*yobj[i]+wid/2, color=2
  for i=0, n_elements(yobj)-1 do atvplot, [0, 4096], [1, 1]*yobj[i]-wid/2, color=4
  for i=0, n_elements(slitarr)-1 do atvxyouts, 0, yobj[i]-6, string(slitarr[i], format='(I3)'), charsize=2, color=2, align=1.05

  return,  flux
end

;
; calling routine to call fink_merge and then output apprpropriate
; fits file, saving the yobj as an extension

pro write_merged, filename, chipno

flux = fink_merge(chipno, yobj)

writefits, filename, flux
a = {yobj: yobj}

mwrfits, a,  filename

return
end



; routine for replotting:
flux = readfits(filename)
  atv, smooth(flux, 3), max=30, min=-5

a =  mrdfits(filename, 1)
yobj = a.yobj
wid = 12
  for i=0, n_elements(yobj)-1 do atvplot, [0, 4096], 2047-[1, 1]*yobj[i]+wid/2, color=2
  for i=0, n_elements(yobj)-1 do atvplot, [0, 4096], 2047-[1, 1]*yobj[i]-wid/2, color=4


end
