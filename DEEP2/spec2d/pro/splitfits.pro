pro splitfits, filename, extno

if n_elements(extno) eq 0 then extno = -1
if extno eq 0 then extno = 1

extin = extno
if extno eq -1 then extin = 7

prefixloc = strpos(filename, 'fits')
prefix = strmid(filename, 0, prefixloc)

errmsg = ''
while(extno lt 11 and extno le extin) do begin
   struct = mrdfits(filename, extno, header)
   fxhclean, header, err=errmsg
   print, errmsg
   if n_tags(struct) ne 0 then begin
      tags = tag_names(struct)
      ntags = n_elements(tags)
      for i=0, ntags-1 do begin
         if tags[i] eq 'LAMBDA' then begin
            if total(tags eq 'LAMBDA0') ne 0 then begin
               s = size(struct.lambda)
               ncol = s[2]
               outarr = struct.lambda+struct.lambda0 # (fltarr(ncol)+1.)
            endif else outarr = struct.lambda
            mkhdr, hdrbase, outarr           
            sxdelpar,hdrbase,'END'
            sxdelpar,hdrbase,''
            headtmp = [hdrbase, header]
            mwrfits, outarr, prefix+tags[i]+'.'+string(extno, form='(i1)')+'.fits', headtmp, /CREATE, /NO_comment
         endif else if tags[i] ne 'LAMBDA0' then begin
            outarr = struct.(i)
            mkhdr, hdrbase, outarr            
            sxdelpar,hdrbase,'END'
            sxdelpar,hdrbase,''
            headmod = header
            if tags[i] ne 'FLUX' then sxdelpar, headmod, 'COMMENT'
            headtmp = [hdrbase, headmod]
            mwrfits, outarr, prefix+tags[i]+'.'+string(extno, form='(i1)')+'.fits', headtmp, /CREATE, /no_comment
            endif
;         endelse
      endfor

    endif
      extno = extno+1
endwhile



return
end
