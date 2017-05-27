pro check_wavesol,outfile

 openw,lun,outfile,/GET_LUN

calibfiles = findfile('calibSlit.*.fits', count=nfiles)

for I=0,nfiles-1 do begin
  cc = mrdfits(calibfiles[i],1,hdr,/SILENT)
  rms = sxpar(hdr,'WAVESIG')
  if rms GT 0.5 then begin
    print,calibfiles[i],rms
    printf,lun,calibfiles[i],rms
  endif
endfor
close,lun

end
