pro mike_mproc, filename, image, invvar=invvar, hdr=hdr, silent=silent

    if NOT keyword_set(silent) then silent = 0

    orig = xmrdfits(filename, 0, hdr, silent=silent)
    if n_elements(orig) EQ 1 then return
    if (size(orig))[0] NE 2 then begin
      print, 'Not a two-dimensional fits image'
      return
    endif

    ocol = (size(orig))[1]
    orow = (size(orig))[2]

    ncol  = long(ocol/ 256)*256L
    nrow  = long(orow/ 512)*512L
    bincol = 2048.0 / ncol
    binrow = 4096.0 / nrow

    if keyword_set(silent) EQ 0 then begin
      print, "Guessing ", ncol, " by", nrow, " science image", $
                 format='(a,i5,a,i6,a)'
      print, "  and so ", bincol, " by", binrow, " binning", $
                 format='(a,f5.2,a,f5.2,a)'
    endif
  

;
;   Throw away overscan rows....
;

    prebias = orig[0:ncol-1,*]
    bias = orig[ncol+1:*,*]
    medbias = djs_median(bias,1)
    postbias = prebias - medbias ## replicate(1, ncol)

    overscan = postbias[*,nrow+1:*]
    medoverscan = djs_median(overscan,2)
    image = postbias[*,0:nrow-1] - medoverscan # replicate(1, nrow)

    ronoise = 4.1
    if ARG_PRESENT(invvar) then begin
      invvar = 1.0/(abs(image) + ronoise^2)
    endif
    return
end
