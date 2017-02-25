
pro bspline_surface, image, image_mask, image_fit, bspline_coeffs, nxbkpt=nxbkpt, nybkpt=nybkpt, nord=nord

    if NOT keyword_set(nord) then nord = 4L
    if NOT keyword_set(nxbkpt) then nxbkpt=4L
    if NOT keyword_set(nybkpt) then nybkpt=4L

    if (size(image))[0] NE 2 then return


    sz = size(image,/dimen)

    if sz[0] mod nxbkpt NE 0 then begin
      print, 'Only works with image columns as a multiple of nxbkpt'
      print, sz[0], nxbkpt
      return
    endif
    if sz[1] mod nybkpt NE 0 then begin
      print, 'Only works with image rows as a multiple of nybkpt'
      print, sz[1], nybkpt
      return
    endif

    x = dindgen(sz[0])
    fullxbkpt = (findgen(nxbkpt-1+2*nord) - (nord-1))*(sz[0]/nxbkpt) - 1
    x_set = create_bsplineset(fullxbkpt, nord)
    x_action = bspline_action(x, x_set, lower=x_lower, upper=x_upper)
    nxc = n_elements(x_set.coeff)



    y = dindgen(sz[1])
    fullybkpt = (findgen(nybkpt-1+2*nord) - (nord-1))*(sz[1]/nybkpt) - 1
    y_set = create_bsplineset(fullybkpt, nord)
    y_action = bspline_action(y, y_set, lower=y_lower, upper=y_upper)
    nyc = n_elements(y_set.coeff)

    bw = nxc*nord
    full_beta = dblarr(nxc*nyc + bw)
    full_alpha = dblarr(nxc*nord, nxc*nyc+bw)
    test_alpha = dblarr(nxc*nyc, nxc*nyc)

    x_sub = (x_upper[0]-x_lower[0]+1)
    y_sub = (x_upper[0]-x_lower[0]+1)
    nd_full = x_sub * y_sub
    full_action = (transpose(x_action[x_lower[0]:x_upper[0], *]))[*] # (y_action[y_lower[0]:y_upper[0], *])[*]
    full_action = reform(transpose(reform(full_action, nord, nd_full, nord), [0,2,1]),nord*nord,nd_full)

    full_corr = matrix_multiply(full_action, full_action, /btranspose)

    fbw = nord*nord
    bi = lindgen(fbw)
    bo = lindgen(fbw)
    for i=1L, fbw-1 do bi = [bi, lindgen(fbw-i)+(fbw+1)*i]
    for i=1L, fbw-1 do bo = [bo, lindgen(fbw-i)+fbw*i]

    coeffs = (lindgen(nord)#replicate(1,nord) + lindgen(nord)##replicate(nxc,nord))[*]
   
    fc = coeffs # coeffs * 0 
    for iy=0,nord-1 do for ix=0,nord-1 do fc[*,ix+iy*nord] = coeffs + (ix + iy*nxc)*(nxc*nord - 1)

stop
    for iy=0,nyc-nord do begin & $
      for ix=0,nxc-nord  do begin & $
        itop = ix + iy*nxc  & $
        full_alpha[fc[bi]+itop*bw] = full_alpha[fc[bi]+itop*bw] + full_corr[bi] & $
        ix_spot = (lindgen(nord)##replicate(nxc,nord) + lindgen(nord) # replicate(1, nord)+itop)[*] # replicate(1,nord*nord)   & $
        iy_spot = transpose(ix_spot) & $
        test_alpha[ix_spot,iy_spot] = test_alpha[ix_spot, iy_spot] + full_corr & $
      endfor & $
    endfor
   
    for iy=0,nyc-nord do begin & $
      for ix=0,nxc-nord do begin & $
        itop = ix + iy*nxc  & $
        full_beta[itop+coeffs] = full_beta[itop+coeffs] + $
             (sqrt(image[x_lower[ix]:x_upper[ix], y_lower[iy]:y_upper[iy]]))[*] ## full_action  & $
      endfor & $
    endfor

    for iiter=1,maxiter do begin

 
      temp_alpha = full_alpha
      temp_beta = full_beta

      bad_pts = where(image_mask LE 0, nbad)
      if nbad GT 0 then begin
        badx = bad_pts mod sz[0]
        bady = bad_pts / sz[0]
        bad_subx = badx mod x_sub
        bad_suby = bady mod y_sub

        for i=0, nbad -1 do begin

          bad_action = full_action[*, bad_subx[i] + bad_suby[i]*x_sub]

          temp_alpha[] = temp_alpha - total(bad_action)^2


      ; First solution 
      err_band = cholesky_band(temp_alpha)
      err_solve = cholesky_solve(temp_alpha, temp_beta)


    image_fit = image * 0.0
    for iy=0,nyc-nord do begin & $
      for ix=0,nxc-nord do begin & $
        itop = ix + iy*nxc  & $
        image_fit[x_lower[ix]:x_upper[ix], y_lower[iy]:y_upper[iy]] = temp_beta[itop+coeffs] # full_action  & $
      endfor & $
    endfor


    return
end      
