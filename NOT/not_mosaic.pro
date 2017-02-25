FUNCTION NOT_MOSAIC, file,  ihdr = ihdr, hdr_arr = hdr_arr, mask = mask

;; 1 lower left
;; 2 lower right
;; 3 upper left
;; 4 upper right

ihdr = headfits(file)
img1 = xmrdfits(file, 1, hdr1, /fscale)
img2 = xmrdfits(file, 2, hdr2, /fscale)
img3 = xmrdfits(file, 3, hdr3, /fscale)
img4 = xmrdfits(file, 4, hdr4, /fscale)

nhdr = n_elements(hdr1)
hdr_arr = strarr(nhdr, 4)
hdr_arr[*, 0] = hdr1
hdr_arr[*, 1] = hdr2
hdr_arr[*, 2] = hdr3
hdr_arr[*, 3] = hdr4

dims = size(img1, /dim)
nx = dims[0]
ny = dims[1]
xgap2 = 28L
xgap4 = 28L
ygap3 = 29L
ygap4 = 15L
ygap2 = 13L
nmx = 2*nx + max([xgap2, xgap4])
nmy = 2*ny + max([ygap2 + ygap4, ygap3])
imag = fltarr(nmx, nmy)
mask = lonarr(nmx, nmy)
;; 1 lower left
imag[0:nx-1L, 0:ny-1L] = img1
mask[0:nx-1L, 0:ny-1L] = 1
;; 3 upper left
imag[0:nx-1L, ny+ygap3:2*ny+ygap3-1L] = img3
mask[0:nx-1L, ny+ygap3:2*ny+ygap3-1L] = 3
;; 2 lower right
imag[nx+xgap2:2*nx+xgap2-1L, ygap2:ny + ygap2 -1L] = img2
mask[nx+xgap2:2*nx+xgap2-1L, ygap2:ny + ygap2 -1L] = 2
;; 4 upper right
imag[nx+xgap4:2*nx+xgap4-1L, ny + ygap2 + ygap4:2*ny+ygap2 + ygap4-1L] = img4
mask[nx+xgap4:2*nx+xgap4-1L, ny + ygap2 + ygap4:2*ny+ygap2 + ygap4-1L] = 4

RETURN, imag
END
