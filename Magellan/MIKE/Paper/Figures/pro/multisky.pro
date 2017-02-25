pro multisky

loadct, 0

file1 = 'data/PKS2000_im1.fits'
image1 = mrdfits(file1,0)   ; this is the "Final red image to be extracted"
obj1   = mrdfits(file1,1)   ; this is the best object image
sky1   = mrdfits(file1,2)   ; this is the best sky image
residual1 = mrdfits(file1,3)  ; this is image - obj - sky

full1 = [image1, obj1, sky1, residual1]

;full1he = -1 * [hist_equal(image1), hist_equal(obj1)-30, $
;         hist_equal(sky1)-30, hist_equal(residual1)]
full1he = -1 * [hist_equal(image1), hist_equal(obj1), $
         hist_equal(sky1), hist_equal(residual1)]


file2 = 'data/PKS2000_im2.fits'
image2 = mrdfits(file2,0)   ; this is the "Final red image to be extracted"
obj2   = mrdfits(file2,1)   ; this is the best object image
sky2   = mrdfits(file2,2)   ; this is the best sky image
residual2 = mrdfits(file2,3)  ; this is image - obj - sky

full2he = -1 * [hist_equal(image2), hist_equal(obj2 <103), $
         hist_equal(sky2), hist_equal(residual2)]

test = [obj1, obj2]

window, 0 
tvscl, full1he
window, 1
tvscl, full2he

x_psopen, 'pks2000_im1.ps', /portrait
tvscl, full1he
x_psclose
x_psopen, 'pks2000_im2.ps', /portrait
tvscl, full2he
x_psclose

end
