pro test23, pixels=pixels, iwave=iwave, wave
if not keyword_set(pixels) then pixels= 3850L
if not keyword_set(iwave) then iwave= 3.5793d
wave= 10^(iwave + findgen(pixels)*(0.0001))-10^(iwave)+10^(iwave)
end
