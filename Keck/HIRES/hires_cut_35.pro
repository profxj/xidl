pro hires_cut_35, fil1
	img = readfits(fil1, head)
	img[0:34,*] = -1.
	writefits, fil1, img, head
return
end
