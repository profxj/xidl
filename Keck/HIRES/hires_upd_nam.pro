pro hires_upd_nam, fil, filnm, signm, contnm
	img = readfits(fil,head)
    	fxaddpar, head, 'FILENAME', filnm
    	fxaddpar, head, 'SIGFILE ', signm
    	fxaddpar, head, 'CONTFILE', contnm
	mwrfits, img, fil, head, /create
return
end
