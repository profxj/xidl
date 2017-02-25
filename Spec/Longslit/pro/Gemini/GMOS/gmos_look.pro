fils = findfile('*.par', count = nfil)
IF nfil EQ 0 THEN message, 'Cannot find plan*.par files' $
ELSE IF nfil EQ 1 THEN planfile = fils $
ELSE planfile =  x_guilist(fils)

planstr = yanny_readone(planfile, hdr = planhdr, /anonymous)

iarc = where(planstr.flavor EQ 'arc', narc)
;;wavefile = 'wave-' + planstr[iarc[0]].filename
islit = where(planstr.flavor EQ 'domeflat', nflat)
slitfile = 'slits-' + planstr[islit[0]].filename
isci = where(planstr.flavor EQ 'science', nsci)

allsci = 'Science/sci-' + planstr[isci].filename
fils = file_test(allsci)
nsci_found = total(fils)

IF nsci_found EQ 0 THEN message, 'Could not find science files' $
ELSE IF nsci_found EQ 1 THEN scifile = allsci[0] $
ELSE scifile = x_guilist(allsci)

slitmask = mrdfits(slitfile, 0)
tset_slits = mrdfits(slitfile, 1)
ystruct = mrdfits(slitfile, 2)

nx = tset_slits[0].DIMS[0]
ny = tset_slits[0].DIMS[1]
nxby3 = nx/3
idim = size(tset_slits[0].COEFF, /dim)
nslit = idim[1]

img_minsky  = mrdfits(scifile, 0, hdr)
ivar        = mrdfits(scifile, 1)
skysub_mask = mrdfits(scifile, 2)
waveimg     = mrdfits(scifile,3)
objstruct = mrdfits(scifile, 4)
;;waveimg  = mrdfits(wavefile, 0)
;; adjust to nod and shuffle region
slitmask = slitmask[nxby3:2*nxby3-1L, *]
;;waveimg  = waveimg[nxby3:2*nxby3-1L, *]

xatv, img_minsky*(skysub_mask)*(ivar GT 0.0) $
      , min = -200.0, max = 200.0, wv = waveimg $
      , sig = slitmask       
IF KEYWORD_SET(PLOTSLIT) THEN BEGIN
    PLOTSYM, 0, 0.3, THICK = 0.5, /FILL
    xatvplot, objstruct.xpos, objstruct.ypos, psym = 8
;xatvplot, left, yrow, psym = 8, color = 1
;xatvplot, right, yrow, psym = 8, color = 2
    objid = strtrim(objstruct.slit, 2) + '-' +  $
      strcompress(objstruct.objno, /rem) + ' '
    nobj = n_elements(objstruct)
    FOR j = 0L, nobj-1L DO BEGIN
        xpix = djs_median(objstruct[j].XPOS) - 5.0
        ypix = 1400
        xatvxyouts, xpix, ypix, strcompress(objid[j], /rem) $
                    , color = 'green', charsize = 2
    ENDFOR
ENDIF

;FOR islit = 0L, nslit-1L DO BEGIN
;   iplot = where(yrow[*, islit] GE ystruct[islit].YSLIT_MIN AND $
;                 yrow[*, islit] LE ystruct[islit].YSLIT_MAX, nip)
;   xatvplot, reform(left[iplot, islit], nip), reform(yrow[iplot, islit], nip)$
;             , psym = 8, color = 2
;   xatvplot, right[iplot, islit], yrow[iplot, islit], psym = 8, color = 1
;ENDFOR


END
