


;path = '/Users/joe/isaac_test_data/'
path = '/Users/joe/gabor_isaac/'
filenames = findfile(path + '/IS*.fits', count = nfile)
FOR ii = 0, nfile-1L DO BEGIN
   hdr = xheadfits(filenames[ii])
   exptime = sxpar(hdr, 'EXPTIME')
   object = sxpar(hdr, 'OBJECT')
   catg = esopar(hdr, 'HIERARCH ESO DPR CATG')
   type = esopar(hdr, 'HIERARCH ESO DPR TYPE')
   tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
   grat = esopar(hdr, 'HIERARCH ESO INS GRAT NAME')
   wcen = esopar(hdr, 'HIERARCH ESO INS GRAT WLEN')
   ob_start = esopar(hdr, 'HIERARCH ESO TPL START')
   ;;IF strmatch(tech, '*SPECTRUM*') AND abs(1.6320-wcen) LE 0.01 THEN
   ;;$
                                ;;print, catg
   forprint, fileandpath(filenames[ii]) $
             , strcompress(object, /rem), type, tech, exptime $
             , grat, wcen, ob_start, textout = 1 $
             , format = '(a40,3x,a15,a10,a10,a10,a10,f6.4,3x,a)'
ENDFOR

END
