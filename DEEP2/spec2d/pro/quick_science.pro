;+
; NAME:
;   QUICK_SCIENCE
; PURPOSE:
;   check for calibration files and start quicklook reductions
; CALLING SEQUENCE:
;   quick_science
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;   run in data directory
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;
;----------------------------------------------------------------------




pro quick_science, level,  newfile
   
  if n_elements(level) eq 0 then level = 2
  if n_elements(newfile) eq 0 then begin
     spawn, 'rsh polo wfi', newfile
     h = headfits('/net/polo'+newfile)
  endif else h = headfits('/net/polo'+newfile)
  
  mask = sxpar(h, 'SLMSKNAM')
  lamps = sxpar(h, 'LAMPS')
  exptime = sxpar(h, 'EXPTIME')
  obstype = sxpar(h, 'OBSTYPE')
  grating = sxpar(h, 'GRATENAM')
  object = sxpar(h, 'OBJECT')
  gratepos = sxpar(h, 'GRATEPOS')
  hatchpos = sxpar(h, 'HATCHPOS')
  if gratepos eq 3 then wave = sxpar(h, 'G3TLTWAV') $
    else wave = sxpar(h, 'G4TLTWAV') 


  isscience =  strpos(lamps, 'Off') ge 0 $ 
    AND exptime gt 250. $
    AND strpos(obstype, 'Object') ge 0 $
    AND wave gt 5000 $
    AND strpos(hatchpos, 'open') ge 0

  mask = strcompress(mask, /REMOVE)

     maskprocessed = n_elements(findfile(mask+'/arc*')) ge 4 $
       AND n_elements(findfile(mask+'/cal*')) ge 10


     stringsep = strsplit(strcompress(newfile, /REMOVE), '/', /extract)
     filestem = stringsep[n_elements(stringsep)-1]

     directory = '/net/polo'
     for i=0, n_elements(stringsep)-2 do directory = directory+'/'+stringsep[i]
     cd,directory
     cd,current=cwd

;     print, 'in directory '+directory
     print,  'processing file ',newfile
;       print, 'mask: ', mask
    

    if isscience AND maskprocessed then begin
       print, 'this is a science file!'
       
       cd, mask
; write quickSlit files
;            string = 'set clobber ; echo "deimos_2dreduce,
;            file='+filestem+', quick='+string(level < 9,
;            format='(i2)')+'" | idl >> science.log'
            
       deimos_2dreduce, file=filestem, quick=level

       quick_sciqa, quicklevel=level, file=filestem
       cd,cwd

    print, 'mask: ', mask, ' completed'
    
    endif


    quick_science, level

return
end


