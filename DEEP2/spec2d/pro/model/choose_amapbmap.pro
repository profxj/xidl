pro choose_amapbmap,header,amapfile,bmapfile
;+
; NAME:
;
;  CHOOSE_AMAPBMAP
;
; PURPOSE:
;
;   Pick proper amap & bmap file for a particular mask
;
; CATEGORY:
;
;   Optical Model
;
; CALLING SEQUENCE:
;
;   chooseamapbmap, header, amapfile, bmapfile
;
; INPUTS:
;
;   header -- a DEIMOS header
;
; OUTPUTS:
;
;   amapfile -- filename (with path) for amap
;   bmapfile -- filename (with path) for bmap
;
; MODIFICATION HISTORY:
;   jan 2003mar05
;-

  mjd = sxpar(header, 'MJD-OBS')
  slider=sxpar(header,'GRATEPOS')

  amapstem='amap.s'+strcompress(string(slider),/remove)+'.*.sav'
  bmapstem='bmap.s'+strcompress(string(slider),/remove)+'.*.sav'

;-------------------------------------------
; choose the right amap
    searchlist=concat_dir(getenv('CALIB_DATA'), amapstem)
    amaps=findfile(searchlist,count=namaps)
    
; exclude dec. 19 superdark - only for oct. 11 data

     if namaps eq 0 then message,'You need an amap file in CALIB_DATA!  I mean it, buster!'

     amapdates=dblarr(namaps)
     for i=0,namaps-1 do amapdates[i]=juldayfromfilename(amaps[i])
   
     minseparation=min(abs(amapdates-mjd),bestdate)
   
     amapfile=amaps[bestdate] 
    
   searchlist=concat_dir(getenv('CALIB_DATA'), bmapstem)
    bmaps=findfile(searchlist,count=nbmaps)
    
; exclude dec. 19 superdark - only for oct. 11 data
     if nbmaps eq 0 then message,'You need an bmap file in CALIB_DATA!  I mean it, buster!'

     bmapdates=dblarr(nbmaps)
     for i=0,nbmaps-1 do bmapdates[i]=juldayfromfilename(bmaps[i])
   
     minseparation=min(abs(bmapdates-mjd),bestdate)
   
     bmapfile=bmaps[bestdate] 
    

return

end
