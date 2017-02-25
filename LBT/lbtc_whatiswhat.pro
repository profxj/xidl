
;procedure that print to screen basic info of fits file. Sort of log 
;reconstructure

;PATH--> path to a folder where images are 
;LOGNAME--> file containing the info  


PRO lbtc_whatiswhat, PATH=path, LOGNAME=logname



  IF ~keyword_set(PATH) THEN path='./'
  IF ~keyword_set(LOGNAME) THEN  logname='log_lbtc.info'
  



;get a list of images
spawn, string("ls ",path,"*.fits*"), namelist
num=N_ELEMENTS(namelist)

;open plan
openw, lun, logname,  /get_lun


for pos=0, num-1 do begin
   
;open headers
   header=headfits(namelist[pos],exten=0)
   
   filt=fxpar(header,"FILTER")
   obj=fxpar(header,"OBJECT")   
   expt=fxpar(header,"EXPTIME")
   
   printf, lun, namelist[pos], " ", obj, " ", filt, " ", expt
   print, namelist[pos], " ", obj, " ", filt, " ", expt
   
endfor


free_lun, lun



END
