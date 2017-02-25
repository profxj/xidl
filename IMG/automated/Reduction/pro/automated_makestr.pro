






PRO automated_makestr,logfile


   IF NOT KEYWORD_SET(path) THEN path = './'

   strname = STRARR(4)

   FOR i=0,3 DO BEGIN

      ;logfile + '_1.str'
      k = i + 1
      k = STRTRIM(k,2)
      strname[i] = logfile + '_' + k + '.str'

   ENDFOR



   ;;read in the log file
   READCOL,logfile,name,date,void,object,filter,texp,airm,FORMAT='A,A,F,A,A,F,F',COMMENT='#'
   n = N_ELEMENTS(name)

   ;;get gain a readnoise for each image in each fits file

   gain = DBLARR(4)
   rn = DBLARR(4)

   FOR i=0,3 DO BEGIN

      header = HEADFITS(name[0],EXTEN=(i+1))
      gain[i] = FXPAR(header,'GAIN')
      rn[i]= FXPAR(header,'RDNOISE')

   ENDFOR

  ; gain1 = REPLICATE(gain[0],n)
  ; rn1 = REPLICATE(rn[0],n)
  ; gain2 = REPLICATE(gain[1],n)
  ; rn2 = REPLICATE(rn[1],n)
  ; gain3 = REPLICATE(gain[2],n)
  ; rn3 = REPLICATE(rn[2],n)
  ; gain4 = REPLICATE(gain[3],n)
  ; rn4 = REPLICATE(rn[3],n)

   gainArr = PTRARR(4)
   rnArr = PTRARR(4)
   
   FOR i=0,3 DO BEGIN

      gainArr[i] = PTR_NEW(REPLICATE(gain[i],n))
      rnArr[i] = PTR_NEW(REPLICATE(rn[i],n))

   ENDFOR


   SPLOG, "Gains are",gain[0],gain[1],gain[2],gain[3]
   SPLOG, "RN's are",rn[0],rn[1],rn[2],rn[3]

   ;;create type from names
   type = STRTRIM(STRMID(name,5,3),2)
   header = HEADFITS(name[0],EXTEN=0)
   date = FXPAR(header,'DATE-OBS')
   yy = STRMID(date,2,2)
   mm=STRMID(date,5,2) 
   dd=STRMID(date,8,2)
   date=yy + mm + dd  

   bia = STRMATCH(type,'BIA')
   fla = STRMATCH(type,'FLA')
   dar = STRMATCH(type,'DAR')
   reindx = bia + fla + dar
   indx = WHERE( reindx LT 1,nn)
   type[indx] = 'IMA'
   

   ;;Prepare structure
   strTmp = {NAME:name,OBJECT:object,FILTER:filter,EXTIME:texp,AIRMASS:airm,GAIN:*gainArr[0],RN:*rnArr[0],TYPE:type}
   str = REPLICATE(strTmp,4)
   FOR i=1,3 DO BEGIN

      k = i + 1
      kString = STRTRIM(k,2)
      str[i].GAIN = *gainArr[i]
      str[i].RN = *rnArr[i]

   ENDFOR
  ; str2 = {NAME:name,OBJECT:object,FILTER:filter,EXTIME:texp,AIRMASS:airm,GAIN:gain2,RN:rn2,TYPE:type}
  ; str3 = {NAME:name,OBJECT:object,FILTER:filter,EXTIME:texp,AIRMASS:airm,GAIN:gain3,RN:rn3,TYPE:type}
  ; str4 = {NAME:name,OBJECT:object,FILTER:filter,EXTIME:texp,AIRMASS:airm,GAIN:gain4,RN:rn4,TYPE:type}


   FOR i=0,3 DO BEGIN

      k = 1 + i
      k = STRTRIM(k,2)
      mwrfits, str[i], strname[i], /create
      splog, "Structure" + k + " written to  ", strname[i]  


   END

  ; mwrfits, str[0], strname[0], /create
  ; splog, "First structure written to  ", strname[0]  
  ; mwrfits, str[1], strname[1], /create  
  ; splog, "Second structure written to  ", strname[1]
  ; mwrfits, str[2], strname[2], /create  
  ; splog, "Third structure written to  ", strname[2]
  ; mwrfits, str[3], strname[3], /create
  ; splog, "Fourth structure written to  ", strname[3]





END
