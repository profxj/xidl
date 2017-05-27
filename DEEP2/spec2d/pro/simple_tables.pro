; based on code from MCC
; takes bintab file name

pro simple_tables, bfile, slitnames=slitnames, mags=mags, slitwid=slitwid, slitlen=slitlen, xmm=xmm, ymm=ymm, objnames=objnames



  header = headfits(bfile)
  
  fits_open, bfile, fcb
  extnames = fcb.extname
  fits_close, fcb
  dex = where(extnames eq 'DesiSlits', cnt)
  if cnt eq 0 then message, 'No DesiSlits Table found!'
  design = mrdfits(bfile, dex[0], /silent)
  dex = where(extnames eq 'ObjectCat', cnt)
  if cnt eq 0 then message, 'No ObjectCat Table found!'
  objcat = mrdfits(bfile, dex[0], /silent)
  dex = where(extnames eq 'SlitObjMap', cnt)
  if cnt eq 0 then message, 'No ObjectCat Table found!'
  objmap = mrdfits(bfile, dex[0], /silent)
   dex = where(extnames eq 'BluSlits', cnt)
  if cnt eq 0 then message, 'No Bluslits Table found!'
  blu_in = mrdfits(bfile, dex[0], /silent)
     temp = {slitno:0L, dslitid:0L, xmm:0.0,  ymm:0.0, $
               pa:0.0, slitlength:0.0,  slitwidth:0.0, $
              xmmb:0.0,ymmb:0.0,xmmt:0.0,ymmt:0.0}
       bluslits = replicate(temp, n_elements(blu_in) )

       bluslits.dslitid = blu_in.dslitid
       bluslits.xmm = (blu_in.slitx1 + blu_in.slitx2 + $
                       blu_in.slitx3 + blu_in.slitx4)/4.
       bluslits.ymm = (blu_in.slity1 + blu_in.slity2 + $
                       blu_in.slity3 + blu_in.slity4 )/4.
       bluslits.xmmt=(blu_in.slitx1 + blu_in.slitx4)/2.
       bluslits.ymmt=(blu_in.slity1 + blu_in.slity4)/2.
       bluslits.xmmb=(blu_in.slitx2 + blu_in.slitx3)/2.
       bluslits.ymmb=(blu_in.slity2 + blu_in.slity3)/2.


       bluslits.pa  = atan(blu_in.slity2 - blu_in.slity1, $
                           blu_in.slitx2 - blu_in.slitx1)*!radeg


       bluslits.pa =  bluslits.pa*(abs(bluslits.pa) lt 90.) + $
                      (bluslits.pa+180.)*(bluslits.pa lt -90.) + $
                      (bluslits.pa-180.)*(bluslits.pa gt  90.)
       bluslits.slitlength = abs(blu_in.slitx1 - blu_in.slitx2)
       bluslits.slitwidth  = float(abs(blu_in.slity2 - blu_in.slity3))
       bluslits = bluslits[sort(bluslits.dslitid)] ;sort on dslitid
  nslits = n_elements(design.slitname)




  for ii=0,nslits-1 do begin
      dex = where(objmap.dslitid eq design[ii].dslitid, cnt)
      bdex = where(bluslits.dslitid eq design[ii].dslitid, bcnt)
      if cnt eq 0 then message, 'No entry in SlitObjMap found!'
      objid = objmap[dex].objectid
;      if cnt gt 1 then print, 'More than 1 object in this slit!' $
      if cnt eq 1 then  begin
; TBD: check that objid wont have multiple entries. if so iterate thru
; them and take the brightest object in the slit or skip it because
; the object might be hard to locate.
          dex = where(objcat.objectid eq objid[0], cnt)
          if cnt eq 0 then message, 'No entry in ObjectCat found!'

          if n_elements(slitnames) gt 0 then $
            slitnames = [slitnames, design[ii].slitname] $
          else slitnames = design[ii].slitname

          if n_elements(mags) gt 0 then $
            mags = [mags, objcat[dex[0]].mag] $
          else mags = objcat[dex[0]].mag

          if keyword_set(objnames) then $
            objnames = [objnames, objcat[dex[0]].object] $
          else objnames = objcat[dex[0]].object

          if keyword_set(slitwid) then $
            slitwid = [slitwid, design[ii].slitwid] $
          else slitwid = design[ii].slitwid

          if n_elements(slitlen) gt 0 then $
            slitlen = [slitlen, design[ii].slitlen] $
          else slitlen = design[ii].slitlen

          if n_elements(xmm) gt 0 then $
            xmm = [xmm, bluslits[bdex[0]].xmm] $
          else xmm = bluslits[bdex[0]].xmm

          if n_elements(ymm) gt 0 then $
            ymm = [ymm, bluslits[bdex[0]].ymm] $
          else ymm = bluslits[bdex[0]].ymm


      endif
  endfor











return
end
