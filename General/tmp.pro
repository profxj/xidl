
From swara@ociw.edu Tue Jan  8 16:52:35 2002
Date: Tue, 8 Jan 2002 16:51:45 -0800 (PST)
From: Swara Ravindranath <swara@ociw.edu>
To: xavier@ociw.edu
Subject: gal.pro

pro gal
;
   ps_open,file='swg.ps',/color,/portrait
;
 image=readfits('b.fits')
;
; im1=image(355:670,360:675)
 im1=rotate(image,90)
 im2=im1+abs(min(im1))+1
 im1=alog10(im2)
 im2=bytscl(im1,MIN=0,MAX=2.9)
 im1=255-im2
; px=!x.window*!d.x_vsize
; py=!y.window*!d.y_vsize
 sz=size(im1)
; erase
 device,/inches,xsize=6.0,scale_factor=1.0
 device,/inches,ysize=6.0,scale_factor=1.0
 xdum=findgen(10)*0.1
 ydum=findgen(10)*0.1
 dec = strarr(6)
 RA = strarr(6)
 dec = ["+29!E!9%!N!317!9'!37!9''!3 ","!317!9'!346!9''!3",$
 "!318!9'!325!9''!3","!319!9'!34!9''!3",$
 "!319!9'!343!9''!3","+29!E!9%!N!320!9'!322!9''!3"]
 RA =["!32!E!17h!N !334!E!17m!N !39!E!17s!N !3",$
 "!32!E!17h!N !334!E!17m!N !312!E!17s!N !3",$
 "!32!E!17h!N !334!E!17m!N !315!E!17s!N !3",$
 "!32!E!17h!N !334!E!17m!N !318!E!17s!N !3",$
 "!32!E!17h!N !334!E!17m!N !321!E!17s!N !3",$
 "!32!E!17h!N !334!E!17m!N !324!E!17s!N !3"]
;px=!x.window*!d.x_vsize
;py=!y.window*!d.y_vsize
 plot,xdum,ydum,TICKLEN=-0.02,xticks=5,xtickname=dec,yticks=5,ytickname=ra,$
 TITLE = '!N!17!N !17',XTITLE='!N!7d!N !5(2000)',$
 YTITLE='!N!7a!N !5(2000)', /nodata
; erase
 tvscl,im1,!x.window(0.21),!y.window(0.21),$
 xsize =!x.window(1.0)-!x.window(0.0),$
 ysize =!y.window(1.0)-!y.window(0.0), /NORM
;
;
 ps_close,/noprint,/noid
return
end
