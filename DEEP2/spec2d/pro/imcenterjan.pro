pro imcenterjan, img, x, y, xcen, ycen, error,  $
              ;badpix = badpix,  $
              bigbox = bigbox, cbox = cbox, iter = iter,  $
              maxpix = maxpix, nomaxpix = nomaxpix,  $
              maxshift = MAXSHIFT,  $
              silent = silent

;+
; program to calculate the center of mass of an image around
; the point (x,y), return the answer in (xcen,ycen).
;
; can roughly correct for the effect of bad pixels
;
; ALGORITHM:
;   1. first finds max pixel value in
;	   a 'bigbox' box around the cursor
;   2. then calculates centroid around the object 
;   3. iterates, recalculating the center of mass 
;      around centroid until the shifts become smaller 
;      than MINSHIFT (0.3 pixels) or until it does 
;      this MAXITER number of times 
;      (** iteration is currently disabled/untested **)
;
; INPUTS:
;   img		image
;   x,y		initial guess for center
;
; OUTPUTS:
;   xcen,ycen	derived centers
;   error	error code for overshifting (0=no,1=yes)
;
; CAVEATS:
;   - works on an integer pixel basis
;   - maxpix algorithm will fail if a hot pixel is nearby
;   - iteration scheme actually seems to lead to larger 
;     differences than not iterating (by about 0.01 pixel) when
;     compared to the results of the Goddard package's CNTRD
;
;
; KEYWORDS:
;   badpix	bad pixel mask; any bad pixels inside centering
;		boxes will be corrected with 'fixpix' routine
; 		before centering REMOVED
;   bigbox	size of box for max pixel finding (default=11)
;   cbox	size of box for centroiding (default=11)
;   iter	max number of iterations (default=1)
;   /maxpix	return location of max pixel value in the box only,
;		instead of then also finding centroid
;   /nomaxpix	don't find max pixel before finding centroid
;      REMOVED
;   maxshift    if shift btwn input & output centroid is larger 
;                 then this, will give an error message and set error=1
;   /silent	don't print out error message
;
;
; Written by : M. Liu  12/12/94
; 06/02/95 (MCL): corrected bug (x and y reversed), 
;	          implemented two step algorithm & maxpix feature
; 07/14/95 (MCL): added iteration scheme 
; 07/26/95 (MCL): added bad pixel fixing option
; 04/16/96 (MCL): if /maxpix and find >1 max pixel, will find centroid
; 12/01/96 (MCL): fixed error in checking if cbox and bigbox are odd 
; 05/26/99 (MCL): big fix - round CBOX and BIGBOX to integers
;                   may screw up the centroiding via BIGBOX step
;                   (b/c may lead to 1 pixel shift?)
;                 changed default 'bigbox' from 11 pix to 1.5*CBOX
; 12/04/99 (MCL): made MAXSHIFT a keyword parm
; 11/27/01 (JAN): major changes
;-

; iteration controls
MAXITER = 1      ; disabled
MINSHIFT = 0.3

; max possible x or y direction shift
if not(keyword_set(MAXSHIFT)) then MAXSHIFT = 30

if keyword_set(cbox) eq 0 then cbox=21
if keyword_set(bigbox) eq 0 then bigbox=1.5*cbox
if keyword_set(iter) eq 0 then iter=1
if n_params() lt 5 then begin
	print,'pro imcenterf,img,x,y,xcen,ycen,(error),'
	print,'              [badpix=],[bigbox='+strc(bigbox)+'],[cbox='+strc(cbox)+'],[iter='+strc(iter)+'],[maxpix],[nomaxpix],[silent]'
        return
endif

error = 0

sz = size(img)
if keyword_set(badpix) eq 0 then badpix = fltarr(sz(1),sz(2))+1.0


; box sizes must be integers
cbox = round(cbox)
bigbox = round(bigbox)


; box size must be odd
if (cbox mod 2 eq 0) then cbox = cbox+1
dc = (cbox-1)/2
if (bigbox mod 2 eq 0) then bigbox = bigbox+1
;if ((bigbox+1)/2 ne bigbox/2.0) then bigbox = bigbox+1
db = (bigbox-1)/2


; need to start with integers
xx = round(x)
yy = round(y)


;
; first find max pixel in box around the cursor
;

;
; then find centroid if desired
;
;if not(keyword_set(maxpix)) then begin
;if not(keyword_set(maxpix)) or (n_elements(xcen) ne 1) then begin
    
    if  (n_elements(xcen) gt 1) then begin
        message, 'multiple max pixels, ignoring /maxpix and ' + $
          'calculating centroid', /info 
        xx = round(total(xcen)/n_elements(xcen)) 
        yy = round(total(ycen)/n_elements(ycen)) 
    endif
    
    done = 0
    niter = 1
;    repeat begin
        
;	cut out relevant portion
        sz = size(img)
        x0 = (xx-db) > 0		; need the ()'s
        x1 = (xx+db) < (sz(1)-1)
        y0 = (yy-db) > 0
        y1 = (yy+db) < (sz(2)-1)
        xs = x1 - x0 + 1
        ys = y1 - y0 + 1
        cut = float(img(x0:x1, y0:y1))
        cutmax = max(cut)
        w=where(cut EQ cutmax)
        cutsize = size(cut)
        my = (floor(w/cutsize[1]))[0]
        mx = (w - my*cutsize[1])[0]


        xx = mx + x0
        yy = my + y0 
        xcen = xx
        ycen = yy

; then find centroid 
if  (n_elements(xcen) gt 1) then begin
    xx = round(total(xcen)/n_elements(xcen)) 
    yy = round(total(ycen)/n_elements(ycen)) 
endif

done = 0
niter = 1
    
;	cut out relevant portion
sz = size(img)
x0 = round((xx-dc) > 0)		; need the ()'s
x1 = round((xx+dc) < (sz[1]-1))
y0 = round((yy-dc) > 0)
y1 = round((yy+dc) < (sz[2]-1))
xs = x1 - x0 + 1
ys = y1 - y0 + 1
cut = float(img[x0:x1, y0:y1])

indx = lindgen(xs, ys)
xcoord = indx mod xs
ycoord =  indx / xs
                     
           ; find x position of center of mass
;cenxx = fltarr(xs, ys, /nozero)
;for i = 0L, (xs-1) do $         ; column loop
;  cenxx[i, *] = cut[i, *] * i
xcen = total(cut*xcoord) / total(cut) + x0

                                ; find y position of center of mass
;cenyy = fltarr(xs, ys, /nozero)
;for i = 0L, (ys-1) do $         ; row loop
;  cenyy[*, i] = cut[*, i] * i
ycen = total(cut*ycoord) / total(cut) + y0




end



