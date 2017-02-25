;+
; NAME: readmhdufits.pro
;
; PURPOSE: 
;       Given a multi-HDU FITS file, this routine will assemble the
;       various components into a single array based on the header
;       keywords describing the layout.
;
; CALLING SEQUENCE: 
;       array = readmhdufits(filename)
;
; INPUTS:
;	filename = name of the multi-HDU FITS file to read
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       NOTRIM   - if set, do not trim pre- and post-scan columns
;       NOBIAS   - if set, do not perform overscan subtraction
;       LINEBIAS - if set, remove overscan line by line [default = use
;                  scalar value]
;       GAINDATA - if set to a structure with gain for each amp, make
;                  gain correction
;       VERBOSE  - if set, give feedback
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;       HEADER  - retrieve the header from the primary HDU and return
;                 it as a string array
;
; OUTPUTS:
;       This function will return a 2-D floating-point array
;       representing the assembled image.
;
; REQUIREMENTS:
;       - Requires the IDL Astronomy Users Library routines 
;       (http://idlastro.gsfc.nasa.gov/)
;
; EXAMPLES:
;       1) Read in a FITS file, with trimming and bias removal,
;       returning the data as "array" and the header as "header":
;               array = readmhdufits( 'lred0001.fits', header=header) 
;
;       2) Read in a FITS file, without trimming or bias removal,
;       returning the data as "array" and the header as "header":
;               array = readmhdufits( 'lred0001.fits', /notrim,
;               /nobias, header=header) 
;
;
;       3) Perform gain correction and use line-by-line bias
;       determination:
;               gain = [1.,2.,3.,4] ; for vidInp1,2,3,4
;               gaindata = build_gaindata(gain)
;               array = readmhdufits( 'lred0001.fits', /linebias, 
;                       gaindata=gaindata)
;	
; PROCEDURE:
;	- Uses the fits keyword in the header extentions 
;        DETSEC  = '[4096:3073,1:4096]' / NOAO mosaic detector section for ds9
;	to piece the mosiac together.
;
; AUTHOR:
;       Marc Kassis, W. M. Keck Obseravtory
;
; MODIFICATION HISTORY:
;	2009-May-27	MKassis v0.0	Original version
;	2009-Jun-02	GWirth	v0.1	- added secparse routine
;					- now returns type FLOAT array
;					- added optional baseline removal
;	2009-Jun-16	GDW	v0.2	- fix problem with NOTRIM mode
;                                       - use MRDFITS instead of FITSREAD
;       2009-Jun-19     GDW     v0.3    - adapt for binned data
;                                       - write lris_read_amp function
;       2009-Jun-19     GDW     v0.5    - fix bug with PRELINE
;       2009-Jun-24     JMS     v0.6    - fix bug with BZERO
;       2009-Jul-02     GDW     v0.7    - add LINEBIAS option
;                                       - add GAINDATA option
;       2009-Jul-06     JMS     v0.8    - multiply gain *before* bias
;                                           subtraction
;	2009-Aug-21     JXP/MK  v0.9    - merged mods from 0.8 with
;					  mods to correctly treat the
;					  pre and post cols when
;					  binned. JXP originally
;					  modified v0.7. 
;	2010-May-20	MK	v1.0	- found 1 pix X offset in all
;					  amplifiers. Corrected. DETSEC
;					  & DATASEC are indexed the same.
;					  DETSEC is 1,1 & DATASEC is 0,0
;-
;------------------------------------------------------------------------

;-----------------------------------------------------------------------
function stringify, value, format
;-----------------------------------------------------------------------
;+
; NAME:
;	STRINGIFY
;
; PURPOSE:
;	This function converts a real value into a string based on the 
;	specified format.
;
; CATEGORY:
;	Formatting.
;
; CALLING SEQUENCE:
;	Result = STRINGIFY( Value, Format)
;
; INPUTS:
;	Value	Value which is to be converted.
;
; OPTIONAL INPUTS:
;	Format	Optional format statement
;	
; OUTPUTS:
;	This function returns a string version of the input array with 
;	leading and trailing spaces removed.
;
; RESTRICTIONS:
;	None
;
; EXAMPLE:
;	1) format a number into a string:
;		label = "X=" + stringify(x)
;
;	1) format a number into a string, with formatting:
;		label = "X=" + stringify(x,'(f15.3)')
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Dec-14	GDW	Original version
;-
;-----------------------------------------------------------------------

; verify input...
np = n_params()
if np lt 1 or np gt 2 then message, 'wrong number of parameters'

; format the value into an appropriate string...
if np gt 1 then begin
    text = string( format=format, value)
endif else begin
    text = string( value)
endelse

; remove leading and trailing spaces...
return, strtrim( text, 2)
end


;------------------------------------------------------------------------
pro secparse, section, x1, x2, y1, y2
;------------------------------------------------------------------------
; Parse a section string of the form [x1:x2,y1:y2], returning
; four values
;------------------------------------------------------------------------
buf = (stregex( section, '\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)\]', $
              /subexp, /extract))[1:4] - 1
x1 = buf[0]
x2 = buf[1]
y1 = buf[2]
y2 = buf[3]
end

;-----------------------------------------------------------------------
function build_gaindata, gain
;-----------------------------------------------------------------------
; convert gain into a structure...
;-----------------------------------------------------------------------
n = n_elements(gain)
vidinp = ['VidInp1','VidInp2','VidInp3','VidInp4']

foo = {GAINDATA, vidinp:'',gain:0.}
gaindata = replicate({GAINDATA}, n)
for i=0,n-1 do begin
    gaindata[i] = {GAINDATA, vidinp:vidinp[i], gain:gain[i]}
endfor 

return, gaindata
end

;-----------------------------------------------------------------------
function gainvalue, gaindata, header
;-----------------------------------------------------------------------
; Return gain value for this amp
;-----------------------------------------------------------------------

;; initial value...
gain = !values.f_nan

;; check what amp this is...
vidinp = sxpar( header, 'extname', count=count)
vidinp = strtrim( vidinp, 2)

;; trap for keyword notfound...
if count ne 1 then begin
    message, 'ERROR: expected 1 EXTNAME keyword in header; got '+stringify(count)
endif 

;; search for matching name...
n = n_elements(gaindata)
for i=0,n-1 do begin
    if gaindata[i].vidinp eq vidinp then begin
        gain = gaindata[i].gain
    endif
endfor 

;; verify match...
if ~ finite(gain) then message, 'ERROR: no gain value given for amp '+vidinp

return, gain
end 


;------------------------------------------------------------------------
function lris_read_amp, filename, ext, $
  linebias=linebias, nobias=nobias, $
  predata=predata, postdata=postdata, header=header, $
  x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
;------------------------------------------------------------------------
; Read one amp from LRIS mHDU image
;------------------------------------------------------------------------

;; Get the pre and post pix values
;; for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
header   = headfits(filename)
precol   = sxpar(header, 'precol')
postpix  = sxpar(header, 'postpix')

;; Deal with binning
BINNING = sxpar(header, 'BINNING')
buf = (stregex( binning, '([0-9]+),([0-9]+)', /subexp, /extract))[1:2]
xbin = buf[0]
ybin = buf[1]
precol = precol/xbin
postpix = postpix/xbin

;; get entire extension...
temp = mrdfits( filename, ext, header, /silent, /unsigned)
tsize = size(temp)

nxt = tsize[1]

;; parse the DETSEC keyword to determine the size of the array.
detsec = sxpar(header, 'DETSEC')
secparse, detsec, x1, x2, y1, y2

;; parse the DATASEC keyword to determine the size of the science region
datasec = sxpar(header, 'DATASEC')
secparse, datasec, xdata1, xdata2, ydata1, ydata2

;; grab the components...
predata  = temp[0:precol-1,*]
; datasec appears to have the x value for the keywords that are zero
; based. This is only true in the image header extensions
; not true in the main header.
;data     = temp[xdata1-1:xdata2-1,*]
data     = temp[xdata1:xdata2,*]
postdata = temp[nxt-postpix:nxt-1,*]

;; flip in X as needed...
if(x1 gt x2) then begin 
    xt=x2
    x2=x1
    x1=xt
    data = reverse(temporary(data),1)
end

;; flip in Y as needed...
if(y1 gt y2) then begin 
    yt=y2
    y2=y1
    y1=yt
    data  = reverse(temporary(data), 2)
    predata  = reverse(temporary(predata), 2)
    postdata = reverse(temporary(postdata),2)
end

;; correct gain if requested...
if keyword_set(GAINDATA) then begin
    gain = gainvalue( gaindata, header)
    data = FLOAT(temporary(data)) * gain
    predata = FLOAT(temporary(predata)) * gain
    postdata = FLOAT(temporary(postdata)) * gain
endif

;; optional bias subtraction...
if ~ keyword_set(NOBIAS) then begin
    if keyword_set( LINEBIAS) then begin
        ;; compute a bias for each line...
        bias = median( postdata, dim=1)

        ;; subtract for data...
        buf = size(data)
        nx = buf[1]
        ny = buf[2]
        data2 = fltarr(nx,ny)
        for i=0,nx-1 do begin
            data2[i,*] = float(data[i,*]) - bias
        endfor 
        data = data2
    endif else begin
        ;; compute a scalar bias....
        bias = median( postdata)
        data -= bias
    endelse
endif

return, data

end

;------------------------------------------------------------------------
function x_readmhdufits, filename, HEADER=header, $
  NOTRIM=notrim, NOBIAS=nobias, LINEBIAS=linebias, $
  GAINDATA=gaindata, VERBOSE=verbose
;------------------------------------------------------------------------

;; init vars
exten=0                         ; number of extensions
n_ext=0                        ; number of extensions in input image
errmsg=''                       ;  error message set to null
xmax = 0
ymax = 0
xmin = 10000
ymin = 10000

;; verify access to file...
if ~ file_test( filename, /read) then begin
    message, 'Cannot open requested file :' + filename
endif

;; count extensions...
fits_info, filename, /SILENT, N_ext=n_ext

;; allocate array to hold the number of columns in each extension...
xcol = intarr(n_ext)

;; Get the pre and post pix values
;; for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
header  = headfits(filename)
PRECOL = sxpar(header, 'PRECOL')
POSTPIX = sxpar(header, 'POSTPIX')
PRELINE = sxpar(header, 'PRELINE')
POSTLINE = sxpar(header, 'POSTLINE')

if keyword_set(verbose) then $
  message, 'PRECOL='     + stringify(PRECOL) $
           + ' POSTPIX='  + stringify(POSTPIX) $
           + ' PRELINE='  + stringify(PRELINE) $
           + ' POSTLINE=' + stringify(POSTLINE), /info

;; get the x and y binning factors...
BINNING = sxpar(header, 'BINNING')
buf = (stregex( binning, '([0-9]+),([0-9]+)', /subexp, /extract))[1:2]
xbin = buf[0]
ybin = buf[1]

;; First read over the header info to determine the size of the output
;; array...
for i=1,n_ext do begin
    theader = headfits(filename, EXTEN=i, ERRMSG=errmsg)
    detsec = sxpar(theader, 'DETSEC')
    if(detsec NE '0') then begin 
        
        ;;parse the DETSEC keyword to determine the size of the array.
        secparse, detsec, x1, x2, y1, y2

        ;; find the range of detector space occupied by the data 
        ;; [xmin:xmax,ymin:ymax]
        xt = x2 > x1
        xmax = xt > xmax 
        yt = y2 > y1
        ymax = yt > ymax 

        ;; find the min size of the array
        xt = x1 < x2
        xmin = xmin < xt
        yt = y1 < y2
        ymin = ymin < yt
        
        xcol[i-1] = xt

    endif 
endfor

;; determine the output array size...
nx = xmax - xmin + 1
ny = ymax - ymin + 1

;; change size for binning...
nx = nx / xbin
ny = ny / ybin

;; Update PRECOL and POSTPIX
precol = precol / xbin
postpix = postpix / xbin

;; change size for pre/postscan...
if keyword_set(NOTRIM) then begin
    nx += n_ext*(precol+postpix)
    ny += preline + postline
endif 

;; allocate output array...
array = fltarr(nx, ny)
if keyword_set(VERBOSE) then $
  message, 'Creating an array of size '+stringify(nx)+ ' by '+stringify(ny), /info

order = sort(xcol)

;; insert extensions into master image...
for i=1,n_ext do begin

    ;; grab complete extension...
    data = lris_read_amp( filename, i, linebias=linebias, nobias=nobias, $
                          predata=predata, postdata=postdata, $
                          x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)

    ;; insert components into output array...
    if keyword_set(NOTRIM) then begin

        ;; insert predata...
        buf = size(predata)
        nxpre = buf[1]
        xs = order[i-1]*PRECOL
        xe = xs + nxpre - 1
        if keyword_set(VERBOSE) then begin
            section = '['+stringify(xs)+':'+stringify(xe)+',*]'
            message, 'inserting extension '+stringify(i)+ $
                     ' predata  in '+section, /info
        endif 
        array[xs:xe,*] = predata

        ;; insert data...
        buf = size(data)
        nxdata = buf[1]
        xs = n_ext*precol + (x1-xmin)/xbin
        xe = xs + nxdata - 1
        if keyword_set(VERBOSE) then begin
            section = '['+stringify(xs)+':'+stringify(xe)+',*]'

            message, 'inserting extension '+stringify(i)+ $
                     ' data     in '+section, /info
        endif 
        array[xs:xe,*] = data

        ;; insert postdata...
        buf = size(postdata)
        nxpost = buf[1]
        xs = nx - n_ext*postpix + order[i-1]*postpix
        xe = xs + nxpost - 1
        if keyword_set(VERBOSE) then begin
            section = '['+stringify(xs)+':'+stringify(xe)+',*]'
            message, 'inserting extension '+stringify(i)+ $
                     ' postdata in '+section, /info
        endif 
        array[xs:xe,*] = postdata

    endif else begin
        buf = size(data)
        nxdata = buf[1]
        nydata = buf[2]

        xs = (x1-xmin)/xbin
        xe = xs + nxdata - 1
        ys = (y1-ymin)/ybin
        ye = ys + nydata - 1 - postline

        yin1 = PRELINE
        yin2 = nydata - POSTLINE - 1

        if keyword_set(VERBOSE) then begin
            section = '['+stringify(xs)+':'+stringify(xe)+ $
                      ','+stringify(ys)+':'+stringify(ye)+']'
            message, 'inserting extension '+stringify(i)+ $
                     ' data     in '+section, /info
        endif 
       
        array[xs:xe,ys:ye] = data[*,yin1:yin2]

    endelse
end

;; make sure BZERO is a valid integer for IRAF
OBZERO = sxpar(header, 'BZERO')
sxaddpar, header, 'O_BZERO', obzero
sxaddpar, header, 'BZERO', 32768-obzero

return, array

end
