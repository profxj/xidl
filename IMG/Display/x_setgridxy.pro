;+ 
; NAME:
; x_setgridxy
;    Version 1.1
;
; PURPOSE:
;    Sets up the tv structure for displaying images.  Mainly
;   fills up the grid and image endpoints.
;
; CALLING SEQUENCE:
;   x_setgridxy, struct, imgreg, /FILL
;
; INPUTS:
;   struct     - Structure including the tv tags
;   imgreg     - Region to display
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FILL       - Fill the screen
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_setgridxy, struct, imgreg
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Nov-2001 Written by JXP
;   06-Feb-2002 Revised to allow FILL (JXP)
;-
;------------------------------------------------------------------------------
pro x_setgridxy, struct, imgreg, FILL=fill

;  Find the tv structure

  struct_nm = tag_names(struct, /STRUCTURE_NAME)
  if struct_nm EQ 'tvstruct' then flg_tv = 1 else begin
      a = tag_names(struct)
      b = where(a EQ 'TV', count)
      if count EQ 0 then begin
          print, 'Not set up for this structure'
          return
      endif
      flg_tv = 0
  endelse
  
  grid = fltarr(2)
  xymnx = fltarr(4)
  
;  Set image size

  imgsz = lonarr(2)
  imgsz[0] = imgreg[2]-imgreg[0]+1
  imgsz[1] = imgreg[3]-imgreg[1]+1

;  Set win
  if flg_tv EQ 1 then win = struct.winsize else win = struct.tv.winsize

;
  ; NO FILL
  if not keyword_set( FILL ) then begin
      if float(win[0])/float(imgsz[0]) LT $
        float(win[1])/float(imgsz[1]) then begin
          grid[0] = win[0]
          grid[1] = round( float(imgsz[1])*float(win[0])/float(imgsz[0]))
          xymnx[0] = float(imgreg[0])
          xymnx[1] = float(imgreg[1])
          xymnx[2] = float(imgreg[2])
          xymnx[3] = xymnx[1]+float(win[1])*float(imgsz[1])/float(grid[1])
      endif else begin
          grid[0] = round( float(imgsz[0])*float(win[1])/float(imgsz[1]))
          grid[1] = win[1]
          xymnx[0] = float(imgreg[0])
          xymnx[1] = float(imgreg[1])
          xymnx[2] = xymnx[0]+float(win[0])*float(imgsz[0])/float(grid[0])
          xymnx[3] = float(imgreg[3])
      endelse
      delvarx, imgsz
  endif else begin ; FILL
      grid[0] = win[0]
      grid[1] = win[1]
      xymnx[0] = float(imgreg[0])
      xymnx[1] = float(imgreg[1])
      xymnx[2] = float(imgreg[2])
      xymnx[3] = float(imgreg[3])
  endelse
      
;  Set xymnx, grid

  if flg_tv EQ 1 then begin
      struct.gridsize = grid
      struct.xymnx = xymnx
  endif else begin
      struct.tv.gridsize = grid
      struct.tv.xymnx = xymnx
  endelse

end
