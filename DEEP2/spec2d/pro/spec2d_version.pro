;+
; NAME:
;   spec2d_version
;
; PURPOSE:
;   Return the version name for the product spec2d
;
; CALLING SEQUENCE:
;   vers = spec2d_version()
;
; INPUTS:
;   <none>
;
; OUTPUTS:
;   vers       - Version name for the product spec2d
;
; COMMENTS:
;   If this version is not tagged by CVS, then we return 'CVS: TOPLEVEL'
;   where TOPLEVEL is the last directory in the environment variable
;   $SPEC2D_DIR.  For example, if you are using a version of the code
;   in the directory '/u/schlegel/spec2d/v0_0', then this returns
;   'CVS:v0_0'.
;
; BUGS:
;
; PROCEDURES CALLED:
;
;; REVISION HISTORY:
;   01-Dec-1999  Written by D. Schlegel, Princeton.
;   22-Feb-2002  Hacked by D. Finkbeiner
;-
;------------------------------------------------------------------------------

function spec2d_version

   ; The following expression in dollar signs is expanded by CVS
   ; and replaced by the tag name for this version.
   name = '$Name: v1_1_4 $'

   words = str_sep(strcompress(name), ' ')

   if (words[0] EQ '$Name:' AND N_elements(words) EQ 3) then begin
      vers = words[1]
   endif else begin
      dirname = getenv('DEEP_DIR')
      if (dirname NE '') then begin
         words = str_sep(dirname,'/')
         nword = N_elements(words)
         vers = 'CVS:' + words[nword-1]
      endif else begin
         vers = 'CVS:Unknown'
      endelse
   endelse

   return, vers
end
;------------------------------------------------------------------------------
