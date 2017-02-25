; + 
; NAME:
; mage_combspec
; Version 0.1
;
; PURPOSE:
;  Combines an array of obj_structs into a mage fspec structure.  This
;  allows the combination of multiple observations and must be done
;  before  making a 1-d spectrum (even of just a single observation).  
;
; CALLING SEQUENCE:
;
;  mage_combspec,objstr,fspec
;
; INPUTS:
;   objstr  - An object structure generated from the mage_script
;             extraction routines, which has been fluxed
;
; RETURNS:
;
; OUTPUTS:
;
;   fspec   - The name of the magefspecstrct produced by combining the
;             input object structures
;
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mage_combspec,obj_strct,fspec
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;    magefspecstrct__define, x_echcombspec
;
; REVISION HISTORY:
;   16-Jun-2008 CLW

pro mage_combspec,objstr,fspec,CHK=CHK

fspec={magefspecstrct}
fspec.nexp=n_elements(objstr)/15  ;because there are 15 orders

;mage_echcombspec,objstr,fspec,[6,20],0
mage_esi_echcombspec, objstr, fspec,CHK=CHK

end
