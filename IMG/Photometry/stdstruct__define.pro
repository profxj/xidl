pro stdstruct__define

;  This routine defines the structure for standard stars

  tmp = {stdstruct, $
         flg_anly: 0,  $    ; Flag to include in analysis (0=No)
         Name: ' ',     $    ; Maps image onto entire window  
         filter: ' ',   $    ; Filter (string)
         Mag: 0.,       $    ; Magnitude
         sig_Mag: 0.,   $    ; Sigma in Magnitude
         AM: 0.,        $    ; Air Mass
         CLR: 0.,       $    ; Color Term (Value)
         sCLR: ' '      $
         }

end
  
         
