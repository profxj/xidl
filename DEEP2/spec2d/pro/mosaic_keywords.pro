;+
; NAME:
;   mosaic_keywords
;
; PURPOSE:
;   Given a FITS header for a given chip, add the proper DATASEC/DETSEC and WCS
;   keywords to allow DEIMOS mosaicing
;
; CALLING SEQUENCE:
;   mosaic_keywords,header,chipnum
; 
; INPUTS:
;   header  - FITS file header
;   chipnum - chip number (so proper keywords are added)
;
; OUTPUTS:
;    header is modified
;
;
; MODIFICATION HISTORY:
;-

pro mosaic_keywords, header, chipnum


    case chipnum of
      1: begin

            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[1:2048,1:4096]', 'NOAO mosaic detector section for ds9', before='COMMENT'
            sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
        end
      2: begin
           sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[2049:4096,1:4096]', 'NOAO mosaic detector section for ds9', before='COMMENT'
           sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  2048.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
         end     
      3: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[4097:6144,1:4096]', 'NOAO mosaic detector section for ds9', before='COMMENT'

           sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
         end
      4: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[6145:8192,1:4096]', 'NOAO mosaic detector section for ds9', before='COMMENT'
           sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  6144.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
        end

      5: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[1:2048,4097:8192]', 'NOAO mosaic detector section for ds9', before='COMMENT'
            sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  0.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME','pane    ', 'UCO/Lick mosaic PANE coordinates'
        end
      6: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[2049:4096,4097:8192]', 'NOAO mosaic detector section for ds9', before='COMMENT'

            sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  2048.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
        end
      7: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[4097:6144,4097:8192]', 'NOAO mosaic detector section for ds9', before='COMMENT'
            sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
        end
      8: begin
            sxaddpar, header, 'DATASEC', '[1:2048,1:4096]', 'NOAO/IRAF data section', before='COMMENT'
            sxaddpar, header, 'DETSIZE', '[1:8192,1:8192]', 'NOAO mosaic detector size for ds9', before='COMMENT'
            sxaddpar, header, 'DETSEC', '[6145:8192,4097:8192]', 'NOAO mosaic detector section for ds9', before='COMMENT'
            sxaddpar, header,'CRPIX1',                  1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRPIX2',                   1., 'reference pixel along FITS axis j'
            sxaddpar, header,'CRVAL1',                  6144.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CRVAL2',                  4096.5, 'coordinate value for WCS axis i at refpix'
            sxaddpar, header,'CD1_1',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD1_2',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar, header,'CD2_1',                   0., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CD2_2',                   1., 'CTM i_j from pixel to WCS'
            sxaddpar,header, 'CTYPE1', 'PANE_X  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CTYPE2', 'PANE_Y  ', 'coordinate/projection type for WCS axis i'
            sxaddpar,header, 'CUNIT1', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'CUNIT2', 'pixel   ', 'physical unit for WCS axis i'
            sxaddpar,header, 'WCSNAME', 'pane    ', 'UCO/Lick mosaic PANE coordinates'
        end
    endcase

  return 
end







