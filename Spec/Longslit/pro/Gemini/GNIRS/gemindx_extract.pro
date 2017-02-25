;---------------------------------------------------------------------------
FUNCTION gemindx_extract, filename
split1 = strsplit(filename, 'S', /extract)
nsplit = n_elements(split1)
split2 = strsplit(split1[nsplit-1L], '.fits*', /extract)
RETURN, long(split2[0])
END

