FUNCTION GNIRS_FILEPREFIX, filename

split1 = strsplit(filename, '/', /extract)
nsplit = n_elements(split1)
file_pref1 = strsplit(split1[nsplit-1L], '\.fits', /extract,/regex)

RETURN, file_pref1[0]
END
