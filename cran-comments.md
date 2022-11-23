
## Resubmission
This is a resubmission. After the manual check by the CRAN's team, I:

* Added a reference in the DESCRIPTION according to the suggested syntax.
* Added the field \value to the Rd files when it was missing
* Wrote TRUE and FALSE instead of T and F
* Set the print call in trophicSDM.R within if (verbose) as suggested
* Reset to user's par() after changing it
* Ensured that I do not use more than two cores in examples, vignettes, etc.
* Precompiled vignette to reduce check time
* Replaced \dontrun with \donttest in documentation
* Remove "this package" and fixed a typo in package description

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION:
  Poggiato (4:558)
  Thuiller (4:594)
  Trophic (4:101)
  oletti (4:574)
  trophic (4:31, 4:158, 4:364)
  
  However, the word 'Trophic' (i.e., relating to feeding and nutrition) is not misspelled and is a common word in ecology.  Poggiato, Thuiller and Andréoletti are the names of the authors of the methodological paper describing the method implemented in the package.
  
Found the following (possibly) invalid DOIs:    DOI: 10.22541/au.166853394.45823739/v1      From: DESCRIPTION      Status: Not Found      Message: 404
  
  The DOI is valid.