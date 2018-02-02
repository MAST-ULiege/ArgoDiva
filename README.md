Welcome to the ArgoDiva toolbox ! 

This toolbox can be used to extract diagnostics from Argo data and prepare files for their interpretation with [DIVA](http://modb.oce.ulg.ac.be/mediawiki/index.php/DIVA).

An example of use is given in [CaseTest.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/CaseTest.R).

The first step is to call [ArgoSelect.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/ArgoSelect.R), to extract a list of variables from a collection of Argo netcdf files.

Second step is to call [ArgoExtract.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/ArgoExtract.R), to compute diagnostics out of compiled profiles. The toolbox allows the user to define new diagnostics in [ArgoVertFunctions.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/ArgoVertFunctions.R).

[ArgoDisplay.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/ArgoDisplay.R) may then be used to generate pdf report on the distribution, sptail and temporal variability of the diagnostics. 

[ArgoPrepareForDiva.R](https://github.com/MAST-ULiege/ArgoDiva/blob/master/ArgoPrepareForDiva.R) generates, on the basis of user options, a data txt file, ready to be used in a DIVA analysis.



