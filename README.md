# TLM-tool
a tool for analyzing TLM assay panels and auto-calling the data set based on a master identification file.

The tool itself is running csv, pandas, and os libraries.  It will analyze a kraken exported data set (csv file) and return an analyzed set complete with Trait/seg/wildtype calls for the N9, Castle, and clubroot panels.


Version 2 notes:

TLM_prod_V2 was added to repo.  This version has the following changes to version 1:

1.  No longer use a master list to add call.  The call itself is now done via a class module (mas_class_tools) in the script itself (no need to update a master list in excel, all changes can be made in the mas_class_tools
2.  tool was updated so that it scans the kraken output file for which panel is run (as it could be any combination of 3 of N9, Castle, or Clubroot panels).  
3.  adding functionality of acyto/gt200 and LepR3 call changes to the tool.  This will work on files that don't have any the TLM panels, or any combination of the 3.  This function will similarily run only if there are acyto, gt200 or lepr3 assays present in the kraken output, and the calls get converted to what they should be in Variety.
4.  added the Glob module to scan the folder for files to run through the tool, so you need only place the file in the folder, and execute the script.  finished files will be saved into the completed folder within the root folder of this tool.


Version 3 Notes:

Native Trait Analaysis Tool amd TLM class added to repo.  This version has the following changes to version 2:
1.  Re-Wrote the tool (now called Native Trait Analysis Tool) to include more trait linked marker panels in the analysis, as well the TLM class takes the place of the mas_class_tools to complete the analysis and conversions in the file.  
2.  Removed some functions no longer needed for the final merge sequence of the files with calls 
3.  The TLM class was re-written to include more allele types as well it no longer draws from a specified set of markers.  Instead it has a dictionary passed into it of a number of assays that draws from the allele types to complete the call.  The number of assays in a panel no longer matters either as it will complete the call based on the ratios of trait/seg/null alleles.
