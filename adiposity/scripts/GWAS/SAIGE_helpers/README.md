Scripts in this folder:

**SAIGE_stepN_wrapper.R** - Script to call models for SAIGE step 1 and step 2, with the functions fitNULLGLMM and SPAGMMATtest respectively, by parsing arguments provided from the command line and concatenating file paths. 
**perform_SAIGE_stepN.sh** - Shell script to execute the wrapper scripts above with command-line arguments submitted by the submit scripts in *../* that loop over different strata and clusters.