library(tidyverse)
library(R.utils)


ms_to_sirius_output = function(folder_with_ms_files,
                               folder_with_SIRIUS,
                               output_dir_name = character()) {
  
  setwd(folder_with_ms_files)
  filenames = dir(pattern="*.ms")
  
  for (file in filenames) {
    print(file)
    setwd(folder_with_SIRIUS)
    system(paste("cd ", folder_with_SIRIUS))
    
    command = str_glue(" sirius ",
                       "-i ", str_glue('"', folder_with_ms_files, "/", file, '"', sep = ""),
                       " -o ", str_glue('"', folder_with_ms_files, "/", output_dir_name, '"', sep = ""),
                       " formula",
                       #" formula -c 10",
                       " -p orbitrap", #qtof
                       " --ppm-max 5",
                       " --ppm-max-ms2 5",
                       " -E CH -e ONP[8]B[11]Si[9]S[12]Cl[18]Se[2]Br[10]FI[6]K[1]Na[1]As[2]",
                       " -d ALL",  #ALL
                       " fingerprint", #added for sirius 5
                       " structure",
                       #" canopus",
                       " write-summaries", #added for sirius 5
                       sep =" ")
    
    
    
    tryCatch(
      
      expr = withTimeout(
        {
          javaOutput = system(command, intern = TRUE) #goes into commant prompt
        },
        timeout = 300
      )
      ,
      error = function(e) {
        return(tibble())
      }
    )
  }
}
