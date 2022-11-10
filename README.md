# paired-seq

#### Notes on demultiplex_plates.R
Dependencies:  
`tidyverse`  
`data.table`  
`future`  
`stringdist`  

The demultiplexing routine is currently run using multisessions in `future`.  
Each session currently has a memory allotment of 4GB but this can be changed in the script.


#### Running the script

##### Args
species = human or mouse  
pattern = pattern for matching to alignment file names. should match A1, A2, B1, B2 plates  
input_path = path to alignment files  
output_prefix = path with prefix ending for writing files to (e.g. ~/my_plate_run)  

##### Cmd
`Rscript demultiplex_plateseq.R species pattern input_path output_prefix`  

#### Files in `data` folder should be in same location as the script when running. Need to fix that.  
