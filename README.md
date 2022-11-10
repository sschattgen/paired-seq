# paired-seq

#### Notes on demultiplex_plates.R
Dependencies:
`tidyverse`
`data.table`
`future`
`stringdist`

The demultiplexing routine is currently run using multisessions in `future`.
Each session currently has a memory allotment of 4GB but this can be changed in the script.
