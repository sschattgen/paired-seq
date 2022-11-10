#parsing results from plateseq mixcr
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(future)){
  install.packages("future")
  library(future)
}
if(!require(stringdist)){
  install.packages("stringdist")
  library(stringdist)
}

plan(multisession, workers = 4)
options(future.globals.maxSize= 4194304000)
# preflight and functions ====

args = commandArgs(trailingOnly=TRUE)
species = args[1]#species, human or mouse
patterns = args[2]#pattern for matching to alignment file names. should match A1, A2, B1, B2 plates
input_path = args[3]#path to alignment files
output_prefix = args[4]#path with prefix ending for writing files to 
#maxreads = args[5] #max number of read, if empty set to 10e6


# if(is.null(maxreads)){
#   print('Setting maxreads to 2 x 10e6')
#   maxreads = 2000000
# } 

nucdict<-c("A","G","T","C")
names(nucdict)<-c("T","C","A","G")
revcomp<-function(str){
  paste0(rev(nucdict[unlist(strsplit(str,split = ""))]),collapse="")
}

#open barcodes
if(species =='human'){
  barcodes<-fread("Hum_plate_map.txt")
} else if(species == 'mouse'){
  barcodes<-fread("Mus_plate_map.txt")
}
barcodes$search_str <- map_chr(barcodes$Seq, ~revcomp(substr(.,4,26)))
barcodes <- barcodes %>%
  select(Name, type, search_str) %>%
  pivot_wider(names_from = 'type',values_from = 'search_str') %>%
  rename(well = Name)

plate_map <- read.csv('384_well_map.csv')
colnames(plate_map) <- c('well','well_384','plate')

find_plates <- function(path, pltmsk){
  fnames<-grep("align.txt",list.files(path,pattern = pltmsk,full.names = T),value=T)
  A1<-grep("A1_S[0-9]{1,2}_L00[1-4]_align.txt",fnames,value = T)
  A2<-grep("A2_S[0-9]{1,2}_L00[1-4]_align.txt",fnames,value = T)
  B1<-grep("B1_S[0-9]{1,2}_L00[1-4]_align.txt",fnames,value = T)
  B2<-grep("B2_S[0-9]{1,2}_L00[1-4]_align.txt",fnames,value = T)
  plate_files <- c(A1, A2, B1, B2)
  return(plate_files)
}

clean_plate <- function(df) {
  
  out_df <- filter_at(df, vars(c("minQualFR3", "minQualCDR3")), all_vars(. > 20)) %>%
    filter(grepl('^C', aaSeqCDR3) & grepl('F$', aaSeqCDR3)) %>%
    mutate(non_functional = case_when(grepl('[*]', aaSeqCDR3) ~ TRUE,
                                      grepl('_', aaSeqCDR3) ~ TRUE, 
                                      TRUE ~ FALSE)
           ) %>%
    mutate(v_gene = str_split(allVHitsWithScore, '[*]', simplify = T)[,1] )%>%
    mutate(j_gene = str_split(allJHitsWithScore, '[*]', simplify = T)[,1] ) %>%
    select_at(vars(-contains("Alignments"), -contains('HitsWithScore'))) %>%
    select(-targetQualities, -refPoints)
  
  return(out_df)

}

demultiplex_plate <- function(df){
  future({ #multisession run
    dclones <- list()
    for(i in seq(nrow(barcodes))){
      if(i == 1){
        #make temp df that gets smaller after each iteration
        tmp_df <- df
      } 
      well_mask <- which(grepl(barcodes[[i,'alpha']], tmp_df$targetSequences) | grepl(barcodes[[i,'beta']], tmp_df$targetSequences))
      dclones[[i]] <- tmp_df[well_mask,] %>% 
        mutate(well = barcodes[[i,'well']])  %>%
        mutate(chain = ifelse(grepl('TRBV',allVHitsWithScore),'beta','alpha'))  
      tmp_df <- tmp_df[-well_mask,]
    }
    out_df <- do.call(rbind, dclones)
    return(out_df)
  })
}

assign_cdr_spot <- function(df, chain){
  if (chain =='alpha') {
    ifelse( df$non_functional == TRUE, 'cdr3a2','cdr3a')
  } else if (chain =='beta') {
    ifelse( df$non_functional == TRUE, 'cdr3b2','cdr3b')
  }
  
}

#log_threshold is intended for us with raw =FALSE but could be used on raw by setting raw to false
call_well_chains <- function(df, raw = TRUE, log_threshold = 1){
  out_df <- df %>% 
    select(-n) %>%
    pivot_wider(values_from = 'log_n', names_from = 'chain', values_fill = 0)
  
  if (raw ==TRUE){
    
    
    out_df <- out_df %>%
    mutate(
      alpha_pass = ifelse(alpha > thresholds['alpha'] , 'TRUE','FALSE' ),
      beta_pass = ifelse(beta > thresholds['beta'] , 'TRUE','FALSE' )
    ) 
  } else {
    
    out_df <- out_df %>%
      mutate(
        alpha_pass = ifelse(alpha > log_threshold , 'TRUE','FALSE' ),
        beta_pass = ifelse(beta  > log_threshold , 'TRUE','FALSE' )
      ) 
    
  }
    
  out_df <- out_df %>%
    mutate(
      well_call = as.factor(
        case_when(
        alpha_pass == T & beta_pass ==T ~ 'paired',
        alpha_pass == T & beta_pass ==F ~ 'alpha_only',
        alpha_pass == F & beta_pass ==T ~ 'beta_only',
        alpha_pass == F & beta_pass ==F ~ 'empty'
      )
      )
    )
  if (raw ==TRUE){
    out_df <- rename(out_df, raw_well_call = well_call)
  } else {
    out_df <- rename(out_df, filtered_well_call = well_call)
  }
  return(out_df)
} 



# running ====

#find and read
plate_alignments <- map(find_plates(input_path, patterns), ~fread(., nThread = getDTthreads()))

for(j in seq_along(plate_alignments)){
  plate_alignments[[j]]$plate <- case_when(
    j == 1 ~ "A1",
    j == 2~  "A2",
    j == 3~  "B1",
    j == 4~  "B2")
}

#demultiplexing reads by well using future for each 96 well quadrant
demultiplexed_futures <- map(plate_alignments, ~demultiplex_plate(.))
# take resolved futures and make into one table and join with 384 well plate map
demultiplexed_clones <- map(demultiplexed_futures, ~value(.)) %>% 
  do.call(bind_rows,. ) %>% 
  right_join(., plate_map)

#demultiplexed_clones_backup <- demultiplexed_clones

#count number of alignments per well
raw_chain_reads_well <- demultiplexed_clones %>%
  select(well_384, chain) %>%
  group_by_all() %>%
  tally() %>%
  distinct_all() %>%
  mutate(row = str_extract(well_384,'[A-Z]'),
         column = str_extract(well_384,'[0-9]{1,2}')) %>%
  mutate(row = factor(row, levels = rev(LETTERS[1:16]))) %>%
  mutate(column = factor(column, levels = seq(24))) %>%
  mutate(log_n = log(n)) %>%
  mutate(filtered = 'unfiltered')
  
# clean plates to remove low-quality alignments, mark nonfuncitonal, call v/j
demultiplexed_clones_clean <-  clean_plate(demultiplexed_clones)

clean_chain_reads_well <- demultiplexed_clones_clean %>%
  select(well_384, chain) %>%
  group_by_all() %>%
  tally() %>%
  distinct_all() %>%
  mutate(row = str_extract(well_384,'[A-Z]'),
         column = str_extract(well_384,'[0-9]{1,2}')) %>%
  mutate(row = factor(row, levels = rev(LETTERS[1:16]))) %>%
  mutate(column = factor(column, levels = seq(24))) %>%
  mutate(log_n = log(n)) %>%
  mutate(filtered = 'filtered')


wide_chain_reads <- bind_rows(raw_chain_reads_well, clean_chain_reads_well) %>%
  select(-c(3:5)) %>%
  pivot_wider(names_from = 'filtered', values_from = 'log_n',values_fill = 0)


# set threshold for cutoff on each chain independently
ccounts <- wide_chain_reads %>% 
  select(chain, unfiltered) %>%
  group_by(chain) %>% 
  group_split() #split counts to chains

thresholds <- c()
for ( z in seq_along(ccounts)){
  d <- density(ccounts[[z]]$unfiltered) #density of read counts per well
  densityY <- d$y
  densityX <- d$x
  
  #find peak 1
  i <- 1
  while(densityY[i+1] > densityY[i] ) {
    i = i+1
  }
  peak1 <- densityX[i]
  
  #find peak 2
  j <- length(densityY)  
  while(densityY[j-1] > densityY[j] ) {
    j = j-1
  }
  peak2 <- densityX[j]
  
  #find min
  threshold <- densityX[which(densityY == min(densityY[densityX > peak1 & densityX < peak2]))]
  
  thresholds[z] <- threshold
  names(thresholds)[z] <- ifelse(z ==1, 'alpha','beta')
}


# call the well chains with cleaned alignments
raw_well_calls <- call_well_chains(raw_chain_reads_well) %>% select(-filtered)
colnames(raw_well_calls)[4:7] <-paste0('raw_',colnames(raw_well_calls)[4:7])
clean_well_calls <- call_well_chains(clean_chain_reads_well, raw = F, log_threshold = 0)
colnames(clean_well_calls)[4:7] <-paste0('clean_',colnames(clean_well_calls)[4:7])

well_calls <- left_join(raw_well_calls, clean_well_calls)

demultiplexed_clones_clean <- left_join(demultiplexed_clones_clean,well_calls) 
chain_reads_well_df <- bind_rows(raw_chain_reads_well, clean_chain_reads_well)

#saving demultiplexed reads
demulti_name <- paste0(output_prefix,'_demultiplexed_raw.tsv')
write_tsv(demultiplexed_clones, demulti_name)

demulti_name2 <- paste0(output_prefix,'_demultiplexed_clean.tsv')
write_tsv(demultiplexed_clones_clean, demulti_name2)


#plotting ====
total_reads <- sum(sapply(plate_alignments, nrow))
demultiplexed_reads <- nrow(demultiplexed_clones)
pct_read <- paste0('(', round( (demultiplexed_reads/total_reads)*100,1),'%)')
p1_name <- paste(output_prefix,'_plate_read_count_heatmap.pdf')
p2_name <- paste(output_prefix,'_plate_read_count_distribution.pdf')
p3_name <- paste(output_prefix,'_plate_chain_calls.pdf')

p1 <- ggplot(filter(chain_reads_well_df, !is.na(chain)), aes(column,row, fill = log_n )) + 
  geom_tile(color = 'black') + 
  scale_fill_viridis_b() +
  facet_wrap(filtered~chain) +
  labs(title='Number aligned chains by chain for each well',
       subtitle = paste(demultiplexed_reads,' / ', total_reads,' ',pct_read,' alignments demultiplexed')) + 
  theme_bw()

chain_pals <- c('red','black')
names(chain_pals) <- c('alpha','beta')
p2 <- ggplot() + 
  geom_density(data = filter(raw_chain_reads_well, chain =='alpha'),aes(log_n, color = chain , fill = chain ), alpha = 0.4) +
  geom_density(data = filter(raw_chain_reads_well, chain =='beta'), aes(log_n, color = chain, fill = chain), alpha = 0.4) +
  geom_vline(xintercept = thresholds['alpha'], color='red', linetype = 'dashed') +
  geom_vline(xintercept = thresholds['beta'],color='black', linetype = 'dashed') +
  labs(title='Desnity of chain counts by well across plate') + 
  scale_x_continuous(limits = c(0,max(raw_chain_reads_well$log_n))) + 
  scale_color_manual(values = chain_pals) + 
  scale_fill_manual(values = chain_pals) +
  theme_bw()

well_cpal <- c('green','red','blue', 'grey80')
names(well_cpal) <- c('paired','alpha_only','beta_only', 'empty')

p3 <- well_calls %>%
  select(1:3,raw_well_call, filtered_well_call) %>%
  pivot_longer(raw_well_call:filtered_well_call,values_to = 'well_call', names_to = 'filtered') %>%
  ggplot(., aes(column,row, fill = well_call )) + 
  geom_tile(color = 'black') + 
  labs(title='Chains detected in wells') +
  scale_fill_manual(values = well_cpal) + 
  theme_bw() +
  facet_wrap(~filtered, ncol =1)


ggsave(p1_name, p1, width = 20, height = 10)
ggsave(p2_name, p2, width = 8, height = 6)
ggsave(p3_name, p3, width = 10, height = 10)


# this is the algo for assigning chains per well and calling first and second chains
# read the comments within the function for more info
call_chains_per_well <- function(demultiplexed_clones_clean, filtered_call =  TRUE) {
  
  well_clones <- list()
  for(i in seq(plate_map$well_384) ){
    
    print(i)
    #i =which(plate_map$well_384 =='B1')
    df <-  demultiplexed_clones_clean %>% 
      filter(well_384 == plate_map$well_384[i]) %>%
      select(well_384,plate,chain,v_gene,j_gene, aaSeqCDR3, nSeqCDR3, non_functional, raw_well_call, filtered_well_call, targetSequences)
    
    if (filtered_call ==TRUE){
      df <- rename(df, well_call = filtered_well_call)
    } else {
      df <- rename(df, well_call = raw_well_call)
    }
    
    
    df <- df %>%
      filter(well_call != 'empty')
    
    if(nrow(df)== 0){
      next
    } else {
      
      chains <- df %>%
        group_by(chain, nSeqCDR3) %>%
        add_tally(name = 'num_reads_nSeqCDR3')  %>%
        group_by(chain, aaSeqCDR3) %>% 
        add_tally(name = 'num_aaSeqCDR3') %>% 
        ungroup %>%
        distinct_at(.vars = c('aaSeqCDR3'), .keep_all = T)  %>%
        group_by(chain) %>%
        group_split()
      
      cnames <- c()
      for(z in seq_along(chains)){
        cnames <-append( cnames, unique(chains[[z]]$chain))
      }
      names(chains) <- cnames
      
      well_call <- unique(df$well_call)
      
      # here we split paired, alpha_only, beta_only to avoid picking a chain we shouldn't
      # really needs to be wrapped :)
      
      if (well_call == 'paired'){
        #some times you only get one clone for a chain, and the reads are much lower then the other
        # we'll take the top chain but let's make a note of it as questionable
        top_chains_S <- map(chains,~arrange(.,desc(num_reads_nSeqCDR3)) %>% slice_head(n=ifelse(nrow(.)==1,1,2))) 
        
        # if we slice two chains we test how similar those cdr3s are and compare their abundances to make some decisions
        # if the edit distance is large, the genes match or the difference in abundace is too great we pick just the top chain
        if (nrow(top_chains_S[['beta']]) == 2) { 
          bdist <- stringdist(top_chains_S[['beta']]$aaSeqCDR3[1],top_chains_S[['beta']]$aaSeqCDR3[2] )
          v_same <- top_chains_S[['beta']]$v_gene[1] == top_chains_S[['beta']]$v_gene[2]
          j_same <- top_chains_S[['beta']]$j_gene[1] == top_chains_S[['beta']]$j_gene[2]
          same_genes <- v_same & v_same
          log_diff <- log(top_chains_S[['beta']]$num_aaSeqCDR3[1] + 1) / log(top_chains_S[['beta']]$num_aaSeqCDR3[2] + 1)
          #if the two top chains are very similar it's probably pcr error, pick the top one
          if ( bdist <= 4 & same_genes == T ) {
            top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
            top_chains_S[['beta']]$confidence <- 'high'
            top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
          } else { #if theyre far apart
            if( log_diff < 1.5 ){
              top_chains_S[['beta']]$confidence <-  'high'
              if (any(top_chains_S[['beta']]$non_functional)){
                if(all(top_chains_S[['beta']]$non_functional)){
                  top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
                }
                top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta')  # make the functional chain as the top choic
              } else{
                top_chains_S[['beta']]$cdr <- c('cdr3b','cdr3b2')
              }
            } else {
              top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
              top_chains_S[['beta']]$confidence <- 'high'
              top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
            }
          }
          # we can kinda count on the noise to provide at least two chains to compare, 
          # otherwise these are fringe cases that eeked over through the plate cleaning and seem to be very low in abundance
          # in other words we don't expect to see a lot of reads very cleanly to one chain, 
          # we'll mark them so we can take these with a grain of salt and look closer further upstream in the demuliplexed
        } else { 
          top_chains_S[['beta']]$confidence <- ifelse(top_chains_S[['beta']]$num_reads_nSeqCDR3 > 2, 'medium','low')
          top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
          
        }
        
        
        if (nrow(top_chains_S[['alpha']]) == 2) {
          adist <- stringdist(top_chains_S[['alpha']]$aaSeqCDR3[1],top_chains_S[['alpha']]$aaSeqCDR3[2] )
          v_same <- top_chains_S[['alpha']]$v_gene[1] == top_chains_S[['alpha']]$v_gene[2]
          j_same <- top_chains_S[['alpha']]$j_gene[1] == top_chains_S[['alpha']]$j_gene[2]
          same_genes <- v_same & v_same
          log_diff <- log(top_chains_S[['alpha']]$num_aaSeqCDR3[1] + 1) / log(top_chains_S[['alpha']]$num_aaSeqCDR3[2] + 1)
          if (adist <= 4 & same_genes == T) {
            top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
            top_chains_S[['alpha']]$confidence <- 'high'
            top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
          } else { 
            if( log_diff < 1.5 ){
              top_chains_S[['alpha']]$confidence <-  'high'
              if (any(top_chains_S[['alpha']]$non_functional)){
                if(all(top_chains_S[['alpha']]$non_functional)){
                  top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
                }
                top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
              } else{
                top_chains_S[['alpha']]$cdr <- c('cdr3a','cdr3a2')
              }
            } else {
              top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
              top_chains_S[['alpha']]$confidence <- 'high'
              top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
            }
          }
        } else { 
          top_chains_S[['alpha']]$confidence <- ifelse(top_chains_S[['alpha']]$num_reads_nSeqCDR3 > 2, 'medium','low')
          top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
        }
        
        top_chains <- top_chains_S %>% 
          do.call(bind_rows,.)
        
      } else if (well_call == 'alpha_only'){
        
        top_chains_S <- map(chains['alpha'],~arrange(.,desc(num_reads_nSeqCDR3)) %>% slice_head(n=ifelse(nrow(.)==1,1,2))) 
        
        if (nrow(top_chains_S[['alpha']]) == 2) {
          adist <- stringdist(top_chains_S[['alpha']]$aaSeqCDR3[1],top_chains_S[['alpha']]$aaSeqCDR3[2] )
          v_same <- top_chains_S[['alpha']]$v_gene[1] == top_chains_S[['alpha']]$v_gene[2]
          j_same <- top_chains_S[['alpha']]$j_gene[1] == top_chains_S[['alpha']]$j_gene[2]
          same_genes <- v_same & v_same
          log_diff <- log(top_chains_S[['alpha']]$num_aaSeqCDR3[1] + 1) / log(top_chains_S[['alpha']]$num_aaSeqCDR3[2] + 1)
          if (adist <= 4 & same_genes == T) {
            top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
            top_chains_S[['alpha']]$confidence <- 'high'
            top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
          } else { 
            if( log_diff < 1.5 ){
              top_chains_S[['alpha']]$confidence <-  'high'
              if (any(top_chains_S[['alpha']]$non_functional)){
                if(all(top_chains_S[['alpha']]$non_functional)){
                  top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
                }
                top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
              } else{
                top_chains_S[['alpha']]$cdr <- c('cdr3a','cdr3a2')
              }
            } else {
              top_chains_S[['alpha']] <- top_chains_S[['alpha']][1,]
              top_chains_S[['alpha']]$confidence <- 'high'
              top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
            }
          }
        } else {
          top_chains_S[['alpha']]$confidence <- ifelse(top_chains_S[['alpha']]$num_reads_nSeqCDR3 > 2, 'medium','low')
          top_chains_S[['alpha']]$cdr <- assign_cdr_spot(top_chains_S[['alpha']], 'alpha') 
        }
        
        top_chains <- top_chains_S[['alpha']]
        
      } else if (well_call == 'beta_only'){
        top_chains_S <- map(chains['beta'],~arrange(.,desc(num_reads_nSeqCDR3)) %>% slice_head(n=ifelse(nrow(.)==1,1,2))) 
        if (nrow(top_chains_S[['beta']]) == 2) {
          bdist <- stringdist(top_chains_S[['beta']]$aaSeqCDR3[1],top_chains_S[['beta']]$aaSeqCDR3[2] )
          v_same <- top_chains_S[['beta']]$v_gene[1] == top_chains_S[['beta']]$v_gene[2]
          j_same <- top_chains_S[['beta']]$j_gene[1] == top_chains_S[['beta']]$j_gene[2]
          same_genes <- v_same & v_same
          log_diff <- log(top_chains_S[['beta']]$num_aaSeqCDR3[1] + 1) / log(top_chains_S[['beta']]$num_aaSeqCDR3[2] + 1)
          if (bdist <= 4 & same_genes == T) {
            top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
            top_chains_S[['beta']]$confidence <- 'high'
            top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
          } else { 
            if( log_diff < 1.5 ){
              top_chains_S[['beta']]$confidence <-  'high'
              if (any(top_chains_S[['beta']]$non_functional)){
                if(all(top_chains_S[['beta']]$non_functional)){
                  top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
                }
                top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
              } else{
                top_chains_S[['beta']]$cdr <- c('cdr3b','cdr3b2')
              }
            } else {
              top_chains_S[['beta']] <- top_chains_S[['beta']][1,]
              top_chains_S[['beta']]$confidence <- 'high'
              top_chains_S[['beta']]$cdr <- assign_cdr_spot(top_chains_S[['beta']], 'beta') 
            }
          }
        } else {
          top_chains_S[['beta']]$confidence <- ifelse(top_chains_S[['beta']]$num_reads_nSeqCDR3 > 2, 'medium','low')
          top_chains_S[['beta']]$cdr <- ifelse( top_chains_S[['beta']]$non_functional == TRUE, 'cdr3b2','cdr3b')
        }
        top_chains <- top_chains_S[['beta']]
      }
      #chains for the well
      well_clones[[i]] <- top_chains
    }
  }
  #chains for the plate
  out_df <- do.call(bind_rows, well_clones)
  return(out_df)
} 


parsed_clones <- call_chains_per_well(demultiplexed_clones_clean, filtered_call = F)

fname <- paste(output_prefix,'_parsed_clones.tsv')
write_tsv(parsed_clones, fname)

# below here is experimenting with correcting for false pairs and filling in missing chains. Not totally sure what it does anymore. 

clone_table <- parsed_clones %>% 
  select(well_384,well_call, cdr, aaSeqCDR3) %>% 
  pivot_wider(values_from = 'aaSeqCDR3', names_from = 'cdr') %>% 
  select( well_384,well_call,  cdr3a, cdr3a2, cdr3b, cdr3b2)


paired_clones <- clone_table %>% 
  filter((!is.na(cdr3a)| !is.na(cdr3a2)) & (!is.na(cdr3b)| !is.na(cdr3b2)))

tallied <- paired_clones %>%
  group_by(cdr3a,cdr3a2,cdr3b, cdr3b2) %>% 
  add_tally(name = 'clone_size') %>%
  distinct_at(.vars = c('cdr3a','cdr3b','cdr3a2','cdr3b2'), .keep_all = T) %>%
  group_by(cdr3a) %>%
  add_tally(name = 'cdr3a_n_clones') %>%
  group_by(cdr3a2) %>%
  add_tally(name = 'cdr3a2_n_clones') %>%
  group_by(cdr3b) %>%
  add_tally(name = 'cdr3b_n_clones')  %>%
  group_by(cdr3b2) %>%
  add_tally(name = 'cdr3b2_n_clones')

alpha_list <- list()
beta_list <- list()

for (i in seq(nrow(clone_table))){
  
  if( !is.na(clone_table[[i,'cdr3b']]) & !is.na(clone_table[[i,'cdr3a']])){
    if (clone_table[[i,'cdr3a']] %notin% names(alpha_list) & !is.na(clone_table[[i,'cdr3b']])){
      alpha_list[[clone_table[[i,'cdr3a']]]] <- clone_table[[i,'cdr3b']]
    } else {
      alpha_list[[clone_table[[i,'cdr3a']]]] <- append( alpha_list[[clone_table[[i,'cdr3a']]]], clone_table[[i,'cdr3b']] )
    }
    if (clone_table[[i,'cdr3b']] %notin% names(beta_list) ){
      beta_list[[clone_table[[i,'cdr3b']]]] <- clone_table[[i,'cdr3a']]
    } else {
      beta_list[[clone_table[[i,'cdr3b']]]] <- append( beta_list[[clone_table[[i,'cdr3b']]]], clone_table[[i,'cdr3a']] )
    }
  }

}
alpha_list<- map(alpha_list,~sort(table(.),decreasing = T))

beta_list<- map(beta_list,~sort(table(.),decreasing = T))

alpha_beta_list <- list()
for(i in seq_along(beta_list)){
  
  if( names(beta_list[[i]][1]) %in% names(alpha_list)) { # if the paired alpha is in the alpha list
    if ( length(names(alpha_list[[names(beta_list[[i]])[1]]])) > 1 ){ # does the largest clone match in each
      if (all(alpha_list[[names(beta_list[[i]])[1]]]) == 1){
        if (any(alpha_list[[names(beta_list[[i]])[1]]]) %in% names(alpha_list) ) {
          alpha_beta_list[[i]] <- c( names(beta_list[[i]][1])[which.max(beta_list[[i]][1])], names(beta_list)[i]) # fix later
        } else {
          alpha_beta_list[[i]] <- c( names(beta_list[[i]][1])[which.max(beta_list[[i]][1])], names(beta_list)[i])
        }
      } else if () {
      
      names(beta_list)[i] == 
        alpha_beta_list[[i]] <- c( beta_list[[i]][1], names(beta_list)[i])
      
  }
  
}


clean_clones <- tallied %>% filter(cdr3a_n_clones == 1 & cdr3b_n_clones == 1)

troublemakers <- tallied %>% filter(cdr3a_n_clones > 1 | cdr3b_n_clones > 1)

contaminated_clones <- tallied %>% filter(xor(cdr3a_n_clones == 1, cdr3b_n_clones == 1))

betas_res <- troblemakers %>% filter( beta_n_clones != 1) %>%
  group_by(cdr3b) %>% 
  group_split() %>% 
  map(~slice_max(., order_by = clone_size) )
