rm(list = ls(all = T));

install.packages(c("tidyverse","stringr","progress"))
library(tidyverse); library(stringr); library(progress);

dta_files = list.files(path = ".", pattern = ".txt"); dta_files %>% length();

## progress bar: elapsed time
pb <- progress::progress_bar$new(
  format = " progressing [:bar] :percent in :elapsed",
  total = length(dta_files), clear = FALSE, width= 60);

for(each_f_n in dta_files){
  raw_data = read_lines(each_f_n, skip = 10);
  
  raw_pro_pep_data = 
    raw_data[(raw_data %>% str_detect("Locus\tSequence")) | # protein header
               (raw_data %>% str_detect("Unique\tFileName")) | # peptide header
               (raw_data %>% str_detect("\\|.*\\|")) | # protein lines
               (raw_data %>% str_detect("_Heart_"))]; # peptide lines
  
  pro_pep_idx = rep(NA, length(raw_pro_pep_data)); pro_pep_idx[c(1,2)] = "search_idx";
  pre_pro_i = 3; # first protein line; static 
  each_protein_idx = 0; # protein idx
  
  for(i in c(3:length(raw_pro_pep_data))){
    each_line = raw_pro_pep_data[i]; # read each line
    
    if(each_line %>% str_detect("\\|.*\\|")){ # each protein line
      if(i - pre_pro_i == 1){ # oh there are multiple protein lines
        pro_pep_idx[i] = each_protein_idx;
      }
      else{ # single protein line
        each_protein_idx = i;
        pro_pep_idx[i] = each_protein_idx;
      }
      pre_pro_i = i; # assign line of previous protein line
    } else{
      pro_pep_idx[i] = each_protein_idx;
    }
  }
  
  raw_pro_pep_data_w_idx = # pro and pep
    raw_pro_pep_data %>% 
    paste(pro_pep_idx, sep = "\t");
  
  pro_data = # pro
    raw_pro_pep_data_w_idx[raw_pro_pep_data_w_idx %>% 
                             str_detect("\\|.*\\|")] %>% # protein lines
    tbl_df() %>% 
    separate(value, into = raw_pro_pep_data_w_idx[1] %>% 
               str_split("\t") %>% 
               .[[1]], 
             sep = "\t");
  
  pep_data = 
    raw_pro_pep_data_w_idx[raw_pro_pep_data_w_idx %>% 
                             str_detect("_Heart_")] %>% # peptide lines
    tbl_df() %>% 
    separate(value, into = raw_pro_pep_data_w_idx[2] %>% 
               str_split("\t") %>% 
               .[[1]], 
             sep = "\t");
  
  dta_data =
    pep_data %>% 
    left_join(pro_data %>% 
                select(Locus, search_idx), 
              by = "search_idx");
  
  pro_data %>% write_csv(paste(str_replace(each_f_n, ".txt",""), "_Pro.csv", sep =""))
  pep_data %>% write_csv(paste(str_replace(each_f_n, ".txt",""), "_Pep.csv", sep =""))
  dta_data %>% write_csv(paste(str_replace(each_f_n, ".txt",""), "_Dta.csv", sep =""))
  
  # show progress
  pb$tick();
  Sys.sleep(1 / 100)
}



