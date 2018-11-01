# Specify the amino acid, mass shit, and output_file name 
# Following the same format as this example
mass_shift = "\\(14\\.9632\\)"
site_mass_shift = "\\K(14\\.9632\\)"
output_file = "optm_Lys2AAA.csv"

library(tidyverse); library(stringr);
dta_files = list.files(path = ".", pattern = "_Dta.csv");

# merge results
dta_data = tibble();
for(each_f_n in dta_files){
  each_data = read_csv(each_f_n);
  dta_data = bind_rows(dta_data,
                       each_data);
}

# remove reversed peptide
dta_data_no_reverse = 
  dta_data %>% 
  filter(!(Locus %>% str_detect("Reverse_sp")));

### to get total redundancy
dat_data_dw = 
  dta_data_no_reverse %>% 
  mutate(Sequence_cleaned = (Sequence %>% str_replace_all(mass_shift, ""))) %>% 
  mutate(Sequence_cleaned = (Sequence_cleaned %>% str_replace_all("\\.", ""))) %>% 
  mutate(Sequence_cleaned = (Sequence_cleaned %>% str_replace_all("-", ""))) %>% 
  mutate(Sequence_cleaned = (Sequence_cleaned %>% str_replace_all("=", "")));

dat_data_ftr = 
  dat_data_dw %>% 
  separate("Locus", into = c("direction", "uniprot", "tax"), sep = "\\|") %>% 
  separate("FileName", into = c("strain", "group", "tissue", "organelle", "timepoint", "fraction"), sep = "_") %>% 
  filter(direction == "sp") %>% 
  unite(uniprot_strain_group_timepoint_Sequence_cleaned, c(uniprot, strain, group, timepoint,  Sequence_cleaned), sep = "_", remove = F);

dat_data_ftr_total_redundancy = 
  dat_data_ftr %>% 
  left_join(dat_data_ftr %>% 
              select(uniprot_strain_group_timepoint_Sequence_cleaned, Redundancy) %>%
              group_by(uniprot_strain_group_timepoint_Sequence_cleaned) %>% 
              summarise(Total_Redundancy = sum(Redundancy)),
            by = "uniprot_strain_group_timepoint_Sequence_cleaned");

optm =
  dat_data_ftr_total_redundancy %>% 
  filter(Sequence %>% str_detect(site_mass_shift));

optm_impt = 
  optm %>% 
  select(uniprot_strain_group_timepoint_Sequence_cleaned, Redundancy, Total_Redundancy, Sequence, 
         uniprot, Sequence_cleaned, strain, group, tissue, timepoint);

optm_impt %>% 
  write_csv(output_file);





