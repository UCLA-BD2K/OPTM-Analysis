# Specify the input file
input_file = "optm_Lys2AAA.csv"

library(tidyverse); library(stringr);

# peptide location search result: PIR - Multiple Peptide Match
# https://research.bioinformatics.udel.edu/peptidematch/batchpeptidematch.jsp
multiple_peptide_match_data = read_tsv("result-file-from-PIR.txt");

# read merged output from step2
ptm_data = read_csv(input_file);

ptm_dta_dw_0 = 
  ptm_data %>% # mark modification to *
  mutate(Sequence_cleaned_mod = (Sequence %>% str_replace_all("\\(14\\.9632\\)", "*"))) %>% 
  mutate(Sequence_cleaned_mod = (Sequence_cleaned_mod %>% str_replace_all("\\.", ""))) %>% 
  mutate(Sequence_cleaned_mod = (Sequence_cleaned_mod %>% str_replace_all("-", ""))) %>% 
  mutate(Sequence_cleaned_mod = (Sequence_cleaned_mod %>% str_replace_all("=", "")));

ptm_dta_dw_1 =
  ptm_dta_dw_0 %>%
  mutate(locates = (Sequence_cleaned_mod %>% str_locate_all("\\*") %>% paste())) %>%
  mutate(locates = locates %>% str_replace("c", "")) %>%
  mutate(locates = locates %>% str_replace("\\(", "")) %>%
  mutate(locates = locates %>% str_replace("\\)", "")) %>%
  separate("locates",
           into = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'),
           sep = ", ") %>%
  mutate(num_ptm = apply(!is.na(cbind(a, b, c, d, e, f, g, h)), MARGIN = 1, FUN = sum)) %>%
  mutate(num_ptm = (num_ptm/2)) %>% 
  mutate(
    a = as.numeric(a),
    b = as.numeric(b),
    c = as.numeric(c),
    d = as.numeric(d),
    e = as.numeric(e),
    f = as.numeric(f),
    g = as.numeric(g),
    h = as.numeric(h));

ptm_dta_dw =
  bind_rows(
    ptm_dta_dw_1 %>%
      filter(num_ptm == 1) %>%
      mutate(a = a - 1,
             b = b - 1),
    ptm_dta_dw_1 %>%
      filter(num_ptm == 2) %>%
      mutate(
        a = a - 1,
        b = b - 2,
        c = c - 1,
        d = d - 2
      ),
    ptm_dta_dw_1 %>%
      filter(num_ptm == 3) %>%
      mutate(
        a = a - 1,
        b = b - 2,
        c = c - 3,
        d = d - 1,
        e = e - 2,
        f = f - 3
      ),
    ptm_dta_dw_1 %>%
      filter(num_ptm == 4) %>%
      mutate(
        a = a - 1,
        b = b - 2,
        c = c - 3,
        d = d - 4,
        e = e - 1,
        f = f - 2,
        g = g - 3,
        h = h - 4
      )) %>%
  gather(key = abcs, value = ptm_positions, a:h) %>%
  filter(!is.na(ptm_positions)) %>%
  select(-abcs) %>%
  # unique() %>%
  unite(peptide_uniprot, Sequence_cleaned, uniprot, sep = "_");

# mutate(ptm_positions = as.numeric(ptm_positions), ones = 1) %>%
# mutate(ptm_positions = ptm_positions - ones);

ptm_dta_location =
  ptm_dta_dw %>%
  left_join(multiple_peptide_match_data %>%
              unite(peptide_uniprot, `#Query Peptide`, `Protein AC`, sep = "_"),
            by = "peptide_uniprot") %>%
  select(
    Sequence,
    peptide_uniprot,
    ptm_peptide = Sequence_cleaned_mod,
    ptm_position = ptm_positions,
    seq_range = `Matched Range(s)`,
    uniprot_strain_group_timepoint_Sequence_cleaned,
    num_ptm,
    ptm_redundancy = Redundancy,
    total_redundancy = Total_Redundancy,
    protein_length = Length,
    protein_name = `Protein Name`
  );

###### ptm_dta_location_dw$i %>% unique()
ptm_dta_location_dw =
  ptm_dta_location %>%
  separate(
    "seq_range",
    into = c('a'),
    # into = c('a', 'b'),
    #into = c('a', 'b'),
    # into = c('a', 'b', 'c'),
    sep = ", "
  ) %>%
  gather(key = abcs, value = seq_range, a) %>%
  # gather(key = abcs, value = seq_range, a:g) %>%
  # gather(key = abcs, value = seq_range, a:m) %>%
  # gather(key = abcs, value = seq_range, a:c) %>%
  filter(!is.na(seq_range)) %>%
  separate(seq_range,
           into = c("start_seq_range", "end_seq_range"),
           sep = "-") %>%
  mutate(
    start_seq_range = as.numeric(start_seq_range),
    ptm_position = as.numeric(ptm_position),
    ones = 1
  ) %>%
  mutate(ptm_location = (start_seq_range - ones + ptm_position)) %>%
  separate("peptide_uniprot",
           into = c("peptide", "uniprot"),
           sep = "_") %>%
  separate("uniprot_strain_group_timepoint_Sequence_cleaned",
           into = c("UUP", "strain", "group", "timepoint", "Sequence_cleaned"),
           sep = "_") %>% 
  # unite("sample", c("strain", "group", "timepoint"), sep = "_" , remove = F) %>% 
  unite(uniprot_subject_replicate_treatment_ptm_location,
        uniprot,
        strain,
        timepoint,
        group,
        ptm_location,
        sep = "_");

# abundance ratio
gygi_abundance_dw =
  ptm_dta_location_dw %>%
  select(uniprot_subject_replicate_treatment_ptm_location,
         ptm_redundancy,
         total_redundancy) %>%
  group_by(uniprot_subject_replicate_treatment_ptm_location) %>%
  summarise(
    ptm_gygi_abundance = sum(ptm_redundancy),
    total_gygi_abundance = sum(total_redundancy)
  ) %>%
  mutate(abundance_ratio = ptm_gygi_abundance / total_gygi_abundance);

gygi_abundance =
  gygi_abundance_dw %>%
  left_join(
    ptm_dta_location_dw %>%
      select(uniprot_subject_replicate_treatment_ptm_location, protein_length, protein_name) %>% 
      unique(),
    by = "uniprot_subject_replicate_treatment_ptm_location"
  )

gygi_abundance = 
  gygi_abundance %>% 
  # select(-ptm_gygi_abundance, -total_gygi_abundance) %>% 
  mutate(ptm_gygi_abundance = ptm_gygi_abundance/2,
         total_gygi_abundance = total_gygi_abundance/2) %>% 
  separate(uniprot_subject_replicate_treatment_ptm_location, into = c("uniprot", "strain", "timepoint", "treatment", "ptm_location"), sep = "_");

gygi_abundance %>%
  write_csv("gygi_optm_result.csv");


