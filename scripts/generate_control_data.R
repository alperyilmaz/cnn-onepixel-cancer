suppressPackageStartupMessages(library(tidyverse))

tibble("control"=rep(0,1024)) %>% 
  mutate(control=ifelse(row_number()==1,1,control)) %>% 
  mutate(control=ifelse(row_number()==1024,1024,control)) %>% 
  mutate(control=ifelse(row_number() %% 20 == 0,row_number() * 1000,control)) %>%   
  write_tsv("processed-data/control_test.tsv")
