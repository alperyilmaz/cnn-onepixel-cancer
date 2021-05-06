suppressPackageStartupMessages(library(tidyverse))
PROCESSED <- "processed-data/"

raw_data <- read_tsv(paste0(PROCESSED,"sample_labels_matching"),col_names=F) %>% 
  mutate(tumor=ifelse(X2=="Tumor",1,0),
         n=row_number()) %>% 
  select(n,sample_id=X1,tumor)

set.seed(1)

train <- raw_data %>% 
  group_by(tumor) %>%
  sample_frac(0.8) %>% 
  ungroup() %>%
  sample_frac(1)  # shuffle samples since tumor and non-tumor are bundled together

test <- raw_data %>% 
  anti_join(train, by="n") %>%
  sample_frac(1) # shuffle samples since tumor and non-tumor are bundled together
  
write_tsv(train,paste0(PROCESSED,"train_sample_labels.tsv"))
write_tsv(test,paste0(PROCESSED,"test_sample_labels.tsv"))
