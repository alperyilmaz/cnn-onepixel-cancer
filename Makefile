SHELL := bash
.ONESHELL:

# GPU mode can be set via command line: make GPU=1
# Default to CPU mode if not specified
GPU ?= 0
# default shuffle mode is off
SHUFFLE ?= 0

# folder assignments
RAW=raw-data
PROCESSED=processed-data
RESULT=results
SCRIPT=scripts
MODEL=model

# Docker configuration based on GPU mode
ifeq ($(GPU),1)
    DOCKER_GPU_FLAGS := --gpus all
    DOCKER_IMAGE := alperyilmaz/one-pixel-attack:gpu
else
    DOCKER_GPU_FLAGS :=
    DOCKER_IMAGE := alperyilmaz/one-pixel-attack:cpu
endif

# coloring the output, taken from https://gist.github.com/rsperl/d2dfe88a520968fbc1f49db0a29345b9
ifneq (,$(findstring xterm,${TERM}))
	GREEN        := $(shell tput bold)$(shell tput -Txterm setaf 2)
	BLUE         := $(shell tput bold)$(shell tput -Txterm setaf 4)
	RESET := $(shell tput -Txterm sgr0)
else
	GREEN        := ""
	BLUE         := ""
	RESET        := ""
endif

all: prepare_data train attack
prepare_data: $(RESULT)/selected_genes_sample_transposed.tsv  $(PROCESSED)/sample_labels_matching  $(PROCESSED)/train_sample_labels.tsv $(RESULT)/tcga_gtex_genes_data.npz
train: $(MODEL)/tcgamodel.h5 $(MODEL)/model_evaluation_report.txt
attack:  $(RESULT)/candidate_genes
.PHONY: prepare_data train attack all archive clean_processed clean_all clean_results check

$(PROCESSED)/genes-deseq2: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Get genes that have deseq2 expression..${RESET}"
	@gunzip -c $< | cut -f1 > $@

$(PROCESSED)/genes_high_deg: $(RAW)/gene_attribute_matrix.txt.gz $(PROCESSED)/genes-deseq2 $(PROCESSED)/entrez2ensembl
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Find genes with many up/down cases..${RESET}"
	@echo -e "                        ${BLUE}Then get Ensembl IDs of selected genes and pick the gene if it has deseq2 data..${RESET}"
	source scripts/functions.sh
	@gunzip -c $< | _count_up_down_cases | _get_ensembl | _has_deseq | sort -k3 -nr | awk '{genes[$$2]++;if (length(genes) <= 1024){print $$0}}' > $@

$(PROCESSED)/entrez2ensembl: $(RAW)/gene_attribute_matrix.txt.gz $(RAW)/biomart_gene_conv.tsv
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Find genes that are ensembl-convertible in LINC data..${RESET}"
	@gunzip -c $<  | awk '$$3^2>0 {printf"%s\t%s\n",$$1,$$3}' | awk 'FNR==NR {if($$3>0){genes[$$3]=$$1}; next} ($$2 in genes){printf"%s\t%s\n",$$1,genes[$$2]}' $(filter-out $<,$^) -  > $@
	@#filter-out is for getting second prerequisite name

$(PROCESSED)/sample_labels: $(RAW)/TcgaTargetGTEX_phenotype.txt.gz
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Get sample labels for selected samples..${RESET}"
	@mkdir -p processed-data
	@gunzip -c $< | awk -F"\t" '$$5=="Metastatic" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' > $@
	@gunzip -c $< | awk -F"\t" '$$5=="Normal Tissue" && $$7=="GTEX" {printf"%s\tNormal\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Primary Blood Derived Cancer - Peripheral Blood" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Primary Tumor" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Solid Tissue Normal" && $$7=="TCGA" {printf"%s\tNormal\n",$$1}' >> $@

$(PROCESSED)/max_gene_exp_per_domain: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz $(PROCESSED)/genes_high_deg
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Calculating maximum gene expression..${RESET}"
	source scripts/functions.sh
	@gunzip -c $<  | _filtergenes |  _pivotlonger | _filtersamples | _get_max_tcga_gtex > $@

$(RESULT)/selected_genes_sample_transposed.tsv: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz $(PROCESSED)/sample_labels $(PROCESSED)/genes_high_deg
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Preparing the matrix of selected 1024 genes and selected tumor/normal samples..${RESET}"
	@echo -e "                        ${BLUE}WARNING: This step requires around 6GB memory, so please close unnecessary programs..${RESET}"
ifeq ($(SHUFFLE),1)
	@echo -e "${BLUE}                        WARNING2: Shuffling the genes..${RESET}"
endif
	source scripts/functions.sh
	@mkdir -p results
ifeq ($(SHUFFLE),1)
	@gunzip -c $< |  _filtergenes |  _pivotlonger | _filtersamples | awk '{printf"%s\t%s\t%.0f\n",$$1,$$2,$$3}' | _pivotwider | python $(SCRIPT)/shuffle_genes.py > $@
else
	@gunzip -c $< |  _filtergenes |  _pivotlonger | _filtersamples | awk '{printf"%s\t%s\t%.0f\n",$$1,$$2,$$3}' | _pivotwider > $@
endif

$(PROCESSED)/sample_labels_matching: $(RESULT)/selected_genes_sample_transposed.tsv $(PROCESSED)/sample_labels
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Remove samples that mismatch between DeSeq and Phenotype datasets..${RESET}"
	@cut -f1 $< | sed 1d | awk 'FNR==NR {matching[$$1]++; next} ($$1 in matching){printf"%s\t%s\n",$$1,$$2}' - $(PROCESSED)/sample_labels > $@

$(PROCESSED)/train_sample_labels.tsv: $(PROCESSED)/sample_labels_matching $(SCRIPT)/train_test_split.R
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Split samples into train and test..${RESET}"
	@Rscript $(SCRIPT)/train_test_split.R

$(PROCESSED)/train_data.tsv: $(RESULT)/selected_genes_sample_transposed.tsv $(PROCESSED)/train_sample_labels.tsv
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Generate train and test data from sample names..${RESET}"
	@awk 'FNR==NR {data[$$1]=$$0; next} NR > 2 && FNR==1 {print data["Sample"];next}  $$2 in data {print data[$$2]}' $<  $(PROCESSED)/train_sample_labels.tsv > $@
	@awk 'FNR==NR {data[$$1]=$$0; next} NR > 2 && FNR==1 {print data["Sample"];next}  $$2 in data {print data[$$2]}' $<  $(PROCESSED)/test_sample_labels.tsv > $(PROCESSED)/test_data.tsv

$(PROCESSED)/control_test.tsv: $(SCRIPT)/generate_control_data.R $(RESULT)/selected_genes_sample_transposed.tsv
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Generate control data..${RESET}"
	@Rscript $<
	@cat $(PROCESSED)/control_test.tsv | tr "\n" "\t" | cat <(head -1 $(filter-out $<,$^)) - | sed -e 's/\t$$/\n/' > $(PROCESSED)/control_sample
	@# filter-out is for second prerequisite filename

$(RESULT)/tcga_gtex_genes_data.npz: $(SCRIPT)/sample_convert.py $(PROCESSED)/train_data.tsv $(PROCESSED)/control_test.tsv
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Preparing npz file..${RESET}"
	@mkdir -p results
	@python $<

$(MODEL)/tcgamodel.h5: $(SCRIPT)/tcga_gtex_cnn_train.py $(RESULT)/tcga_gtex_genes_data.npz
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Starting the training process..${RESET}"
	@mkdir -p figures model
	@docker run $(DOCKER_GPU_FLAGS) -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf $(DOCKER_IMAGE) python $<
	#@docker run --gpus all -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf alperyilmaz/one-pixel-attack:gpu python $<

$(MODEL)/model_evaluation_report.txt: $(SCRIPT)/evaluate_model.py $(MODEL)/tcgamodel.h5 $(RESULT)/tcga_gtex_genes_data.npz
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Preparing ROC curve and model evaluation report..${RESET}"
	@docker run $(DOCKER_GPU_FLAGS) -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf $(DOCKER_IMAGE) python $<
	#@docker run --gpus all -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf alperyilmaz/one-pixel-attack:gpu python $<
	@touch $@

$(RESULT)/attack_complete: $(SCRIPT)/attack_tcga_all.py
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Starting one pixel attack..${RESET}"
	@[ -d $(RESULT)/attack_results ] || mkdir -p $(RESULT)/attack_results
	@docker run $(DOCKER_GPU_FLAGS) -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf $(DOCKER_IMAGE) python $<
	#@docker run --gpus all -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf alperyilmaz/one-pixel-attack:gpu python $<
	@touch $@

$(RESULT)/attack_summary: $(RESULT)/attack_complete
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Summarizing attack results..${RESET}"
	@echo -e "${GREEN}[ $$(date +'%Y-%m-%d %H:%M:%S') ] WARNING: All prior attack results are combined into one file! ${RESET}"
	@awk -F"," -v OFS="\t" '$$7=="True" {print $$4,$$2,$$5,$$6,$$9,$$10,$$11,$$12,$$13,$$14,$$15,$$16,$$17}' $(RESULT)/attack_results/attack_results* | sort -n -k1,1 | uniq > $@

$(RESULT)/attack_summary_annotated: $(SCRIPT)/extract_attack.py $(RESULT)/attack_summary $(MODEL)/tcgamodel.h5 $(RESULT)/tcga_gtex_genes_data.npz
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Annotating attack results..${RESET}"
	@mkdir -p results/attack_images/
	@docker run $(DOCKER_GPU_FLAGS) -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf $(DOCKER_IMAGE) python $<
	#@docker run --gpus all -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf -w /tf alperyilmaz/one-pixel-attack:gpu python $<

$(PROCESSED)/min_max_gene_exp_per_domain: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz 
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Calculating min and max expression per gene..${RESET}"
	source scripts/functions.sh
	@gunzip -c $<  | _filtergenes |  _pivotlonger | _filtersamples | awk -f scripts/print_min_max.awk > $@

$(RESULT)/candidate_genes: $(PROCESSED)/min_max_gene_exp_per_domain $(RESULT)/attack_summary_annotated
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Checking if successfully attacked gene is within expression limits..${RESET}"
	@awk -f $(SCRIPT)/check_gene_expression_boundary.awk $< $(filter-out $<,$^) > $(PROCESSED)/attack_summary_annotated_boundary
	@awk '$$5 != $$6 && $$7 ~/Within/' $(PROCESSED)/attack_summary_annotated_boundary | sort -u > $@
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] All steps are completed..${RESET}"

archive:
	@echo "Archiving results.."
	@echo "Running make will start new attack.."
	@tar czvf archive/archive_"$$(date +'%Y%m%d')".tar.gz figures/* misc/* model/* results/attack_images/* results/attack_results/* results/attack_summary* results/candidate_genes results/tcga_gtex_genes_data.npz
	@rm results/attack_complete results/attack_images/* results/attack_results/* results/attack_summary* results/candidate_genes

check:
	@echo -e "${BLUE}[ $$(date +'%Y-%m-%d %H:%M:%S') ] Check if necessary programs are available..${RESET}"
	@gawk --version >/dev/null 2>&1 || (echo "ERROR: gawk is required (standard awk is not enough)."; exit 1)
	@Rscript -e "library(tidyverse)" >/dev/null 2>&1 || (echo "ERROR: R with tidyverse required."; exit 1)
	@docker --version >/dev/null 2>&1 || (echo -e "ERROR: docker is required for running Tensorflow scripts.\n       If you prefer not to use docker, please install Python packages listed under misc/requirements.txt along with tensorflow"; exit 1)
	@echo "All dependencies are met. You can run 'make'.."
clean_processed:
	@echo "Cleaning processed data.."
	@echo rm -r assets processed_data/*  variables/ saved_model.pb

clean_results:
	@echo rm

clean_all:
	@echo rm
