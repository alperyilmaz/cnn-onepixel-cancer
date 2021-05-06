SHELL=bash
# folder assignments
RAW=raw-data
PROCESSED=processed-data
RESULT=results
SCRIPT=scripts
MODEL=model

all: prepare_data train
prepare_data: $(RESULT)/selected_genes_sample_transposed.tsv  $(PROCESSED)/sample_labels_matching  $(PROCESSED)/train_sample_labels.tsv $(RESULT)/tcga_gtex_genes_data.npz
train: $(MODEL)/tcgamodel.h5
.PHONY: prepare_data train all

$(PROCESSED)/genes-deseq2: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Get genes that have deseq2 expression.."
	@gunzip -c $< | cut -f1 > $@

$(PROCESSED)/genes_high_deg: $(RAW)/gene_attribute_matrix.txt.gz $(PROCESSED)/genes-deseq2 $(PROCESSED)/entrez2ensembl
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Find genes with many up/down cases.."
	@echo "                        Then get Ensembl IDs of selected genes and pick the gene if it has deseq2 data.."
	@gunzip -c $< | _count_up_down_cases | _get_ensembl | _has_deseq | sort -k3 -nr | awk '{genes[$$2]++;if (length(genes) <= 1024){print $$0}}' > $@

$(PROCESSED)/entrez2ensembl: $(RAW)/gene_attribute_matrix.txt.gz $(RAW)/biomart_gene_conv.tsv
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Find genes that are ensembl-convertible in LINC data"
	@gunzip -c $<  | awk '$$3^2>0 {printf"%s\t%s\n",$$1,$$3}' | awk 'FNR==NR {if($$3>0){genes[$$3]=$$1}; next} ($$2 in genes){printf"%s\t%s\n",$$1,genes[$$2]}' $(filter-out $<,$^) -  > $@
	@#filter-out is for getting second prerequisite name

$(PROCESSED)/sample_labels: $(RAW)/TcgaTargetGTEX_phenotype.txt.gz
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Get sample labels for selected samples"
	@gunzip -c $< | awk -F"\t" '$$5=="Metastatic" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' > $@
	@gunzip -c $< | awk -F"\t" '$$5=="Normal Tissue" && $$7=="GTEX" {printf"%s\tNormal\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Primary Blood Derived Cancer - Peripheral Blood" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Primary Tumor" && $$7=="TCGA" {printf"%s\tTumor\n",$$1}' >> $@
	@gunzip -c $< | awk -F"\t" '$$5=="Solid Tissue Normal" && $$7=="TCGA" {printf"%s\tNormal\n",$$1}' >> $@

$(PROCESSED)/max_gene_exp_per_domain: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz $(PROCESSED)/genes_high_deg
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Calculating maximum gene expression.."
	@gunzip -c $<  | _filtergenes |  _pivotlonger | _filtersamples | _get_max_tcga_gtex > $@

$(RESULT)/selected_genes_sample_transposed.tsv: $(RAW)/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz $(PROCESSED)/sample_labels $(PROCESSED)/genes_high_deg
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Preparing the matrix of selected 1024 genes and selected tumor/normal samples.."
	@echo "                        WARNING: This step requires around 6GB memory, so please close unnecessary programs.."
	@gunzip -c $< |  _filtergenes |  _pivotlonger | _filtersamples | awk '{printf"%s\t%s\t%.0f\n",$$1,$$2,$$3}' | _pivotwider > $@

$(PROCESSED)/sample_labels_matching: $(RESULT)/selected_genes_sample_transposed.tsv $(PROCESSED)/sample_labels
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Remove samples that mismatch between DeSeq and Phenotype datasets.."
	@cut -f1 $< | sed 1d | awk 'FNR==NR {matching[$$1]++; next} ($$1 in matching){printf"%s\t%s\n",$$1,$$2}' - $(PROCESSED)/sample_labels > $@

$(PROCESSED)/train_sample_labels.tsv: $(PROCESSED)/sample_labels_matching $(SCRIPT)/train_test_split.R
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Split samples into train and test.." 
	@Rscript $(SCRIPT)/train_test_split.R

$(PROCESSED)/train_data.tsv: $(RESULT)/selected_genes_sample_transposed.tsv $(PROCESSED)/train_sample_labels.tsv
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Generate train and test data from sample names.."
	@awk 'FNR==NR {data[$$1]=$$0; next} NR > 2 && FNR==1 {print data["Sample"];next}  $$2 in data {print data[$$2]}' $<  $(PROCESSED)/train_sample_labels.tsv > $@
	@awk 'FNR==NR {data[$$1]=$$0; next} NR > 2 && FNR==1 {print data["Sample"];next}  $$2 in data {print data[$$2]}' $<  $(PROCESSED)/test_sample_labels.tsv > $(PROCESSED)/test_data.tsv

$(PROCESSED)/control_test.tsv: $(SCRIPT)/generate_control_data.R $(RESULT)/selected_genes_sample_transposed.tsv
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Generate control data.."
	@Rscript $<
	@cat $(PROCESSED)/control_test.tsv | tr "\n" "\t" | cat <(head -1 $(filter-out $<,$^)) - | sed -e 's/\t$$/\n/' > $(PROCESSED)/control_sample
	@# filter-out is for second prerequisite filename

$(RESULT)/tcga_gtex_genes_data.npz: $(SCRIPT)/sample_convert.py $(PROCESSED)/train_data.tsv $(PROCESSED)/control_test.tsv
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Preparing npz file.."
	@python $<

$(MODEL)/tcgamodel.h5: $(SCRIPT)/tcga_gtex_cnn_train.py $(RESULT)/tcga_gtex_genes_data.npz
	@echo "[ $$(date +'%Y-%m-%d %H:%M:%S') ] Starting the training process.."
	@docker run -it --rm -u $$(id -u):$$(id -g) -v $$(pwd):/tf alperyilmaz/one-pixel-attack python $<
