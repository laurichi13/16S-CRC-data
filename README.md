# 16S rRNA sequencing gene datasets for CRC data

###### tags: `16S rRNA` `COST Action` `ML`

### Used datasets 


| Dataset                                                           | 16S rRNA Region    | Control (n) | Adenoma (n) | CRC (n) | Available metadata                                                |
| ----------------------------------------------------------------- | --- | ----------- | ----------- | ------- | ----------------------------------------------------------------- |
| [Baxter](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4823848/)   | V4    | 172         | 198         | 120     | Gender, age, weight, height, BMI, country, race                   |
| [Zackular](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4221363/) |   V4  | 30          | 30          | 30      | Gender, age, weight, height, BMI, country, race, FOBT, medication |
| [Zeller](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299606/)   |    V4 | 75          | 13          | 41      | Gender, age,  BMI, country, FOBT                                  |
| **TOTAL**                                                            |    V4 | 277         | 241         | 191     | *All of the above*                                                |
# Data processing & sharing

All datasets were processed using [qiime2](https://docs.qiime2.org/2021.11/) pipeline with [DADA2](https://benjjneb.github.io/dada2/) for Sequence quality control and feature table construction, and [SILVA](https://www.arb-silva.de/) database for taxonomic assignment, then a phyloseq object was constructed.  


## Download qiime2 artifacts 

[Here](https://drive.google.com/file/d/1ME5o9vIZ1opihPtkVFlc0gk--0lSTkjq/view?usp=share_link) 

## Download processed data :open_file_folder: 

- Abundance table at genus level is [here.](https://drive.google.com/file/d/1w5sFdGfuwug4cQuEqgHcjKlfLbNF6iVn/view?usp=share_link)
:warning: <font color = 'gray'>Sample counts with NO filtering.</font> 
- Clean metadata is [here](https://drive.google.com/file/d/1yvMj_t-XHAa3tIM-Qd59kPCAhBlgKRQh/view?usp=sharing). 
:memo: <font color = 'gray'> Countries: CA - Canada. USA - United States of America. FRA - France. </font> 
- Phyloseq object is [here](https://drive.google.com/file/d/1OZBXD1XB_5N4PbL8At1aksnG0Wns0F65/view?usp=sharing). 

## Steps for processing data

Qiime version: 2021.4 :eyes:

### ASV table construction
```bash=s
conda activate qiime2-2021.4 

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifiest.csv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33 \

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
  
  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-n-threads 6 \
  --p-trim-left-f 22 \
  --p-trim-left-r 22 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```

### Taxonomic assigmnet 

:warning:  Download SILVA database classifier used with Qiime version 2021.4 


```bash=s
  wget https://data.qiime2.org/2021.4/common/silva-138-99-nb-classifier.qza
```
Taxonomic assigmnet
```bash=s
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file list.tsv \
  --o-visualization taxa-bar-plots.qzv
```
### Phyloseq 

```R=1
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(viridis)

#ASV table####
ASV<-read_qza("DATA/table.qza")
#Taxonomy table
taxonomyq <-read_qza("DATA/taxonomy.qza")
#Tax table transformations
taxtableq <- taxonomyq$data %>% as.tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtableq$Kingdom <- gsub("d__", "",taxtableq$Kingdom)
taxtableq$Phylum <- gsub("p__", "",taxtableq$Phylum)
taxtableq$Class <- gsub("c__", "",taxtableq$Class)
taxtableq$Order <- gsub("o__", "",taxtableq$Order)
taxtableq$Family <- gsub("f__", "",taxtableq$Family)
taxtableq$Genus <- gsub("g__", "",taxtableq$Genus)
taxtableq$Species <- gsub("s__", "",taxtableq$Species)

#Clean taxa: Remove Kingdom Unassigned & Eukaryota
taxtableq <- taxtableq[grep("Eukaryota|Unassigned", taxtableq$Kingdom, invert=T),]

#Tree
tree<-read_qza("DATA/rooted-tree.qza")

#Metadata
metadata<-read_csv("DATA/metadata.csv")

#Construct phyloseq object
ASVs <- otu_table(ASV$data, taxa_are_rows = T) 
tree <- phy_tree(tree$data) 
TAXq <- tax_table(as.data.frame(taxtableq) %>% select_("-Confidence") %>% column_to_rownames("Feature.ID") %>% as.matrix()) #moving the taxonomy to the way phyloseq wants it
sample_metadata <- sample_data(metadata %>% as.data.frame()) 
sample_names(sample_metadata) <- paste(metadata$SampleID)
physeq<- merge_phyloseq(ASVs, tree, TAXq,sample_metadata)

#Remove samples with no diagnosis - Some extra samples were inlcuded by mistake-  
physeq <- subset_samples(physeq, Diagnosis!="NA")

#Agglomerate at genus level 
physeq.gen <- tax_glom(physeq, taxrank="Genus")

#Write data for sharing
data <- psmelt(physeq.gen) 
ASV_tab <- data %>% select(SampleID,Abundance,Genus) 
genus <- pivot_wider(ASV_tab, id_cols=SampleID, names_from=Genus, values_from=Abundance, values_fn = mean)
write_csv(genus,"genus.csv")

#Save phyloseq object for sharing
saveRDS(physeq.gen, "physeq.RDS")
```

### Metadata 
Complete metadata found for the three studies is [here](https://docs.google.com/spreadsheets/d/1iJ2JeL-FQrs4Aj9dGJuNusGvHf3-_rWs/edit?usp=sharing&ouid=104693150197955528399&rtpof=true&sd=true). 
:memo: <font color = 'gray'> In "Clean metadata" is only the metadata that is common for the three studies. </font>
