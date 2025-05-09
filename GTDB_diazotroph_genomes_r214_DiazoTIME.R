# Data analyses for DiazoTIME genome database
# Julian Damashek, Hamilton College
# jdamashe@hamilton.edu, juliandamashek@gmail.com
# First created: June 21, 2022
# Last edited: May 2, 2025


# Load packages -----------------------------------------------------------

# install, load packages - un-comment lines below to install
# install.packages("remotes","tidyverse","egg","naniar","RColorBrewer","patchwork", "R.utils", "openxlsx")
# remotes::install_github("KarstensLab/microshades", dependencies = T)

library(tidyverse)
library(egg)
library(naniar)
library(RColorBrewer)
library(microshades)
library(patchwork)
library(R.utils) # For unzipping the gz files
library(openxlsx) # For opening raw METABOLIC file from Zenodo
theme_set(theme_classic())



# Load data ---------------------------------------------------------------

# METABOLIC output from Zenodo
raw.metabolic <- read.xlsx("data/METABAOLIC_raw_outputs.xlsx", sheet = "FunctionHit")

#load METABOLIC gene metadata

# #load genome taxonomy data from GTDB r214
# input files are the genome metadata file downloaded from GTDB, parsed to just have the 'accession' and 'taxonomy' columns

# Archaea metadata

url.arc <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/ar53_metadata_r214.tsv.gz"
destination.arc <- "data/ar53_metadata_r214.tsv.gz"
download.file(url.arc, destination.arc) 

gunzip("data/ar53_metadata_r214.tsv.gz", remove=FALSE, overwrite=TRUE)
gtdb.meta.arch <- read_tsv("data/ar53_metadata_r214.tsv")

# Bacteria metadata, downloaded from: 
url.bac <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/bac120_metadata_r214.tsv.gz"
destination.bac <- "data/bac120_metadata_r214.tsv.gz"
download.file(url.bac, destination.bac)

gunzip("data/bac120_metadata_r214.tsv.gz", remove=FALSE, overwrite=TRUE)
gtdb.meta.bac <- read_tsv("data/bac120_metadata_r214.tsv")

# combine Archaea and Bacteria metadata together 
gtdb.meta.combined <- rbind(gtdb.meta.arch, gtdb.meta.bac)



# Format files ------------------------------------------------------------

# Edit METABOLIC gene metadata
gene.path <- raw.metabolic %>%
  select("Category", "Function", "Gene.abbreviation") %>%
  arrange((.[["Category"]]))


# GTDB metadata: parse to just taxonomy data, then expand each taxonomy level into its own column
genome.tax <- gtdb.meta.combined %>% dplyr::select(c(accession,gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","species"), sep="(;s__)") %>% 
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","genus"), sep="(;g__)") %>%
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","family"), sep="(;f__)") %>% 
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","order"), sep="(;o__)") %>%
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","class"), sep="(;c__)") %>% 
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","phylum"), sep="(;p__)") %>%
  separate(gtdb_taxonomy, into=c("gtdb_taxonomy","domain"), sep="(d__)") %>% 
  dplyr::select(-c(gtdb_taxonomy)) %>%
  # replace Pseudomonadota with class names so we can see proteobacterial types
  mutate(phylum.no.proteo=phylum, phylum=case_when(phylum=="Pseudomonadota"~class, TRUE~phylum), .before=3) %>%
  # add "Genome" column (= accession column but without the first 3 characters)
  mutate(Genome=gsub('.[^.]*$', '', accession), .before=2, Genome=str_sub(Genome,4,-1))


#combine the different subgroups of Bacillota, Desulfobacterota, and Nitrospirota into single groups
genome.tax.fixed <- genome.tax %>% 
  mutate(phylum=case_when(str_detect(phylum,"^Bacillota")~"Bacillota", TRUE~phylum), 
         phylum=case_when(str_detect(phylum,"^Desulfobacterota")~"Desulfobacterota", TRUE~phylum),
         phylum=case_when(str_detect(phylum,"^Nitrospirota")~"Nitrospirota", TRUE~phylum))




# Make full METABOLIC file for downstream analysis ------------------------

# analyze full METABOLIC output from "Function.Hit" tab of spreadsheet
gtdb.met <- raw.metabolic %>% select(!c("Category", "Function")) %>%
  pivot_longer(cols = !c(Gene.abbreviation), names_to= 'Genome', values_to= 'gene') %>%
  pivot_wider(names_from = Gene.abbreviation, values_from = gene) %>%
  mutate(accession=Genome, 
         Genome=gsub('.[^.]*$', '', accession), .before=1) %>%
  separate_wider_delim(Genome, ".", names = c("Genome", "extra", "extra2")) %>%
  dplyr::select(-c(accession, extra, extra2)) %>% 
  gather(key="Gene.abbreviation", value="presence", -c(Genome)) %>%
  # add METABOLIC pathway metadata and GTDB taxonomy
  left_join(gene.path, by="Gene.abbreviation") %>% 
  left_join(genome.tax.fixed, by="Genome") %>%
  # Remove Fe/Mn reduction genes, which seem to give a lot of false positives. They seem to be hitting a common cytochrome or something. 
  # And remove N2 fixation, which they all have by definition
  # And remove ammonia oxidation, which is hitting C1-oxidation genes
  filter(Function!="Metal (Iron/Manganese) reduction" & 
           Function!="N2 fixation" & 
           Function!="Ammonia oxidation") %>%
  # combined "Methane metabolism" into "C1 metabolism", and "Urea utilization" into "Nitrogen cycling"
  mutate(Category=case_when(Category=="Methane metabolism"~"C1 metabolism", 
                            Category=="Urea utilization"~"Nitrogen cycling", 
                            TRUE~Category))



# Summaries of most abundant genomes at different taxonomic levels --------

# count number of genomes per genus
gtdb.count.by.gen <- gtdb.met %>% 
  distinct(Genome, .keep_all=T) %>% 
  group_by(genus) %>% count(name="n_genus") %>% 
  arrange(desc(n_genus)) %>% ungroup() 

n.genomes <- sum(gtdb.count.by.gen$n_genus) #2798
genus_key <- genome.tax.fixed %>% dplyr::select(c(genus,phylum)) %>% distinct()
gtdb.top.gen <- gtdb.count.by.gen %>% filter(n_genus >= 10) %>% left_join(genus_key, by="genus") %>%
  mutate(genus=reorder(factor(genus),-n_genus), 
         phylum=reorder(factor(phylum),-n_genus))

#count number of genomes per family
gtdb.count.by.fam <- gtdb.met %>% 
  distinct(Genome, .keep_all=T) %>% 
  group_by(family) %>% 
  count(name="n_family") %>% 
  arrange(desc(n_family)) %>% 
  ungroup() 

family_key <- genome.tax.fixed %>% 
  select(c(family,phylum)) %>% 
  distinct()

gtdb.top.fam <- gtdb.count.by.fam %>% 
  filter(n_family >= 20) %>% 
  left_join(family_key, by="family") %>% 
  mutate(family=reorder(factor(family),-n_family), 
         phylum=reorder(factor(phylum),-n_family))



#count number of genomes per phylum
gtdb.count.by.phy <- gtdb.met %>% 
  distinct(Genome, .keep_all=T) %>% 
  group_by(phylum) %>% 
  count(name="n_phylum") %>% 
  arrange(desc(n_phylum)) %>% 
  ungroup() %>% 
  mutate(phylum=reorder(factor(phylum),-n_phylum))

#count number of genomes with each metabolic category present 
gtdb.met.pres <- gtdb.met %>% 
  filter(presence=="Present") %>% 
  mutate_if(is.character, str_replace_all, pattern="As cycling", replacement="Arsenic cycling") %>%
  mutate_if(is.character, str_replace_all, pattern="Oxygen metabolism - c", replacement="C") %>% 
  mutate_if(is.character, str_replace_all, pattern="Methane oxidation - ", replacement="") %>%
  mutate_if(is.character, str_replace_all, pattern=", QoxABCD", replacement="") %>% 
  group_by(Category) %>% 
  mutate(N=n()) %>% 
  ungroup()

#count functions (smaller scale) and categories (larger scale) present per genome
gtdb.func.pres.gen.count <- gtdb.met.pres %>% 
  group_by(Genome,Function) %>% 
  add_count(Genome,Function,name="n_func") %>% #n_func is number of FUNCTIONS in each genome
  distinct(Genome,Function,n_func,.keep_all=T) %>% 
  ungroup() %>% 
  mutate(func_pres=ifelse(n_func>=1,1,0)) %>% #func_pres is 1 if there are ≥1 gene from function in genome, 0 otherwise
  group_by(Genome,Category) %>% 
  add_count(Genome,Category,name="n_cat") %>% 
  mutate(cat_pres=ifelse(n_cat>=1,1,0)) %>% 
  distinct(Genome,Category,cat_pres,.keep_all=T) %>% 
  ungroup() %>% #n_cat is the number of CATEGORIES (broader)
  left_join(gtdb.count.by.phy, by="phylum") %>% 
  arrange(Genome,Category)

#now use these counts to summarize across levels of taxonomy

#phyla
#count number of genomes in each phylum with a pathway function present and calculate proportion of genomes in phylum with that function
gtdb.func.pres.phy <- gtdb.func.pres.gen.count %>% 
  group_by(phylum,Function) %>% 
  add_tally(name="n_genomes_function") %>% 
  ungroup() %>% 
  mutate(perc_phylum_function=100*n_genomes_function/n_phylum) %>% 
  group_by(phylum,Function) %>% 
  ungroup %>% 
  distinct(phylum,Function,.keep_all=T) %>%
  select(c(Category,Function,phylum,n_phylum,n_genomes_function,perc_phylum_function)) %>% 
  mutate(phylum=reorder(factor(phylum),n_phylum)) %>% 
  arrange(Category) %>% 
  arrange(phylum) #switched order of arranges here

#grouped by category (broader)
gtdb.cat.pres.phy <- gtdb.func.pres.gen.count %>% 
  group_by(phylum,Category) %>% 
  add_tally(name="n_genomes_category") %>% 
  ungroup() %>% 
  mutate(perc_phylum_category=100*n_genomes_category/n_phylum) %>% 
  group_by(phylum,Category) %>% 
  ungroup %>% 
  distinct(phylum,Category,.keep_all=T) %>%
  select(c(Category,Function,phylum,n_phylum,n_genomes_category,perc_phylum_category)) %>% 
  mutate(phylum=reorder(factor(phylum),n_phylum)) 

# #count number of GENOMES within each category
gtdb.met.genome.count <- gtdb.func.pres.gen.count %>% 
  group_by(Category) %>% 
  mutate(cat.count=sum(cat_pres)) %>% 
  ungroup() %>%
  group_by(Function) %>% 
  mutate(func.count=sum(func_pres)) %>% 
  ungroup() %>% 
  distinct(Function,func.count,.keep_all=T) %>%
  mutate(cat.perc=100*cat.count/n.genomes) %>% 
  mutate(func.perc=100*func.count/n.genomes) %>% 
  arrange(desc(cat.count))

# colculate % genomes in database but taking out some metabolisms not based on energy
gtdb.met.genome.count.trim <- gtdb.met.genome.count %>% 
  filter(!Category %in% c("Nitrile hydration","Thermophilic specific"))



# Plot METABOLIC categories across genomes --------------------------------

#plot categories with stacked bars for functions

#reorder functions
gtdb.met.genome.count.trim$Function <- factor(gtdb.met.genome.count.trim$Function, levels=c(
  "Acetogenesis",
  "Acetate to acetyl-CoA",
  "Sulfur oxidation",
  "Sulfide oxidation",
  "Sulfite reduction",
  "Sulfate reduction",
  "Thiosulfate oxidation",
  "Thiosulfate disproportionation",
  "Cytochrome c oxidase, caa3-type",
  "Cytochrome c oxidase, cbb3-type",
  "Cytochrome (quinone) oxidase, bd type",
  "Cytochrome (quinone) oxidase, bo type",
  "Ni-Fe Hydrogenase",
  "FeFe hydrogenase",
  "Fe hydrogenase",
  "Nitrate reduction",
  "Nitrite reduction to ammonia",
  "Urease",
  "Nitrite reduction",
  "Nitric oxide reduction",
  "Nitrous oxide reduction",
  "Nitrite oxidation",
  "Formaldehyde oxidation",
  "Methanol oxidation",
  "Soluble methane monoxygenase",
  "Methane production",
  "Partculate methane monooxygenase",
  "Methyl amine -> formaldehyde",
  "CBB cycle - Rubisco (Form I)",
  "CBB cycle - Rubisco (Form II)",
  "Reverse TCA cycle",
  "Arsenate reduction",
  "Arsenite oxidation",
  "Selenate reduction",
  "Perchlorate reduction",
  "Chlorite reduction"),
  labels=c(
    "Acetogenesis",
    "Acetate to acetyl-CoA",
    "Sulfur oxidation",
    "Sulfide oxidation",
    "Sulfite reduction",
    "Sulfate reduction",
    "Thiosulfate oxidation",
    "Thiosulfate disproportionation",
    "Cytochrome c oxidase, caa3-type",
    "Cytochrome c oxidase, cbb3-type",
    "Cytochrome (quinone) oxidase, bd type",
    "Cytochrome (quinone) oxidase, bo type",
    "Ni-Fe hydrogenase",
    "Fe-Fe hydrogenase",
    "Fe hydrogenase",
    "Nitrate reduction",
    "Nitrite reduction to ammonia",
    "Urease",
    "Nitrite reduction",
    "Nitric oxide reduction",
    "Nitrous oxide reduction",
    "Nitrite oxidation",
    "Formaldehyde oxidation",
    "Methanol oxidation",
    "Soluble methane monoxygenase",
    "Methane production",
    "Particulate methane monooxygenase",
    "Methyl amine -> formaldehyde",
    "CBB cycle - Rubisco (Form I)",
    "CBB cycle - Rubisco (Form II)",
    "Reverse TCA cycle",
    "Arsenate reduction",
    "Arsenite oxidation",
    "Selenate reduction",
    "Perchlorate reduction",
    "Chlorite reduction"))

func.cols <- c("#6baed6","#deebf7", #2 fermentation
               "#67000d","#cb181d","#ef3b2c","#fb6a4a","#fcbba1","#fff5f0", #6 S cycling
               "#feb24c","#fed976","#ffeda0","#ffffcc", #4 oxidative phosphorylation
               "#74c476","#c7e9c0","#f7fcf5", #3 hydrogenases
               "#3f007d","#6a51a3","#807dba","#9e9ac8","#bcbddc","#dadaeb","#efedf5", #7 N cycling
               "#d94801","#fd8d3c", "#fdae6b", "#a63603","#f16913","#fdd0a2", #3 C1 metabolism, 3 methane metabolism
               "#00441b","#238b45","#74c476", #3 C fixation
               "#d94801","#fdae6b", #2 arsenic cycling #    "#08306b","#c6dbef"
               "#252525","#969696","#d9d9d9" #3 other random reductions
)

gtdb.met.cat <- ggplot(data=gtdb.met.genome.count.trim, 
                       aes(x=fct_inorder(Category),y=func.perc,fill=Function)) + 
  geom_bar(stat="identity",color="black") + 
  scale_fill_manual(values=func.cols) + 
  scale_y_continuous(name="Genomes with metabolic categories\n(% diazotroph database)", 
                     expand=c(0,0), 
                     limits=c(0,100)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(angle=45,
                                 hjust=1,
                                 #vjust=1.05,
                                 color="black"), 
        axis.ticks.x=element_blank(),
        legend.title=element_blank(), 
        legend.margin=margin(60,0,0,40), 
        legend.text=element_text(size=7, margin=margin(0,30,0,0)))


ggsave(gtdb.met.cat, file="figures/metabolic functions GTDB percent genomes stacked.pdf", height=6.5, width=13.5)

rm(list=ls())

