library(tidyverse)
theme_set(theme_classic())
library(lemon)
library(naniar)
library(RColorBrewer)
library(microshades)
library(ggupset)
library(ComplexUpset)
library(patchwork)


#########################################
#
# Add GTDB taxonomy and METABOLIC information to IMG sequences
#
#########################################

#GTDB metadata
gtdb.meta <- read_csv("example_database_use/DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv") %>% 
  #combine the different subgroups of Bacillota, Desulfobacterota, and Nitrospirota into single groups
  mutate(phylum=case_when(str_detect(phylum,"^Bacillota")~"Bacillota", TRUE~phylum), 
         phylum=case_when(str_detect(phylum,"^Desulfobacterota")~"Desulfobacterota", TRUE~phylum),
         phylum=case_when(str_detect(phylum,"^Nitrospirota")~"Nitrospirota", TRUE~phylum))
#remove the ".1" etc from the end of the genome accessions
gtdb.meta <- gtdb.meta %>% mutate(Genome=gsub('.[^.]*$', '', gtdb.meta$Genome))

#IMG metadata
img.meta <- read_csv("OmicsMetadata_Revised_030123.csv")
img.meta$FastaID <- as.character(img.meta$FastaID)

#nifH data
diamond.nifH <- read_tsv("IMG_nifH_diamond.txt", col_names=F)
colnames(diamond.nifH) <- c("qseqid","Genome","evalue","bitscore")
diamond.nifH <- diamond.nifH %>% mutate(gene="NifH", FastaID=as.character(gsub('_.*', '', diamond.nifH$qseqid)))

#nifK data
diamond.nifK <- read_tsv("IMG_nifK_diamond.txt", col_names=F)
colnames(diamond.nifK) <- c("qseqid","Genome","evalue","bitscore")
diamond.nifK <- diamond.nifK %>% mutate(gene="NifK", FastaID=as.character(gsub('_.*', '', diamond.nifK$qseqid)))

#nifD data
diamond.nifD <- read_tsv("IMG_nifD_diamond.txt", col_names=F)
colnames(diamond.nifD) <- c("qseqid","Genome","evalue","bitscore")
diamond.nifD <- diamond.nifD %>% mutate(gene="NifD", FastaID=as.character(gsub('_.*', '', diamond.nifD$qseqid)))

#combine all 3 nif genes together
diamond.all.nif <- rbind(diamond.nifH, diamond.nifK, diamond.nifD) 
diamond.all.nif <- diamond.all.nif %>% left_join(gtdb.meta, by="Genome") %>% left_join(img.meta, by="FastaID", relationship="many-to-many") %>% 
  #get rid of garbage ecosystem subtypes
  filter(!str_detect(Ecosystem_Subtype_NEW,"include"))

#identify low abundance phyla 
nif.count <- nrow(diamond.all.nif)

#group anything < 1% of all nif sequences  to "Low abundance" phylum
phy.counts.nif <- diamond.all.nif %>% group_by(phylum) %>% summarize(n=n()) %>% filter(n>=nif.count*0.0095) %>% arrange(desc(n))
phy.counts.nif.abund <- c(unique(phy.counts.nif$phylum),"Low abund.")
phy.abund.order <- as.data.frame(cbind(phy.counts.nif.abund, factor(phy.counts.nif.abund, levels=phy.counts.nif.abund)))
colnames(phy.abund.order) <- c("phylum","phy_factor")
phy.abund.order$phy_factor <- as.numeric(phy.abund.order$phy_factor)
diamond.all.nif.phy <- diamond.all.nif %>% mutate(phylum=case_when(phylum %in% phy.abund.order$phylum[-max(length(phy.abund.order[,1]))]~phylum, TRUE~"Low abund.")) %>% 
  ungroup()

# img.phy.by.eco.sub <- img.phy.by.eco.sub %>% group_by(Ecosystem_Subtype_NEW) %>% count(name="n_eco_sub") %>% arrange(desc(n_eco_sub)) %>% ungroup() %>% 
# mutate(Ecosystem_Subtype_NEW=reorder(factor(Ecosystem_Subtype_NEW),-n_eco_sub))

#count by subtypes, splitting sample type (aquatic, benthic, etc.)
img.eco.sub.counts.nif <- diamond.all.nif.phy %>% group_by(Ecosystem_Subtype_NEW,Sample_Type) %>% count(name="n_eco_sub") %>% arrange(desc(n_eco_sub)) %>% ungroup() 
img.eco.sub.nif.abund <- unique(img.eco.sub.counts.nif$Ecosystem_Subtype_NEW)
img.eco.sub.nif.order <- as.data.frame(cbind(img.eco.sub.nif.abund, factor(img.eco.sub.nif.abund, levels=img.eco.sub.nif.abund)))
colnames(img.eco.sub.nif.order) <- c("Ecosystem_Subtype_NEW","eco_sub_order")
img.eco.sub.nif.order$eco_sub_order <- as.numeric(img.eco.sub.nif.order$eco_sub_order)

#incorporate phylum, ecosystem subtype, sample type into dataset
img.phy.by.eco.sub <- diamond.all.nif.phy %>% group_by(gene,Ecosystem_Subtype_NEW,Sample_Type,phylum) %>% summarize(Count=n()) %>% 
  left_join(phy.abund.order, by="phylum") %>% left_join(img.eco.sub.nif.order, by="Ecosystem_Subtype_NEW") %>% ungroup() %>% 
  filter(!is.na(eco_sub_order))
#set order of variables for plot
img.phy.by.eco.sub$gene <- factor(img.phy.by.eco.sub$gene, levels=c("NifH","NifD","NifK"))
img.phy.by.eco.sub <- img.phy.by.eco.sub %>% mutate(phylum=reorder(phylum,phy_factor))
img.phy.by.eco.sub <- img.phy.by.eco.sub %>% mutate(Ecosystem_Subtype_NEW=reorder(Ecosystem_Subtype_NEW,eco_sub_order),
                                                    cyanos=case_when(phylum=="Cyanobacteriota"~"Cyanobacteria", TRUE~"Other")) %>% group_by(gene,Ecosystem_Subtype_NEW,Sample_Type) %>%
  mutate(Total=sum(Count)) %>% ungroup()
img.phy.by.eco.sub$Sample_Type <- factor(img.phy.by.eco.sub$Sample_Type, levels=c("Aquatic","Sediment","Biofilm","Plant"))

# #add up nif totals for each ecosystem subtype
# img.nif.totals <- img.phy.by.eco.sub %>% group_by(gene,Ecosystem_Subtype_NEW) %>% summarize(Total=sum(Count))

#plot phyla by ecosystem subtype
#15 colors needed for phyla
col.phy <- c(rev(microshades_palette("micro_cvd_green")), rev(microshades_palette("micro_cvd_orange")[c(2:5)]), 
             rev(microshades_palette("micro_cvd_purple")),microshades_palette("micro_cvd_gray")[c(2)])

img.eco.type.totals <- ggplot(data=distinct(img.phy.by.eco.sub,Ecosystem_Subtype_NEW,Sample_Type,gene,.keep_all=T), aes(x=Ecosystem_Subtype_NEW, y=Total)) + 
  facet_rep_grid(rows=vars(gene),cols=vars(Sample_Type),repeat.tick.labels=T,scales="free_y") + geom_bar(stat="identity") + 
  scale_y_continuous(name="Nif proteins (count)", expand=c(0,0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), axis.title.x=element_blank())
pdf(file="IMG total Nifs by phylum eco and sample type.pdf",width=8, height=8)
print(img.eco.type.totals)
dev.off()

img.by.eco.type <- ggplot(data=img.phy.by.eco.sub, aes(x=Ecosystem_Subtype_NEW, y=100*Count/Total, fill=phylum, color=phylum)) + 
  facet_rep_grid(rows=vars(gene),cols=vars(Sample_Type),repeat.tick.labels=T,scales="free_y") +
  geom_bar(stat="identity", position="stack") + scale_color_manual(values=col.phy, name=NULL) + scale_fill_manual(values=col.phy, name=NULL) +
  scale_y_continuous(name="Environmental proteins (count)", expand=c(0,0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), axis.title.x=element_blank())
pdf(file="IMG nif genes by phylum eco and sample type.pdf",width=8, height=8)
print(img.by.eco.type)
dev.off()

# img.eco.type.cyanos <- ggplot(data=img.phy.by.eco.sub, aes(x=Ecosystem_Subtype_NEW, y=100*Count/Total, fill=cyanos, color=cyanos)) + facet_rep_wrap(~gene) +
#   geom_bar(stat="identity", position="stack") + scale_color_manual(values=c("#66c2a4","#525252"), name=NULL) + scale_fill_manual(values=c("#66c2a4","#525252"), name=NULL) +
#   scale_y_continuous(name="Environmental proteins (count)", expand=c(0,0)) + 
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), axis.title.x=element_blank())
img.eco.type.cyanos <- ggplot(data=filter(img.phy.by.eco.sub,cyanos=="Cyanobacteria"), aes(x=Ecosystem_Subtype_NEW, y=100*Count/Total)) + geom_bar(stat="identity",fill="#66c2a4",color="black") +
  facet_rep_grid(rows=vars(gene),cols=vars(Sample_Type),repeat.tick.labels=T,scales="free_y") +
  scale_y_continuous(name="Cyanobacterial proteins (% Nif count)", expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), axis.title.x=element_blank())

pdf("IMG nif cyanos by phylum eco and sample type.pdf",height=8,width=8)
print(img.eco.type.cyanos)
dev.off()

#plot by METABOLIC data
#sulfur cycling
gtdb.genomes.s <- gtdb.met.pres %>% filter(Category=="Sulfur cycling") %>% select(c(Genome,phylum,class,Function,presence))
gtdb.genomes.s.list <- gtdb.genomes.s %>% distinct(Genome) 
gtdb.genomes.s.list <- gtdb.genomes.s.list %>% mutate(Genome=gsub('.[^.]*$', '', gtdb.genomes.s.list$Genome))

img.s.cycling <- diamond.all.nif %>% filter(Genome %in% gtdb.genomes.s.list$Genome) %>% group_by(gene,Ecosystem_Subtype_NEW,Genome) %>% summarize(Count=n()) %>% 
  left_join(img.eco.sub.nif.order, by="Ecosystem_Subtype_NEW") %>% ungroup() %>% filter(!is.na(eco_sub_order)) %>% 
  left_join(img.nif.totals, by=c("gene","Ecosystem_Subtype_NEW")) %>% mutate(perc_gene=100*Count/Total, Category="S cycling")
# img.s.cycling$gene <- factor(img.s.cycling$gene, levels=c("NifH","NifD","NifK"))
# img.s.cycling <- img.s.cycling %>% mutate(Ecosystem_Subtype_NEW=reorder(Ecosystem_Subtype_NEW,eco_sub_order))

#C1 metabolism
gtdb.genomes.c1 <- gtdb.met.pres %>% filter(Category=="C1 metabolism") %>% select(c(Genome,phylum,class,Function,presence))
gtdb.genomes.c1.list <- gtdb.genomes.c1 %>% distinct(Genome) 
gtdb.genomes.c1.list <- gtdb.genomes.c1.list %>% mutate(Genome=gsub('.[^.]*$', '', gtdb.genomes.c1.list$Genome))

img.c1.cycling <- diamond.all.nif %>% filter(Genome %in% gtdb.genomes.c1.list$Genome) %>% group_by(gene,Ecosystem_Subtype_NEW,Genome) %>% summarize(Count=n()) %>% 
  left_join(img.eco.sub.nif.order, by="Ecosystem_Subtype_NEW") %>% ungroup() %>% filter(!is.na(eco_sub_order)) %>% 
  left_join(img.nif.totals, by=c("gene","Ecosystem_Subtype_NEW")) %>% mutate(perc_gene=100*Count/Total, Category="C1 metabolism")
# img.s.cycling$gene <- factor(img.s.cycling$gene, levels=c("NifH","NifD","NifK"))
# img.s.cycling <- img.s.cycling %>% mutate(Ecosystem_Subtype_NEW=reorder(Ecosystem_Subtype_NEW,eco_sub_order))

img.metabolic <- cbind(rbind(img.s.cycling, img.c1.cycling))
img.metabolic.sums <- img.metabolic %>% group_by(gene,Ecosystem_Subtype_NEW,Category) %>% summarise(perc_gene=sum(perc_gene))
img.metabolic.sums$gene <- factor(img.metabolic.sums$gene, levels=c("NifH","NifD","NifK"))
img.metabolic.sums <- img.metabolic.sums %>% mutate(Ecosystem_Subtype_NEW=reorder(Ecosystem_Subtype_NEW,eco_sub_order))

img.by.eco.type.metabolic <- ggplot(data=img.metabolic.sums, aes(x=Ecosystem_Subtype_NEW, y=perc_gene, fill=Category)) + facet_rep_wrap(~gene) +
  # geom_bar(stat="identity", position="stack",color="#bdbdbd",fill="#bdbdbd",width=0.5) + #scale_color_manual(values=col.phy, name=NULL) + scale_fill_manual(values=col.phy, name=NULL) +
  geom_point(shape=21) + scale_y_continuous(name="% Nif subunit proteins") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), axis.title.x=element_blank(), axis.title.y=element_text(size=10))
pdf(file="IMG nif genes metabolic eco subtype.pdf",width=8, height=3)
print(img.by.eco.type.metabolic)
dev.off()