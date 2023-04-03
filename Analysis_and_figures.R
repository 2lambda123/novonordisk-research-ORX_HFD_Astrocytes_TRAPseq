
##### Loading libraries and functions #####

# Load libraries

# General
library(tidyverse)
library(openxlsx)
#library(data.table)

# Annotation
#library(biomaRt)
library(annotables)

# RNA-seq processing & analysis
library(tximeta)
library(wasabi)
library(sleuth)

# Single cell deconvolution
library(SeuratDisk)
library(Seurat)
library(MuSiC)

# Gene set enrichment analysis
library(fgsea)
library(msigdf)

# Plotting
library(lemon)
library(ggradar)
library(ggvenn)
library(cowplot)
library(RColorBrewer)
#require(ggiraphExtra)

set.seed(111)

# For reporting
sessionInfo()



##### Loading data #####



t2g <- left_join(annotables::grcm38_tx2gene,annotables::grcm38)[c("enstxp","ensgene","symbol")] %>%
  dplyr::rename(target_id = enstxp, 
         ens_gene = ensgene,
         ext_gene = symbol)
dim(t2g)
sum(is.na(t2g$ext_gene))
sum(t2g$ext_gene == "")
sum(t2g$ext_gene == " ")
sum(duplicated(t2g$target_id))

t2g <- t2g %>%
  filter(!duplicated(t2g$target_id)) 
dim(t2g)
sum(duplicated(t2g$target_id))

# t2g <- t2g %>%
#   filter(!is.na(ext_gene),
#          ext_gene != "",
#          ext_gene != " ") 

head(t2g)


path_to_files <- ".//Salmon_homemade_mm10_bootstraped_update/"
files <- list.files(path=path_to_files, pattern="quant.sf", full.names = T, recursive = T)
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) %>%
  gsub("_S[0-9]*_.*","",.)
#Note: Tximeta seems to work with a named vector as for Tximport, but in theory requires a dataframe. adding the below.
files <- files %>%
  as.data.frame() %>%
  magrittr::set_colnames("files") %>%
  rownames_to_column("names")
se <- tximeta(files)
gse <- summarizeToGene(se, countsFromAbundance = "no")


path_to_files <- ".//Salmon_homemade_mm10_bootstraped_update/"
files <- list.dirs(path=path_to_files, full.names = T, recursive = F)
files <- files[grepl("trimmed",files)]
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) %>%
  gsub("_S[0-9]*_.*","",.)
prepare_fish_for_sleuth(files)



#### Load and process metadata ####

# Load metadata 
sample_metadata <- read.csv(".//Single_folder_Output_mm/SampleSheet.csv",skip = 20) %>%
  as.data.frame()
head(sample_metadata)
sample_metadata <- sample_metadata %>% 
  dplyr::select(Sample_Name) %>%
  mutate(Sample_Type = gsub("_.*","",Sample_Name),
         Diet = gsub(".*-","",Sample_Name),
         Treatment = gsub(".*_|-.*","",Sample_Name),
         Sample_N = gsub("SN_|IP_","",Sample_Name),
         Sample_N = gsub("_.*","",Sample_N)) %>%
  mutate(Sample_Group = paste0(Diet,"_",Treatment))

head(sample_metadata)
dim(sample_metadata)



# Extract MultiQC stats
sample_metadata <- read.table(".//Salmon_homemade_mm10_bootstraped_update/multiqc_data/multiqc_general_stats.txt",sep="\t",header = T) %>%
  as.data.frame() %>%
  filter(!grepl("Undetermined",Sample))  %>%
  dplyr::rename(Sample_Name = Sample) %>%
  mutate(Sample_Name = gsub(".*aux_info \\| |_S[0-9].*|\\|","",Sample_Name)) %>%
  rename_with(~gsub("Salmon_mqc.generalstats.salmon.","",.x)) %>%
  left_join(sample_metadata,.) %>%
  mutate(rownames = Sample_Name) %>%
  tibble::column_to_rownames("rownames")

# Stats for paper
summary(sample_metadata$percent_mapped)
sd(sample_metadata$percent_mapped)
summary(sample_metadata$num_mapped)/1e6
sd(sample_metadata$num_mapped)/1e6






#### MUSIC: cell-type deconvolution ####


# Load the count data 
sample_counts <- SummarizedExperiment::assays(gse)[["counts"]] %>%
  as.data.frame() %>%
  rownames_to_column("ensgene") %>%
  mutate(ensgene = gsub("\\..*","",ensgene)) %>%
  left_join(annotables::grcm38[,c("ensgene","symbol")]) %>%
  dplyr::select(-ensgene) %>%
  filter(!duplicated(symbol),
         !is.na(symbol)) %>%
  column_to_rownames("symbol")  %>%
  dplyr::select(matches(sample_metadata$Sample_Name)) %>%
  distinct() %>%
  as.matrix()
dim(sample_counts)



#sanity check 
mean(colnames(sample_counts) %in% sample_metadata$Sample_Name)
all(colnames(sample_counts) == sample_metadata$Sample_Name)


# sample_counts_GEO <- sample_counts %>%
#   as.data.frame() %>%
#   rownames_to_column("symbol")
# head(sample_counts_GEO)
# write.table(sample_counts_GEO,
#             "./GEO_submission_Dec22_update/GSE214683_gene_level_raw_counts_matrix_for_deconvolution.txt",
#             row.names = F,
#             col.names = T,
#             quote = F)




# Load the single-cell RNAseq hippocampi data
# https://www.sciencedirect.com/science/article/pii/S009286741830789X?via%3Dihub
# http://mousebrain.org/loomfiles_level_L1.html

l1_hippocampus <- SeuratDisk::Connect(filename = "./MouseBrainOrg/l1_hippocampus.loom", mode = "r")
l1_hippocampus.seurat <- as.Seurat(l1_hippocampus)


# Set identities
Idents(l1_hippocampus.seurat) <- l1_hippocampus.seurat$Class

# Remove unwanted cluster
cell.use <- row.names(l1_hippocampus.seurat@meta.data)[-which(l1_hippocampus.seurat@meta.data$Class %in% "Excluded")]
l1_hippocampus.seurat_subset <- subset(l1_hippocampus.seurat, cells = cell.use)
table(l1_hippocampus.seurat_subset@meta.data$Class)


l1_hippocampus_ExpressionSet <- ExpressionSet(assayData=as.matrix(l1_hippocampus.seurat_subset@assays$RNA@counts),
                                              phenoData=new("AnnotatedDataFrame",l1_hippocampus.seurat_subset@meta.data))



TrapSeq_ExpressionSet_temp <- ExpressionSet(assayData=as.matrix(sample_counts),
                                            phenoData=new("AnnotatedDataFrame",sample_metadata))

Est.prop.TRAPseq = music_prop(bulk.eset = TrapSeq_ExpressionSet_temp, 
                              sc.eset = l1_hippocampus_ExpressionSet, 
                              clusters = 'Class',
                              samples = 'DonorID', 
                              select.ct = NULL,
                              verbose = T)

((Est.prop.TRAPseq$Est.prop.weighted + Est.prop.TRAPseq$Est.prop.allgene) / 2 ) %>%
  as.data.frame() %>%
  rowSums()



sample_metadata <- ((Est.prop.TRAPseq$Est.prop.weighted + Est.prop.TRAPseq$Est.prop.allgene) / 2 ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample_Name") %>%
  left_join(sample_metadata,.)

sample_metadata %>%
  dplyr::select(-Sample_Name,-Diet,-Treatment,-Sample_N,-Sample_Group,-percent_mapped,-num_mapped) %>%
  pivot_longer(-Sample_Type, names_to = "celltype",values_to = "percent") %>%
  group_by(Sample_Type,celltype) %>%
  summarise(mean=mean(percent)*100,
            sd = sd(percent)*100)

#### Cell proportion plots ####


radar_data <- sample_metadata %>%
  dplyr::select(-Sample_Name,-Diet,-Treatment,-Sample_N,-Sample_Group,-percent_mapped,-num_mapped) %>%
  mutate(Sample_Type = case_when(Sample_Type == "IP" ~ "Immunoprecipitate",
                                 Sample_Type == "SN" ~ "Supernatant")) %>%
  group_by(Sample_Type) %>%
  summarise_all(mean) 
G1 <- ggradar(
  radar_data, 
  group.point.size = 1, 
  axis.label.size	= 3.5,
  legend.text.size = 10,
  legend.position = "bottom",
  grid.label.size = 4,
  gridline.label.offset = 0.3,
  gridline.mid.colour = "grey"
)



col_pal <- brewer.pal(7,"Dark2")

G2 <- sample_metadata %>%
  dplyr::select(Sample_Type,Astrocytes, Neurons, Oligos, Blood, Ependymal, Immune, Vascular) %>%
  pivot_longer(-Sample_Type, names_to = "Cell_Type", values_to = "Percentage")  %>%
  mutate(Sample_Type = case_when(Sample_Type == "IP" ~ "Immunoprecipitate",
                                 Sample_Type == "SN" ~ "Supernatant")) %>%
  ggplot(aes(x=Cell_Type,y=Percentage,colour=Cell_Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  facet_wrap(~Sample_Type,scales = "free") + 
  theme_classic()  +
  xlab(element_blank()) + 
  ylab('Percentage of cells\n ') +
  geom_hline(yintercept = 0, colour="black", linetype="dotted") +
  guides(colour=guide_legend(title="Cell type")) +
  scale_color_manual(name = " ",
                     values= col_pal)  



col_pal <- brewer.pal(7,"Paired")


G3 <- sample_metadata %>%
  mutate(Sample_Type = case_when(Sample_Type == "IP" ~ "Immunoprecipitate",
                                 Sample_Type == "SN" ~ "Supernatant")) %>%
  mutate(Sample_Group = gsub("_"," ",Sample_Group)) %>%
  mutate(Sample_Group = factor(Sample_Group,levels=c("Chow Sham", "HFD Sham", "HFD GDX"), ordered = T)) %>%
  dplyr::select(Sample_Type, Sample_Group,Astrocytes, Neurons, Oligos, Blood, Ependymal, Immune, Vascular) %>%
  pivot_longer(-c(Sample_Group,Sample_Type), names_to = "Cell_Type", values_to = "Percentage")  %>%
  ggplot(aes(x=Cell_Type,y=Percentage,colour=Sample_Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width=0.3)) +
  facet_rep_wrap(Sample_Type~Cell_Type,scales="free",nrow=2) +
  theme_classic()  +
  xlab(element_blank()) + 
  ylab('Percentage of cells\n ')+ 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_color_manual(name = " ",
                     values= c(col_pal[2],col_pal[4],col_pal[6]))  







# test if any cell % is significantly different between groups


for (var in c("Astrocytes", "Neurons", "Oligos", "Blood", "Ependymal", "Immune", "Vascular")){
  print(var)
  formula = as.formula(paste0(var," ~ Sample_Group"))
  
  print("IP")
  
  sample_metadata %>%
    filter(Sample_Type == "IP") %>%
    rstatix::pairwise_t_test(formula, p.adjust.method = "bonferroni")  %>%
    print()
  
  print("SN")
  
  sample_metadata %>%
    filter(Sample_Type == "SN") %>%
    rstatix::pairwise_t_test(formula, p.adjust.method = "bonferroni") %>%
    print()
  
}








#### Sleuth: prepare data ####


sleuth_data <- files %>%
  as.data.frame() %>%
  magrittr::set_colnames("path") %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_metadata,.) %>%
  mutate(sample=Sample_Name) %>%
  mutate(DE_group = ifelse(Sample_Group == "HFD_GDX", "HFD_GDX","BothDiet_Sham")) %>%
  dplyr::select(sample, Sample_Type, DE_group, Diet, Sample_Group, Astrocytes, Neurons, Oligos, Blood, Ependymal, Immune, Vascular, path)



# Save TPM matrix for GEO


# sleuth_for_tpm <- sleuth_data %>%
#   sleuth_prep(extra_bootstrap_summary = TRUE,
#               read_bootstrap_tpm = TRUE,
#               aggregation_column = 'ext_gene',
#               target_mapping = t2g,
#               transformation_function = function(x) log2(x + 0.5),
#               gene_mode = T)
# 
# 
# tpm_matrix <- sleuth_to_matrix(sleuth_for_tpm, 'obs_norm', 'tpm')
# 
# 
# tpm_matrix <- tpm_matrix %>%
#   as.data.frame() %>%
#   rownames_to_column("symbol")
# 
# write.table(tpm_matrix,
#             "./GEO_submission_Dec22_update/GSE214683_gene_level_TPM_matrix.txt",
#             row.names = F,
#             col.names = T,
#             quote = F)
# 
# 
# 





#### Sleuth: DE analysis ####





list_results <- list()
list_sleuth <- list()

for (subset_name in c("IP","SN")){
  
  print(subset_name)
  
  list_sleuth[[paste0(subset_name,"_GDXeffect")]] <- sleuth_data %>%
    filter(Sample_Type == subset_name) %>%
    filter(Sample_Group %in% c("HFD_Sham","HFD_GDX")) %>%
    mutate(Pairwise_GDXeffect = factor(Sample_Group,levels=c("HFD_Sham","HFD_GDX"))) %>%
    sleuth_prep(extra_bootstrap_summary = TRUE, 
                read_bootstrap_tpm = TRUE,
                aggregation_column = 'ext_gene', 
                target_mapping = t2g,
                transformation_function = function(x) log2(x + 0.5),
                gene_mode = T) #https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html
  list_sleuth[[paste0(subset_name,"_GDXeffect")]] <- sleuth_fit(list_sleuth[[paste0(subset_name,"_GDXeffect")]], ~ Pairwise_GDXeffect  + Astrocytes + Neurons + Oligos + Blood + Ependymal + Immune , 'full') # removed "vascular" to have a valid model
  list_sleuth[[paste0(subset_name,"_GDXeffect")]] <- sleuth_fit(list_sleuth[[paste0(subset_name,"_GDXeffect")]], ~ Astrocytes  + Neurons + Oligos + Blood + Ependymal + Immune , 'reduced')
  list_sleuth[[paste0(subset_name,"_GDXeffect")]] <- sleuth_lrt(list_sleuth[[paste0(subset_name,"_GDXeffect")]], 'reduced', 'full')
  list_sleuth[[paste0(subset_name,"_GDXeffect")]] <- sleuth_wt(list_sleuth[[paste0(subset_name,"_GDXeffect")]], 'Pairwise_GDXeffectHFD_GDX')
  
  list_results[[paste0(subset_name,"_GDXeffect")]] <- sleuth_results(list_sleuth[[paste0(subset_name,"_GDXeffect")]], 'reduced:full', 'lrt', show_all = FALSE)
  list_results[[paste0(subset_name,"_GDXeffect")]] <- sleuth_results(list_sleuth[[paste0(subset_name,"_GDXeffect")]], 'Pairwise_GDXeffectHFD_GDX', 'wt', show_all = FALSE) %>% 
    dplyr::select(target_id, b) %>% 
    left_join(list_results[[paste0(subset_name,"_GDXeffect")]],.) %>%
    drop_na()
  
  
  
  
  
  list_sleuth[[paste0(subset_name,"_forLinear")]] <- sleuth_data %>%
    filter(Sample_Type == subset_name) %>%
    mutate(Ordered_grouping = factor(Sample_Group,levels=c("Chow_Sham", "HFD_Sham", "HFD_GDX"), ordered = T)) %>%
    sleuth_prep(extra_bootstrap_summary = TRUE, 
                read_bootstrap_tpm = TRUE,
                aggregation_column = 'ext_gene', 
                target_mapping = t2g,
                transformation_function = function(x) log2(x + 0.5),
                gene_mode = T)
  
  list_sleuth[[paste0(subset_name,"_Ordinal")]] <- list_sleuth[[paste0(subset_name,"_forLinear")]]
  list_sleuth[[paste0(subset_name,"_Ordinal")]] <- sleuth_fit(list_sleuth[[paste0(subset_name,"_Ordinal")]], ~ Ordered_grouping  + Astrocytes + Neurons + Oligos + Blood + Ependymal + Immune , 'full') # removed "vascular" to have a valid model
  list_sleuth[[paste0(subset_name,"_Ordinal")]] <- sleuth_fit(list_sleuth[[paste0(subset_name,"_Ordinal")]], ~ Astrocytes  + Neurons + Oligos + Blood + Ependymal + Immune , 'reduced')
  list_sleuth[[paste0(subset_name,"_Ordinal")]] <- sleuth_lrt(list_sleuth[[paste0(subset_name,"_Ordinal")]], 'reduced', 'full')
  list_sleuth[[paste0(subset_name,"_Ordinal")]] <- sleuth_wt(list_sleuth[[paste0(subset_name,"_Ordinal")]], 'Ordered_grouping.L')
  
  list_results[[paste0(subset_name,"_Ordinal")]] <- sleuth_results(list_sleuth[[paste0(subset_name,"_Ordinal")]], 'reduced:full', 'lrt', show_all = FALSE)
  list_results[[paste0(subset_name,"_Ordinal")]] <- sleuth_results(list_sleuth[[paste0(subset_name,"_Ordinal")]], 'Ordered_grouping.L', 'wt', show_all = FALSE) %>% 
    dplyr::select(target_id, b) %>% 
    left_join(list_results[[paste0(subset_name,"_Ordinal")]],.) %>%
    drop_na()
  
  
}


for (subset_name2 in names(list_results)){
  
  list_results[[subset_name2]]$Fraction <- gsub("_.*","",subset_name2)
  list_results[[subset_name2]]$Analysis <- gsub(".*_","",subset_name2)
  list_results[[subset_name2]]$Subset <- subset_name2
  
}


results_all <- do.call(rbind,list_results) %>%
  as.data.frame() %>% 
  dplyr::select(Subset,Fraction, Analysis, ens_gene, target_id, pval, qval, test_stat, b) %>%
  magrittr::set_colnames(c("subset","fraction","analysis","ensgene","symbol","p_val","p_adj","stat","log2_fc")) %>%
  mutate(stat = stat * sign(log2_fc))


table(results_all$subset, results_all$p_adj < 0.05)

results_all %>%
  dplyr::select(subset,symbol,p_adj,log2_fc) %>%
  distinct() %>%
  group_by(subset) %>%
  tally(p_adj <= 0.05)

results_all %>%
  dplyr::select(subset,symbol,p_adj,log2_fc) %>%
  distinct() %>%
  mutate(sign = sign(log2_fc)) %>%
  group_by(subset,sign) %>%
  tally(p_adj <= 0.05)


results_all %>%
  filter(subset == "IP_Ordinal" & p_adj <= 0.05) %>%
  filter(log2_fc > 0) %>%
  arrange(p_adj)


results_all %>%
  filter(subset == "IP_Ordinal" & p_adj <= 0.05) %>%
  filter(log2_fc < 0) %>%
  arrange(p_adj)

results_all %>%
  filter(subset == "SN_Ordinal" & p_adj <= 0.05) %>%
  filter(log2_fc > 0) %>%
  arrange(p_adj)


results_all %>%
  filter(subset == "SN_Ordinal" & p_adj <= 0.05) %>%
  filter(log2_fc < 0) %>%
  arrange(p_adj)



#### DE plots - Venn diagram ####


list_sig_genes <- list(
  IP_Pairwise = results_all %>% filter(subset == "IP_GDXeffect" & p_adj <= 0.05) %>% pull(symbol) %>% unique(), 
  IP_Ordinal = results_all %>% filter(subset == "IP_Ordinal" & p_adj <= 0.05) %>% pull(symbol) %>% unique(), 
  SN_Ordinal = results_all %>% filter(subset == "SN_Ordinal" & p_adj <= 0.05) %>% pull(symbol) %>% unique(),
  SN_Pairwise = results_all %>% filter(subset == "SN_GDXeffect" & p_adj <= 0.05) %>% pull(symbol) %>% unique()
)
names(list_sig_genes) <- gsub("_"," ",names(list_sig_genes))


list_sig_genes_down <- list(
  IP_Pairwise_down = results_all %>% filter(subset == "IP_GDXeffect" & p_adj <= 0.05 & log2_fc < 0) %>% pull(symbol) %>% unique(), 
  IP_Ordinal_down = results_all %>% filter(subset == "IP_Ordinal" & p_adj <= 0.05 & log2_fc < 0) %>% pull(symbol) %>% unique(), 
  SN_Ordinal_down = results_all %>% filter(subset == "SN_Ordinal" & p_adj <= 0.05 & log2_fc < 0) %>% pull(symbol) %>% unique(),
  SN_Pairwise_down= results_all %>% filter(subset == "SN_GDXeffect" & p_adj <= 0.05 & log2_fc < 0) %>% pull(symbol) %>% unique()
)
names(list_sig_genes_down) <- gsub("_","\n",names(list_sig_genes_down))


list_sig_genes_up <- list(
  IP_Pairwise_up = results_all %>% filter(subset == "IP_GDXeffect" & p_adj <= 0.05 & log2_fc > 0) %>% pull(symbol) %>% unique(), 
  IP_Ordinal_up = results_all %>% filter(subset == "IP_Ordinal" & p_adj <= 0.05 & log2_fc > 0) %>% pull(symbol) %>% unique(), 
  SN_Ordinal_up = results_all %>% filter(subset == "SN_Ordinal" & p_adj <= 0.05 & log2_fc > 0) %>% pull(symbol) %>% unique(),
  SN_Pairwise_up = results_all %>% filter(subset == "SN_GDXeffect" & p_adj <= 0.05 & log2_fc > 0) %>% pull(symbol) %>% unique()
)
names(list_sig_genes_up) <- gsub("_"," ",names(list_sig_genes_up))




col_pal <- brewer.pal(8,"Set1")
ggven_colors <- c(col_pal[2],col_pal[5],col_pal[3],col_pal[1])



ggven_all <- ggvenn(list_sig_genes, 
                    fill_color = ggven_colors,
                    show_percentage = FALSE,
                    stroke_size = 0.3, 
                    fill_alpha = 0.3,
                    set_name_size = 3,
                    text_size = 3)  +
  scale_fill_manual(breaks=c("A","B","C","D"),
                    values= ggven_colors,
                    labels = names(list_sig_genes)) + 
  guides(fill = guide_legend(title = "", override.aes = list(colour = "black"))) +
  theme(legend.position = "right") + 
  ggtitle("All significant genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


ggven_down <- ggvenn(list_sig_genes_down, 
                     fill_color = ggven_colors,
                     show_percentage = FALSE,
                     stroke_size = 0.3, 
                     fill_alpha = 0.3,
                     set_name_size = 3,
                     text_size = 3)  +
  scale_fill_manual(breaks=c("A","B","C","D"),
                    values= ggven_colors,
                    labels = names(list_sig_genes))  + 
  ggtitle("Significantly downregulated genes")  +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


ggvenn_up <- ggvenn(list_sig_genes_up, 
                    fill_color = ggven_colors,
                    show_percentage = FALSE,
                    stroke_size = 0.3, 
                    fill_alpha = 0.3,
                    set_name_size = 3,
                    text_size = 3)  +
  scale_fill_manual(breaks=c("A","B","C","D"),
                    values= ggven_colors,
                    labels = names(list_sig_genes))  + 
  ggtitle("Significantly upregulated genes",) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))  + 
  theme(legend.text = element_text(size=12))



ggven_all$layers[[3]] <- NULL
ggven_down$layers[[3]] <- NULL
ggvenn_up$layers[[3]] <- NULL

legend_ggven <- get_legend(ggven_all + theme(legend.text = element_text(size=11)))


y_limits <- c(-1.8, 1.2)
ggven_all <- ggven_all + guides(fill = "none") + scale_y_continuous(limits = y_limits) + coord_fixed(ratio = 1.2)
ggven_down <- ggven_down + scale_y_continuous(limits = y_limits) + coord_fixed(ratio = 1.2)
ggvenn_up <- ggvenn_up + scale_y_continuous(limits = y_limits) + coord_fixed(ratio = 1.2)




#### DE plots - Volcanot plots ####


data_for_plot <- results_all %>%
  filter(analysis == "GDXeffect") %>%
  dplyr::select(-ensgene) %>%
  distinct() %>%
  mutate(color_group = case_when(
    symbol %in% list_sig_genes[["IP Ordinal"]] & fraction == "IP" & log2_fc < 0 ~ "Sig_Ordinal_Neg",
    symbol %in% list_sig_genes[["IP Ordinal"]] & fraction == "IP" & log2_fc > 0 ~ "Sig_Ordinal_Pos",
    symbol %in% list_sig_genes[["SN Ordinal"]] & fraction == "SN" & log2_fc < 0 ~ "Sig_Ordinal_Neg",
    symbol %in% list_sig_genes[["SN Ordinal"]] & fraction == "SN" & log2_fc > 0~ "Sig_Ordinal_Pos",
    p_adj <= 0.05 & log2_fc < 0 ~ "Sig_Neg",
    p_adj <= 0.05 & log2_fc > 0 ~ "Sig_Pos"
  ) ) %>%
  mutate(fraction = gsub("IP", "Immunoprecipitate", fraction)) %>%
  mutate(fraction = gsub("SN", "Supernatant", fraction))

table(data_for_plot$color_group)

subset_for_text <- data_for_plot %>%
  group_by(fraction) %>%
  slice_min(order_by=p_adj, n=10, with_ties=F)

subset_for_text <- data_for_plot %>%
  filter(grepl("Cont|Ord",color_group)) %>%
  rbind(subset_for_text) %>%
  distinct()


col_pal <- brewer.pal(8,"Paired")

volcano_plot <- ggplot(data_for_plot, aes(x=log2_fc,y=-log10(p_adj))) + 
  geom_point(alpha = 0.5,
             aes(color = color_group)) + 
  theme_classic() + 
  labs(x="Log2 Fold Change", y="-log10(p_adj)") + 
  geom_point(data=subset(subset_for_text,grepl("Cont|Ord",color_group)),
             aes(x=log2_fc,y=-log10(p_adj),color = color_group)) +
  ggrepel::geom_label_repel(data = subset_for_text, 
                            aes(label=symbol,
                                color = color_group),
                            label.size = NA,  
                            label.padding=.1, 
                            na.rm=TRUE,
                            fill = alpha(c("white"),0.8),
                            force= 8, 
                            max.overlaps = 40,
                            show.legend  = FALSE) +
  facet_wrap(~ fraction, ncol = 2) +
  scale_color_manual(name = " ",
                     breaks=c("Sig_Ordinal_Neg","Sig_Ordinal_Pos","Sig_Neg", "Sig_Pos"),
                     values= c(col_pal[8],col_pal[4],col_pal[2],col_pal[6]),
                     labels = c("Sig. Downregulated\n(Pairwise + Ordinal)",
                                "Sig. Upregulated\n(Pairwise + Ordinal)",
                                "Sig. Downregulated\n(Pairwise only)",
                                "Sig. Upregulated\n(Pairwise only)"))  + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=9))


volcano_plot






#### GSEA ####




#gene sets
c2cp <- msigdf.mouse %>%  
  filter(category_code=="c2") %>% 
  filter(grepl("cp.",category_subcode)) %>% 
  dplyr::select(geneset, mouse.symbol) %>% 
  group_by(geneset) %>% 
  summarize(mouse.symbol=list(mouse.symbol)) %>% 
  deframe()
c5bp <- msigdf.mouse %>%  
  filter(category_subcode=="go.bp") %>% 
  dplyr::select(geneset, mouse.symbol) %>% 
  group_by(geneset) %>% 
  summarize(mouse.symbol=list(mouse.symbol)) %>% 
  deframe()
c5mf <- msigdf.mouse %>%  
  filter(category_subcode=="go.mf") %>% 
  dplyr::select(geneset, mouse.symbol) %>% 
  group_by(geneset) %>% 
  summarize(mouse.symbol=list(mouse.symbol)) %>% 
  deframe()
hallmark <- msigdf.mouse %>%  
  filter(category_code=="h") %>% 
  dplyr::select(geneset, mouse.symbol) %>% 
  group_by(geneset) %>% 
  summarize(mouse.symbol=list(mouse.symbol)) %>% 
  deframe()


gene_set_list <- list(hallmark,c2cp,c5bp,c5mf)
names(gene_set_list) <- c("hallmark","c2cp","c5bp","c5mf")

gsea_list <- list()
gsea_list_collapsed <- list()

N_cores <- 10


for (subset_name in c("IP_GDXeffect","IP_Ordinal","SN_GDXeffect","SN_Ordinal")){
  
  print(subset_name)
  
  
  
  temp_data <- results_all %>%
    filter(subset == subset_name) %>%
    arrange(stat) %>%
    dplyr::select(symbol,stat,p_val) %>%
    drop_na() %>%
    distinct()
  
  rnk_list <- setNames(object = temp_data$stat, nm = temp_data$symbol)

  
  for (geneset_name in names(gene_set_list)){
    
    
    subset_name2 <- paste0(subset_name,"_",geneset_name)
    
    print(geneset_name)
    gsea_list[[subset_name2]] <- fgsea(pathways = gene_set_list[[geneset_name]],
                                       stats = rnk_list,
                                       minSize = 20,
                                       eps = 1e-30,
                                       nproc = N_cores) 
    
    
    
    gsea_list[[subset_name2]] %>%
      filter(padj <= 0.05) %>%
      as.data.frame() %>%
      dplyr::select(pathway,NES,padj,size) %>%
      arrange(-abs(NES)) %>%
      head(40) %>%
      print()
    
    
    gsea_list_collapsed[[subset_name2]] <- collapsePathways(fgseaRes = filter(gsea_list[[subset_name2]], padj <= 0.05),
                                                            pathways = gene_set_list[[geneset_name]],
                                                            stats = rnk_list) #,
    #nperm = 1e6)
    
    
    gsea_list_collapsed[[subset_name2]] <- gsea_list[[subset_name2]] %>%
      filter(pathway %in% gsea_list_collapsed[[subset_name2]]$mainPathways)
    
    
    gsea_list_collapsed[[subset_name2]] %>%
      filter(padj <= 0.05) %>%
      as.data.frame() %>%
      dplyr::select(pathway,NES,padj,size) %>%
      arrange(-abs(NES)) %>%
      head(40) %>%
      print()
    
    
    
    
    
    gsea_list[[subset_name2]]$Fraction <- gsub("_.*","",subset_name2)
    gsea_list[[subset_name2]]$Analysis <-  gsub("IP|SN|_|c2cp|c5bp|c5mf|hallmark|c8","",subset_name2)
    gsea_list[[subset_name2]]$GeneSet <- gsub(".*_","",subset_name2)
    gsea_list[[subset_name2]]$Subset <- subset_name2
    
    gsea_list_collapsed[[subset_name2]]$Fraction <- gsub("_.*","",subset_name2)
    gsea_list_collapsed[[subset_name2]]$Analysis <-  gsub("IP|SN|_|c2cp|c5bp|c5mf|hallmark|c8","",subset_name2)
    gsea_list_collapsed[[subset_name2]]$Geneset <- gsub(".*_","",subset_name2)
    gsea_list_collapsed[[subset_name2]]$Subset <- subset_name2
    
    
    
  }
  
}



gsea_list_collapsed <- do.call(rbind,gsea_list_collapsed) %>%
  as.data.frame() 

head(gsea_list_collapsed)




#### GSEA plots ####


list_GSEA_plot <- list()




set_name <- "hallmark"
N_sig = 10
col_pal <- brewer.pal(8,"Paired")

top_min <- gsea_list_collapsed %>%
  filter(Analysis %in% c("GDXeffect","Ordinal")) %>%
  filter(Geneset == set_name) %>%
  group_by(Analysis,Fraction) %>%
  slice_min(order_by = NES, n=N_sig) 
top_max <- gsea_list_collapsed %>%
  filter(Analysis %in% c("GDXeffect","Ordinal")) %>%
  filter(Geneset == set_name) %>%
  group_by(Analysis,Fraction) %>%
  slice_max(order_by = NES, n=N_sig)
top_all <- rbind(top_min,top_max) %>%
  pull(pathway)

list_GSEA_plot[["hallmark"]] <- gsea_list_collapsed %>%
  filter(Analysis %in% c("GDXeffect","Ordinal"))  %>%
  filter(pathway %in% top_all) %>%
  dplyr::select(pathway,padj,NES,Fraction,Analysis) %>%
  mutate(pathway = gsub("_"," ",pathway)) %>%
  mutate(pathway = gsub("HALLMARK ","",pathway)) %>%
  mutate(pathway = stringr::str_wrap(pathway, 40)) %>%
  filter(padj <= 0.05) %>%
  mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>%
  mutate(group_name = paste0(Fraction," ",Analysis)) %>%
  mutate(group_name = factor(group_name, levels = c("SN Ordinal","SN Pairwise","IP Ordinal","IP Pairwise"))) %>%
  mutate(pathway = reorder(pathway, NES, mean)) %>%
  ggplot(aes(x = pathway, y=NES, fill=group_name)) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  coord_flip()  + 
  theme_classic() +
  theme(panel.grid = element_blank()) + # remove grid lines
  geom_vline(xintercept = seq(0.5, length(top_all), by = 1), color="gray", size=.5, alpha=.5) + # set vertical lines between x groups 
  scale_fill_manual(name = "",
                    breaks = c("SN Ordinal","SN Pairwise","IP Ordinal","IP Pairwise"),
                    values= c(col_pal[3],col_pal[4],col_pal[1],col_pal[2]),
                    guide = guide_legend(reverse = TRUE, override.aes = list(colour = "black"))) +
  xlab("Hallmark pathway") +
  ylab("Normalized Enrichment Score (NES)") 


list_GSEA_plot[["hallmark"]] 




N_sig =5


for (set_name in  c("c2cp","c5bp","c5mf")){
  
  
  top_min <- gsea_list_collapsed %>%
    filter(Analysis %in% c("GDXeffect","Ordinal")) %>%
    filter(Geneset == set_name) %>%
    group_by(Analysis,Fraction) %>%
    slice_min(order_by = NES, n=N_sig) 
  top_max <- gsea_list_collapsed %>%
    filter(Analysis %in% c("GDXeffect","Ordinal")) %>%
    filter(Geneset == set_name) %>%
    group_by(Analysis,Fraction) %>%
    slice_max(order_by = NES, n=N_sig)
  top_all <- rbind(top_min,top_max) %>%
    pull(pathway)
  
  list_GSEA_plot[[set_name]] <- gsea_list_collapsed %>%
    filter(Analysis %in% c("GDXeffect","Ordinal"))  %>%
    filter(pathway %in% top_all) %>%
    dplyr::select(pathway,padj,NES,Fraction,Analysis) %>%
    mutate(pathway = gsub("_"," ",pathway)) %>%
    mutate(pathway = gsub("GOBP |GOMF ","",pathway)) %>%
    # mutate(pathway = stringr::str_wrap(pathway, 60)) %>%
    filter(padj <= 0.05) %>%
    mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>%
    mutate(group_name = paste0(Fraction," ",Analysis)) %>%
    mutate(group_name = factor(group_name, levels = c("SN Ordinal","SN Pairwise","IP Ordinal","IP Pairwise"))) %>%
    mutate(pathway = reorder(pathway, NES, mean)) %>%
    ggplot(aes(x = pathway, y=NES, fill=group_name)) +
    geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
    coord_flip()  + 
    theme_classic() +
    theme(panel.grid = element_blank()) + # remove grid lines
    geom_vline(xintercept = seq(0.5, length(top_all), by = 1), color="gray", size=.5, alpha=.5) + # set vertical lines between x groups 
    scale_fill_manual(name = "",
                      breaks = c("SN Ordinal","SN Pairwise","IP Ordinal","IP Pairwise"),
                      values= c(col_pal[3],col_pal[4],col_pal[1],col_pal[2]),
                      guide = guide_legend(reverse = TRUE, override.aes = list(colour = "black"))) +
    ylab("Normalized Enrichment Score (NES)") 
  
}





#### Paper figures ####

# Figure 1
grid_volc <- plot_grid(NULL,
                       volcano_plot + guides(colour = guide_legend(nrow = 2)) + theme(legend.text = element_text(size=10)), 
                       ncol= 1,
                       rel_heights = c(0.05,1))
grid_row1 <- plot_grid(G1 + coord_cartesian(clip = "off") + guides(colour = guide_legend(nrow = 2))  + theme(legend.text = element_text(size=10)),
                       grid_volc, 
                       nrow = 1,
                       labels = c("A","B"),
                       rel_widths = c(0.4,0.6))
grid_row2 <- plot_grid(ggven_all ,
                       ggven_down,
                       ggvenn_up, 
                       legend_ggven,
                       rel_widths = c(1,1,1,0.5),
                       nrow = 1)
grid_row3 <- plot_grid(list_GSEA_plot[["hallmark"]] + guides(fill = "none") ,
                       get_legend(list_GSEA_plot[["hallmark"]] + theme(legend.text = element_text(size=11))),
                       rel_widths = c(3,0.5),
                       nrow = 1)
Figure1 <- plot_grid(grid_row1,
                     grid_row2,
                     NULL,
                     grid_row3,
                     ncol=1,
                     labels = c("","C","D",""),
                     rel_heights = c(1,0.8,0.1,1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

Figure1

ggsave(filename="Figure1.png", 
       plot=Figure1, 
       device="png", 
       path="./",
       dpi=1000)



# Figure S1
FigureS1 <- cowplot::plot_grid(NULL,
                               G2 + theme(legend.text = element_text(size=11), 
                                          axis.text.x = element_text(size=11), 
                                          axis.text.y = element_text(size=11),
                                          axis.title.y = element_text(size=12),
                                          strip.text = element_text(size=12)),
                               NULL,
                               G3 + theme(legend.text = element_text(size=11),
                                          axis.title.y = element_text(size=12)),
                               ncol=1,
                               rel_heights = c(0.1,0.8,0.1,1),
                               labels = c("A","","B","")) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

FigureS1

ggsave(filename="FigureS1.png", 
       plot=FigureS1, 
       device="png", 
       path="./",
       dpi=1000)



# Figure S22

list_GSEA_plot[["c2cp"]] <- list_GSEA_plot[["c2cp"]] + xlab("Canonical Pathway") + theme(axis.title.y = element_text(size=10))
list_GSEA_plot[["c5bp"]] <- list_GSEA_plot[["c5bp"]] + xlab("Gene Ontology pathway (Biologial Process)") + theme(axis.title.y = element_text(size=10))
list_GSEA_plot[["c5mf"]] <- list_GSEA_plot[["c5mf"]] + xlab("Gene Ontology pathway (Molecular Function)") + theme(axis.title.y = element_text(size=10))

FigureS2 <- plot_grid(NULL,
                      list_GSEA_plot[["c2cp"]],
                      NULL,
                      list_GSEA_plot[["c5bp"]],
                      NULL,
                      list_GSEA_plot[["c5mf"]],
                      ncol=1,
                      rel_heights = c(0.1,1,0.1,1,0.1,1),
                      labels = c("A","","B","","C",""),
                      align = "v") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

FigureS2

ggsave(filename="FigureS2.png", 
       plot=FigureS2, 
       device="png", 
       path="./",
       dpi=1000)

#### Supplementary data ####




gsea_all <- do.call(rbind,gsea_list) %>%
  as.data.frame() %>%
  arrange(-abs(NES))


list_results_xslx <- list(
  IP_Pairwise_DE = results_all %>% mutate(analysis = gsub("GDXeffect","Pairwise",analysis)) %>% filter(subset == "IP_GDXeffect") %>% arrange(p_adj) %>% dplyr::select(-ensgene,-subset) %>% distinct() ,
  IP_Ordinal_DE = results_all %>% filter(subset == "IP_Ordinal") %>% arrange(p_adj) %>% dplyr::select(-ensgene,-subset) %>% distinct() ,
  SN_Pairwise_DE = results_all %>% mutate(analysis = gsub("GDXeffect","Pairwise",analysis)) %>% filter(subset == "SN_GDXeffect") %>% arrange(p_adj) %>% dplyr::select(-ensgene,-subset) %>% distinct() ,
  SN_Ordinal_DE = results_all %>% filter(subset == "SN_Ordinal") %>% arrange(p_adj) %>% dplyr::select(-ensgene,-subset) %>% distinct() ,
  
  IP_Pairwise_GSEA_hallmark = gsea_all  %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "IP_GDXeffect_hallmark") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  IP_Ordinal_GSEA_hallmark = gsea_all  %>% filter(Subset == "IP_Ordinal_hallmark") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Pairwise_GSEA_hallmark = gsea_all  %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "SN_GDXeffect_hallmark") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Ordinal_GSEA_hallmark = gsea_all  %>% filter(Subset == "SN_Ordinal_hallmark") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  
  IP_Pairwise_GSEA_c2cp = gsea_all  %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "IP_GDXeffect_c2cp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  IP_Ordinal_GSEA_c2cp = gsea_all  %>% filter(Subset == "IP_Ordinal_c2cp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Pairwise_GSEA_c2cp = gsea_all  %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "SN_GDXeffect_c2cp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Ordinal_GSEA_c2cp = gsea_all %>% filter(Subset == "SN_Ordinal_c2cp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  
  IP_Pairwise_GSEA_c5bp = gsea_all %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "IP_GDXeffect_c5bp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  IP_Ordinal_GSEA_c5bp = gsea_all %>% filter(Subset == "IP_Ordinal_c5bp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Pairwise_GSEA_c5bp = gsea_all  %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "SN_GDXeffect_c5bp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Ordinal_GSEA_c5bp = gsea_all  %>% filter(Subset == "SN_Ordinal_c5bp") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  
  IP_Pairwise_GSEA_c5mf = gsea_all %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "IP_GDXeffect_c5mf") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  IP_Ordinal_GSEA_c5mf = gsea_all %>% filter(Subset == "IP_Ordinal_c5mf") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Pairwise_GSEA_c5mf = gsea_all %>% mutate(Analysis = gsub("GDXeffect","Pairwise",Analysis)) %>% filter(Subset == "SN_GDXeffect_c5mf") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() ,
  SN_Ordinal_GSEA_c5mf = gsea_all  %>% filter(Subset == "SN_Ordinal_c5mf") %>% dplyr::relocate(Fraction, Analysis, GeneSet, .before=1) %>% dplyr::select(-Subset) %>% arrange(-abs(NES)) %>% distinct() 
  
)



description_sheet <-  cbind(names(list_results_xslx),
                            c("Differential expression results for the pairwise comparison `HFD Sham vs HFD ORX` (IP fraction)",
                              "Differential expression results for the pairwise comparison `HFD Sham vs HFD ORX` (SN fraction)",
                              "Differential expression results for the ordinal regression analysis (IP fraction)",
                              "Differential expression results for the ordinal regression analysis (SN fraction)",
                              
                              "Gene-set enrichment analysis results (MisgDB hallmark pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB hallmark pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (SN fraction)",
                              "Gene-set enrichment analysis results (MisgDB hallmark pathways) for the ordinal regression analysis (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB hallmark pathways) for the ordinal regression analysis (SN fraction)",
                              
                              "Gene-set enrichment analysis results (MisgDB canonical pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB canonical pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (SN fraction)",
                              "Gene-set enrichment analysis results (MisgDB canonical pathways) for the ordinal regression analysis (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB canonical pathways) for the ordinal regression analysis (SN fraction)",
                              
                              "Gene-set enrichment analysis results (MisgDB GO biological process pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO biological process pathways) for the pairwise comparison `HFD Sham vs HFD ORX` (SN fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO biological process pathways) for the ordinal regression analysis (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO biological process pathways) for the ordinal regression analysis (SN fraction)",
                              
                              "Gene-set enrichment analysis results (MisgDB GO molecular function pathways) for the pairwise comparison `HFD Sham vs HFD GDX` (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO molecular function pathways) for the pairwise comparison `HFD Sham vs HFD GDX` (SN fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO molecular function pathways) for the ordinal regression analysis (IP fraction)",
                              "Gene-set enrichment analysis results (MisgDB GO molecular function pathways) for the ordinal regression analysis (SN fraction)"
                            )) %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("Sheetname","Description"))

list_results_xslx <- c(list(description_sheet),list_results_xslx)
names(list_results_xslx)[[1]] <- "Description"

write.xlsx(list_results_xslx,"./Supplementary_data.xlsx")



#### END ####
#### END ####
#### END ####
#### END ####
#### END ####



