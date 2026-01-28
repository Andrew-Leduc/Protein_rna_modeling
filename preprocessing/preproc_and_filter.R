library(Seurat)
library()

Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path,set.attributes = T,whole.header = T)
  convert_mouse <- names(convert_mouse)
  parse_row<-grep("GN=",convert_mouse, fixed=T)
  split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
  gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
  prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
  prot_parse <- grep("|",prot, fixed=T)
  gene_parse <- grep(" ",gene, fixed=T)
  split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
  split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
  split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
  split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))
  
  return(convert_mouse)
}

mRNA_raw_path <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/'
protein_dat_path <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


rna_seq <- readRDS(paste0(mRNA_raw_path,'seurat_sc_trachea.rds'))

p_all_abs <- read.csv(paste0(protein_dat_path,'/04_Gene_X_SingleCell_and_annotations/sc_protein_absolute.csv'),row.names = 1)
p_all <- read.csv(paste0(protein_dat_path,'/04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv'),row.names = 1)
p_all <- as.matrix(p_all)

meta_data <- read.csv(paste0(protein_dat_path,'/04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)
mRNA_meta <- read.csv(paste0(protein_dat_path,'/04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'),row.names = 1)




meta_data_sample <- meta_data %>% filter(Cell_Type == 'Basal')

hist(meta_data_sample$prot_total_adj)

meta_data_sample <- meta_data_sample %>% filter(prot_total_adj > - 0.2 & prot_total_adj < 0.2)

hist(meta_data_sample$prot_total_adj)

p_all_filt <- p_all[rowSums(is.na(p_all_filt)==F) > 400, meta_data_sample$ID]

convert = Proc_fasta(paste0(protein_dat_path,'Mouse.fasta'))
convert = convert %>% filter(split_prot %in% rownames(p_all_abs))
convert = convert %>% filter(split_gene %in% rownames(rna_seq@assays$RNA$counts))

convert_lim <- convert %>% filter(split_prot %in% rownames(p_all_filt)) 
p_all_filt <- p_all_filt[convert_lim$split_prot, ]
rownames(p_all_filt) <- convert_lim$split_gene

write(p_all_filt,'/Users/andrewleduc/Desktop/Github/rna_prot_mod/dat/norm_log2_prot.csv')

