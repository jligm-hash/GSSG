library(data.table)
library(R.utils)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

#! set the pwd of the input set
option_list <- list(
  make_option("--genescore_dir", type="character", default="sclinker_genescores/Alzheimers/", 
              help="Directory where you have sclinker output files"),
  make_option("--bed_dir", type="character", default="sclinker_beds/Alzheimers/", 
              help="Name of the file prefix"),
  make_option("--annot_name", type="character", default="Disease_Endothelial_L2", 
              help="Output directory with gene program files"),
  make_option("--enhancer", type="character", default="BLD", 
              help="Tissue for which enhancer-gene link is considered")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)


#! input is the score of the genes
genescore_dir <- opt$genescore_dir
bed_dir <- opt$bed_dir
annot_name <- opt$annot_name
enhancer_tissue <- opt$enhancer #! e.g. use the "BLD" tissue for the enhancer-gene link

score_file = paste0(genescore_dir, "/", annot_name, ".txt")
gene_scores = read.delim(score_file, header=F)

if(!dir.exists(paste0(bed_dir, "/", annot_name))){
  dir.create(paste0(bed_dir, "/", annot_name))
}


#! define the function
source("all_bedgraph_methods.R")


gene_scores = read.delim(score_file, header=F)

#! read the score of the genes
scores = gene_scores[,2]
names(scores) = gene_scores[,1] #! name as the gene name


out = ABC_Road_bedgraph_calc(scores,
                             output_cell = paste0(bed_dir, "/", annot_name), #! output pwd
                             tissuename = enhancer_tissue,
                             output_bed = paste0("ABC_Road_", enhancer_tissue, ".bed")) #! output file for the pwd


#! gene coordinate with 100k window
df = read.table("../../processed_data/Gene_100kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

#! add the gene score to the gene with the 100k window
temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", annot_name, "/", "100kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)





