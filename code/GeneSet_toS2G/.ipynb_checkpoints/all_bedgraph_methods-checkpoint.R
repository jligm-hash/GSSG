

#! define the function
PCHiC_bedgraph_calc <- function(scores,
                                output_cell,
                                output_bed = "temp.bed"){
  library(data.table)
  pchic_data = data.frame(fread("../../processed_data/PCHiC_peak_matrix_cutoff5.tsv", header=T))

  names1 = lapply(pchic_data$baitName, function(x) return(strsplit(as.character(x), ";")[[1]]))
  num_genes_per_connect = unlist(lapply(names1, function(x) return(length(x))))

  genes_all_connect = unlist(names1)
  indices_per_connect = unlist(lapply(1:nrow(pchic_data), function(x) return(rep(x, num_genes_per_connect[x]))))

  common_genes = intersect(unique(genes_all_connect), names(scores))
  idx = which(!is.na(match(genes_all_connect, common_genes)))

  new_grades = scores[match(genes_all_connect[idx], names(scores))]
  df3 = cbind(pchic_data[indices_per_connect[idx], c("oeChr", "oeStart", "oeEnd")], new_grades)
  df3[,1] = paste0("chr", df3[,1])

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

}


ABC_intergenic_bedgraph_calc <- function(scores,
                                         output_cell,
                                         output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("../../processed_data/",
                                   "ABCpaper_NasserFulcoEngreitz2020_Blood_AvgHiC.txt.gz")))
  df = df_pre[unique(c(grep("intergenic", df_pre$name))),]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  common_genes = intersect(names(scores), df2[,4])
  idx = which(!is.na(match(df2[,4], common_genes)))
  df3 = df2[idx,]
  grades = scores[match(df3[,4], names(scores))]
  temp = cbind(df3, grades)
  final_bed = temp[, c(1:3, 5)]

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ABC_Road_bedgraph_calc <- function(scores, #! gene score calculated by scRNA
                                      output_cell, #! output cell name
                                      tissuename,  ## "BLD" / "BRN"
                                      output_bed = "temp.bed"){ #! output bed pwd

  df_pre = data.frame(fread(paste0("../../processed_data/",
                                   "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))

    
  df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
    
  if(tissuename != "ALL"){

    if(tissuename == "BRN"){
      tissuenames2 = c("neur", "Neur", "astro", "spinal", "Brain", "brain")
    }
    if(tissuename == "PANC"){
      tissuenames2 = c("pancrea")
    }
    if(tissuename == "BLD"){
      tissuenames2 = as.character(read.table("../../processed_data/ABC.listbloodQC.txt", header=F)[,1])
        #! get the tissue name2
        #! 69 tissue names
        
    }
    if(tissuename == "LNG"){
      tissuenames2 = c("lung")
    }
    if(tissuename == "GI"){
      tissuenames2 = c("stomach", "intestine", "colon", "esophagus", "Stomach", "Intestine", "Colon", "Rectal", "rectal",
                       "Duodenum", "Gastric")
    }
    if(tissuename == "HRT"){
      tissuenames2 = c("heart", "cardiac", "coronary")
    }
    if(tissuename == "FAT"){
      tissuenames2 = c("adipo", "Adipo")
    }
    if(tissuename == "LIV"){
      tissuenames2 = c("liver", "hepato", "HepG2")
    }
    if(tissuename == "KID"){
      tissuenames2 = c("kidney")
    }
    if(tissuename == "SKIN"){
      tissuenames2 = c("skin", "keratin")
    }
    tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
    #! get the id of the selcted tissue, e.g. 69 tissue name from BLD
                                          
  }else if (tissuename == "ALL"){
    tissue_ids = 1:nrow(df_pre)
  }

  #! df_pre = AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz, may from ABC paper?        
  df = df_pre[tissue_ids, ] #! subset select the 69 tissue names
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)

  #! extract the enhancer cooridinate -> target gene name for selected tissues
  colnames(df2) = c("chr", "start", "end", "TargetGene") 
                                          
                                          
  matched_ids = match(df2$TargetGene, names(scores))
    #! get score array of the target gene in the selected tissues, e.g. 69 tissues         
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(df2[,c(1:3)], temp) #! cbind the enhancer-gene score for the selected tissue

#! read the road map data
  roadmap_meta = read.delim("../../processed_data/Roadmap_map_EID_names.txt",
                            header=F)

  if (tissuename != "ALL"){
    #! get the road id for the selected tissue
    Road_ids =  unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                               function(x) return(grep(x, roadmap_meta[,2])))), 1]))

    print(unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                         function(x) return(grep(x, roadmap_meta[,2])))), 3])))
  }else{
    Road_ids = unique(as.character(roadmap_meta[,1]))
  }

  Enhancer = c()
  for(ee in Road_ids){
    temp = read.table(paste0("../../processed_data//RoadmapLinks/",
                             "links_", ee, "_7_2.5.txt"))
    temp2 = read.table(paste0("../../processed_data//RoadmapLinks/",
                              "links_", ee, "_6_2.5.txt"))
      
      #! get the enhancer-target df from the road map
    Enhancer = rbind(Enhancer, temp[,1:4], temp2[,1:4]) #! from the file, 1:4 is chr, start, end, target
    cat("We processed file:", ee, "\n")
  }

  #! map the ENSG to gene symbol
  geneanno = read.csv("../../processed_data/gene_anno_unique_datefix.txt", sep = "\t")
  ff = geneanno$symbol[match(Enhancer[,4], geneanno$id)] #! n by 1 array/df
  tmp = cbind.data.frame(Enhancer, ff, 1) #! sames to be 1 in the last column?
                                                         
  dff = tmp[,c(1:3, 5, 6)] #! there is no 6 column? -> the 6th column is 1

  matched_ids = match(dff[,4], names(scores))
  temp = as.numeric(scores)[matched_ids] #! get the score array for the target gene array in the selected tissues
  temp[is.na(temp)] = 0

  final_bed2 = cbind(dff[,1:3], temp) #! add the enhancer chr, start, end, target gene score to the df
  colnames(final_bed1) = c("V1", "V2", "V3", "V4")
  colnames(final_bed2) = c("V1", "V2", "V3", "V4")
  final_bed = rbind(final_bed1, final_bed2)

  #! combine and write out the enhancer-gene score targets
  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


Yoshida_bedgraph_calc <- function(scores,
                                  output_cell,
                                  output_bed = "temp.bed"){

  df = data.frame(fread(paste0("../../processed_data",
                               "/", "Yoshida_correlated_summits_per_gene_inside100KB_sort.bed")))
  unique_genes = as.character(unique(df[,4]))
  gene_orthologs = read.table("../../processed_data/Orthologs_Yoshida_Mouse_Human.txt", header=T)

  temp = gene_orthologs[match(df[,4], gene_orthologs[,1]), 2]
  df2 = cbind.data.frame(df[which(!is.na(temp)), 1:3], temp[which(!is.na(temp))])
  colnames(df2) = c("chr", "start", "end", "gene")

  common_genes = intersect(names(scores), df2[,4])
  idx = which(!is.na(match(df2[,4], common_genes)))
  df3 = df2[idx,]
  grades = scores[match(df3$gene, names(scores))]
  bed1 = cbind(df3, grades)
  bed2 = bed1[,c(1:3, 5)]
  write.table(bed2, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

eQTLCPP_bedgraph_calc = function(scores,
                                 output_cell,
                                 output_bed = "temp.bed"){
  options(scipen = 999)

  Enhancer = read.table("../../processed_data/eQTL_Blood_CPP.txt",
                        header=F)
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

Roadmap_Enhancer_bedgraph_calc = function(scores,
                                          output_cell,
                                          output_bed = "temp.bed"){
  options(scipen = 999)

  Enhancer = read.table("../../processed_data/Roadmap_Enhancers_Blood.txt",
                        header=F)
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


