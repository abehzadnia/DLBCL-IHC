# ====================   ===
# Script 4: Heatmaps     ===
# AB                     ===
# Ver1; 28-02-21         ===                
# ====================   ===
# Library               ====
# ====================   ===
# packages
library(tidyverse)
library(readxl)
library(pheatmap)

# functions
# global
source("./functions/immune_score1.R")
source("./functions/ihc_score.R")
source("./functions/ihc_percentage.R")
source("./functions/immune_score_kappa.R")
source("./functions/pen_score.R")

# local
get_ihc_score<-function(stain="x",batch="y",score="z" ){
        df1<-ihc_score(read_qp(stain,batch),stain ,score)
        df2<-merge(dplyr::filter(IDs, BATCH == batch), df1)
        return(df2)
}

get_ihc_percentage<-function(stain="x",batch="y",score="z" ){
        df1<-ihc_percentage(read_qp(stain,batch),stain ,score)
        df2<-merge(dplyr::filter(IDs, BATCH == batch), df1)
        return(df2)
}

# ========================= =
# QP IHC SCORE Data      ====
# ======================-== =
IDs<-readRDS("./outputs/IDs.RDS")

df_QP_ihc_scores <- 
        rbind(
        plyr::join_all(list(
        #RCHOP-2
                get_ihc_score("CD3","RCHOP-2","HostResponse"),
                get_ihc_score("CD8","RCHOP-2","HostResponse"),
                get_ihc_score("CD68","RCHOP-2","HostResponse"),
                get_ihc_score("CD163","RCHOP-2","HostResponse"),
                get_ihc_score("PD1","RCHOP-2","Tumour"),
                get_ihc_score("PD1","RCHOP-2","HostResponse"),
                get_ihc_score("PDL1","RCHOP-2","Tumour"),
                get_ihc_score("PDL1","RCHOP-2","HostResponse"),
                get_ihc_score("IRF1","RCHOP-2","Tumour"),
                get_ihc_score("IRF1","RCHOP-2","HostResponse"),
                get_ihc_score("STAT1P","RCHOP-2","Tumour"),
                get_ihc_score("STAT1P","RCHOP-2","HostResponse")
        )),
        
        plyr::join_all(list(
        #RCHOP-3
                get_ihc_score("CD3","RCHOP-3","HostResponse"),
                get_ihc_score("CD8","RCHOP-3","HostResponse"),
                get_ihc_score("CD68","RCHOP-3","HostResponse"),
                get_ihc_score("CD163","RCHOP-3","HostResponse"),
                get_ihc_score("PD1","RCHOP-3","Tumour"),
                get_ihc_score("PD1","RCHOP-3","HostResponse"),
                get_ihc_score("PDL1","RCHOP-3","Tumour"),
                get_ihc_score("PDL1","RCHOP-3","HostResponse"),
                get_ihc_score("IRF1","RCHOP-3","Tumour"),
                get_ihc_score("IRF1","RCHOP-3","HostResponse"),
                get_ihc_score("STAT1P","RCHOP-3","Tumour"),
                get_ihc_score("STAT1P","RCHOP-3","HostResponse")
        ))
)


# ========================= =
# QP IHC PERCENTAGE Data ====
# ======================-== =

df_QP_ihc_percentages <- 
        rbind(
        plyr::join_all(list(
        #RCHOP-2
                get_ihc_percentage("CD3","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD8","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD68","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD163","RCHOP-2","HostResponse"),
                get_ihc_percentage("PD1","RCHOP-2","Tumour"),
                get_ihc_percentage("PD1","RCHOP-2","HostResponse"),
                get_ihc_percentage("PDL1","RCHOP-2","Tumour"),
                get_ihc_percentage("PDL1","RCHOP-2","HostResponse"),
                get_ihc_percentage("IRF1","RCHOP-2","Tumour"),
                get_ihc_percentage("IRF1","RCHOP-2","HostResponse"),
                get_ihc_percentage("STAT1P","RCHOP-2","Tumour"),
                get_ihc_percentage("STAT1P","RCHOP-2","HostResponse")
                )),
                
        plyr::join_all(list(
        #RCHOP-3
                get_ihc_percentage("CD3","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD8","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD68","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD163","RCHOP-3","HostResponse"),
                get_ihc_percentage("PD1","RCHOP-3","Tumour"),
                get_ihc_percentage("PD1","RCHOP-3","HostResponse"),
                get_ihc_percentage("PDL1","RCHOP-3","Tumour"),
                get_ihc_percentage("PDL1","RCHOP-3","HostResponse"),
                get_ihc_percentage("IRF1","RCHOP-3","Tumour"),
                get_ihc_percentage("IRF1","RCHOP-3","HostResponse"),
                get_ihc_percentage("STAT1P","RCHOP-3","Tumour"),
                get_ihc_percentage("STAT1P","RCHOP-3","HostResponse")
                ))
        )



# ========================= =
# GENE Exp               ====
# ======================-== =
path_new <- "./data/2021.01.25_IHC_ExpressionData/ihc_expression_new_gene_Final.txt"
#path_GSE <- "./data/2021.01.25_IHC_ExpressionData/ihc_expression_GSE32918_gene_Final.txt"

dat_new <- read.delim(path_new)
#dat_GSE <- read.delim(path_GSE) NOT YET!


# Gene expressions only
df_gep <- dat_new[-c(2:26),]

# transpose, df and colnames
df_gep <- t(df_gep)
cols <- df_gep[1,]
df_gep <- data.frame(df_gep[-1,])
names(df_gep) <- cols

# selecting cases
df_gep_HMDS <- df_gep[df_gep$HMDS_request_ID 
                      %in% levels(df_QP_ihc_scores$HMDS_ID),]
#correct names and remove rownames
colnames(df_gep_HMDS)[1] <- "HMDS_ID"
rownames(df_gep_HMDS) <- NULL


# ========================= =
# Heatmap - IHC Score    ====
# ======================-== =

# normalisation
scale2<-function(x){
        m <- apply(x,2,mean, na.rm=T)
        sd <- apply(x,2, sd, na.rm=T)
        z <- scale(x,m,sd)
        return(z)
}

df_QP_ihc_scores_normalised <- 
        df_QP_ihc_scores %>%
        select(-TMA_ID,-BATCH) %>%
        mutate_at(c(2:13), function(x){as.numeric(as.character(x))}) %>%
        group_by(HMDS_ID) %>%
        summarise_all("mean",na.rm=T)

df_QP_ihc_scores_normalised <- data.frame(df_QP_ihc_scores_normalised[,1],
                scale2(df_QP_ihc_scores_normalised[,-1]))

for(i in 2:dim(df_QP_ihc_scores_normalised)[2]){
        colnames(df_QP_ihc_scores_normalised)[i] <- paste0("QP_",colnames(df_QP_ihc_scores_normalised)[i])
}

# combining IHC and gene data
df_QP_ihc_scores_normalised <- df_QP_ihc_scores_normalised[complete.cases(df_QP_ihc_scores_normalised),]

df_heatmap <- merge(df_QP_ihc_scores_normalised, 
                    df_gep_HMDS[,c(1,57:20345)] , by = "HMDS_ID")

# # # after email on Monday 1/3/21 # # # # 

df_heatmap_ihc <- as.matrix(df_heatmap[,2:13])
mode(df_heatmap_ihc) <- "numeric"

df_heatmap_gene <- as.matrix(df_heatmap[,-c(1:13)])
mode(df_heatmap_gene) <- "numeric"

cor_matrix <- cor(df_heatmap_ihc, df_heatmap_gene)

write.csv(cor_matrix, "./outputs/cor_matrix_IHC-SCORE.csv")



###### heatmap
cor_matrix_trim <- cor_matrix
#cor_matrix_trim[!lower.tri(cor_matrix_trim)] <- NA # remove diagonal
cor_matrix_trim <- 
        data.frame(cor_matrix_trim) %>%
        rownames_to_column() %>%
        gather(key="variable", value="correlation", -rowname) %>%
        filter(abs(correlation) > 0.5)

cor_matrix_trim <- spread(cor_matrix_trim, variable, correlation)
row_names <- cor_matrix_trim[,1]
cor_matrix_trim <- as.matrix(cor_matrix_trim[,-1])
rownames(cor_matrix_trim) <- row_names
cor_matrix_trim[is.na(cor_matrix_trim)] <- 0


# heatmap
tiff(filename="./figures/heatmap_IHC_SCORE_COR0.5.tiff", width = 20000, height = 2000,
     res=400)
pheatmap(cor_matrix_trim, scale = "none", fontsize = 5)
#pheatmap(cor_matrix,fontsize = 1, scale = "none")
dev.off()





# ========================= =
# Heatmap -IHCPercentage ====
# ======================-== =

# normalisation
df_QP_ihc_percentages_normalised <- 
        df_QP_ihc_percentages %>%
        select(-TMA_ID,-BATCH) %>%
        mutate_at(c(2:13), function(x){as.numeric(as.character(x))}) %>%
        group_by(HMDS_ID) %>%
        summarise_all("mean",na.rm=T)

df_QP_ihc_percentages_normalised <- data.frame(df_QP_ihc_percentages_normalised[,1],
                                               scale2(df_QP_ihc_percentages_normalised[,-1]))

for(i in 2:dim(df_QP_ihc_percentages_normalised)[2]){
        colnames(df_QP_ihc_percentages_normalised)[i] <- paste0("QP_",colnames(df_QP_ihc_percentages_normalised)[i])
}



# combining IHC and gene data
df_QP_ihc_percentages_normalised <- df_QP_ihc_percentages_normalised[complete.cases(df_QP_ihc_percentages_normalised),]

df_heatmap_perc <- merge(df_QP_ihc_percentages_normalised, 
                    df_gep_HMDS[,c(1,57:20345)] , by = "HMDS_ID")

# # # after email on Monday 1/3/21 # # # # 

df_heatmap_ihc_perc <- as.matrix(df_heatmap_perc[,2:13])
mode(df_heatmap_ihc_perc) <- "numeric"

df_heatmap_gene <- as.matrix(df_heatmap[,-c(1:13)])
mode(df_heatmap_gene) <- "numeric"

cor_matrix_perc <- cor(df_heatmap_ihc_perc, df_heatmap_gene)

write.csv(cor_matrix_perc, "./outputs/cor_matrix_IHC-PERCENTAGE.csv")



###### heatmap
cor_matrix_perc_trim <- cor_matrix_perc
#cor_matrix_trim[!lower.tri(cor_matrix_trim)] <- NA # remove diagonal
cor_matrix_perc_trim <- 
        data.frame(cor_matrix_perc_trim) %>%
        rownames_to_column() %>%
        gather(key="variable", value="correlation", -rowname) %>%
        filter(abs(correlation) > 0.5)

cor_matrix_perc_trim <- spread(cor_matrix_perc_trim, variable, correlation)
row_names <- cor_matrix_perc_trim[,1]
cor_matrix_perc_trim <- as.matrix(cor_matrix_perc_trim[,-1])
rownames(cor_matrix_perc_trim) <- row_names
cor_matrix_perc_trim[is.na(cor_matrix_perc_trim)] <- 0


# heatmap
tiff(filename="./figures/heatmap_IHC_PERCENTAGE_COR0.5.tiff", 
     width = 20000, height = 2000,
     res=400)
pheatmap(cor_matrix_trim, scale = "none", fontsize = 5)

dev.off()




# PIR GENES

pir_genes <- c("TRAT1","CD3D","CD3G","GZMK","ITM2A",
               "GIMAP6","CD2","CLEC2B","GZMA","TC2N","FGL2","IFNG",
               "BCL11B","UBASH3A")

df_gep_pir <- df_gep_HMDS %>% select("HMDS_ID",paste(pir_genes))

df_gep_pir <- merge(df_QP_ihc_percentages_normalised, 
                         df_gep_pir , by = "HMDS_ID")

#### 
df_pir_ihc_perc <- as.matrix(df_gep_pir[,2:13])
mode(df_pir_ihc_perc) <- "numeric"

df_pir_gene <- as.matrix(df_gep_pir[,-c(1:13)])
mode(df_pir_gene) <- "numeric"

cor_matrix_pir <- cor(df_pir_ihc_perc, df_pir_gene)


###### heatmap
cor_matrix_pir_trim <- cor_matrix_pir
#cor_matrix_trim[!lower.tri(cor_matrix_trim)] <- NA # remove diagonal
cor_matrix_perc_trim <- 
        data.frame(cor_matrix_perc_trim) %>%
        rownames_to_column() %>%
        gather(key="variable", value="correlation", -rowname) %>%
        filter(abs(correlation) > 0.5)

cor_matrix_perc_trim <- spread(cor_matrix_perc_trim, variable, correlation)
row_names <- cor_matrix_perc_trim[,1]
cor_matrix_perc_trim <- as.matrix(cor_matrix_perc_trim[,-1])
rownames(cor_matrix_perc_trim) <- row_names
cor_matrix_perc_trim[is.na(cor_matrix_perc_trim)] <- 0


# heatmap
tiff(filename="./figures/heatmap_IHC_PIR.tiff", 
     width = 2000, height = 2000,
     res=400)
pheatmap(cor_matrix_pir_trim, scale = "none", fontsize = 5)
dev.off()



