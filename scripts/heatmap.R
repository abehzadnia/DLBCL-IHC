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
                get_ihc_score("CD3","RCHOP-2","Tumour"),
                get_ihc_score("CD3","RCHOP-2","HostResponse"),
                get_ihc_score("CD8","RCHOP-2","Tumour"),
                get_ihc_score("CD8","RCHOP-2","HostResponse"),
                get_ihc_score("CD68","RCHOP-2","Tumour"),
                get_ihc_score("CD68","RCHOP-2","HostResponse"),
                get_ihc_score("CD163","RCHOP-2","Tumour"),
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
                get_ihc_score("CD3","RCHOP-3","Tumour"),
                get_ihc_score("CD3","RCHOP-3","HostResponse"),
                get_ihc_score("CD8","RCHOP-3","Tumour"),
                get_ihc_score("CD8","RCHOP-3","HostResponse"),
                get_ihc_score("CD68","RCHOP-3","Tumour"),
                get_ihc_score("CD68","RCHOP-3","HostResponse"),
                get_ihc_score("CD163","RCHOP-3","Tumour"),
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
                get_ihc_percentage("CD3","RCHOP-2","Tumour"),
                get_ihc_percentage("CD3","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD8","RCHOP-2","Tumour"),
                get_ihc_percentage("CD8","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD68","RCHOP-2","Tumour"),
                get_ihc_percentage("CD68","RCHOP-2","HostResponse"),
                get_ihc_percentage("CD163","RCHOP-2","Tumour"),
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
               get_ihc_percentage("CD3","RCHOP-3","Tumour"),
                get_ihc_percentage("CD3","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD8","RCHOP-3","Tumour"),
                get_ihc_percentage("CD8","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD68","RCHOP-3","Tumour"),
                get_ihc_percentage("CD68","RCHOP-3","HostResponse"),
                get_ihc_percentage("CD163","RCHOP-3","Tumour"),
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
        mutate_at(c(2:17), function(x){as.numeric(as.character(x))}) %>%
        group_by(HMDS_ID) %>%
        summarise_all("mean",na.rm=T)

df_QP_ihc_scores_normalised <- data.frame(df_QP_ihc_scores_normalised[,1],
                scale2(df_QP_ihc_scores_normalised[,-1]))

for(i in 2:dim(df_QP_ihc_scores_normalised)[2]){
        colnames(df_QP_ihc_scores_normalised)[i] <- paste0("QP_",colnames(df_QP_ihc_scores_normalised)[i])
}

# combining IHC and gene data
df_heatmap <- merge(df_QP_ihc_scores_normalised, 
                    df_gep_HMDS[,c(1,57:20345)] , by = "HMDS_ID")

heatmap_matrix <- as.matrix(df_heatmap[,-1])
mode(heatmap_matrix)<-"numeric"

# heatmap

tiff(filename="./figures/heatmap_IHC_SCORES.tiff", width = 40000, height = 10000,
     res=1000)
heatmap(heatmap_matrix, scale = "none")
dev.off() # ?? doesn't look right!!


# heatmap 2
tiff(filename="./figures/pheatmap_SCORES_10000.tiff", 
     res=1000, width = 40000, height = 10000)
pheatmap(heatmap_matrix[,c(1:10000)],fontsize = 0.5, scale = "none")
dev.off()


# ========================= =
# Heatmap -IHCPercentage ====
# ======================-== =

# normalisation
df_QP_ihc_percentages_normalised <- 
        df_QP_ihc_percentages %>%
        select(-TMA_ID,-BATCH) %>%
        mutate_at(c(2:17), function(x){as.numeric(as.character(x))}) %>%
        group_by(HMDS_ID) %>%
        summarise_all("mean",na.rm=T)

df_QP_ihc_percentages_normalised <- data.frame(df_QP_ihc_percentages_normalised[,1],
                                          scale2(df_QP_ihc_percentages_normalised[,-1]))

for(i in 2:dim(df_QP_ihc_percentages_normalised)[2]){
        colnames(df_QP_ihc_percentages_normalised)[i] <- paste0("QP_",colnames(df_QP_ihc_percentages_normalised)[i])
}


# combining IHC and gene data
df_heatmap <- merge(df_QP_ihc_percentages_normalised, 
                    df_gep_HMDS[,c(1,57:20345)] , by = "HMDS_ID")

heatmap_matrix <- as.matrix(df_heatmap[,-1])
mode(heatmap_matrix)<-"numeric"

# heatmap

tiff(filename="./figures/heatmap_IHC_PERCENTAGE.tiff", width = 40000, height = 10000,
     res=1000)
heatmap(heatmap_matrix, scale = "none")
dev.off()


