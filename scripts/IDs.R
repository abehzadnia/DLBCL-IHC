# ====================   ===
# Script 4: Heatmaps     ===
# AB                     ===
# 30-01-21               ===                
# ====================   ===
# Library               ====
# ====================   ===
library(tidyverse)
library(readxl)


# ====================   ===
# Data                  ====
# ====================   ===

path <- "../Prof Tooze - Onedrive/Prof Tooze/IHC_Scoring_TMAs.2_05-16.3_29.06_15.11.16_RT.xlsx"
rchop1 <- read_excel(path,sheet = "RCHOP 1", col_names = TRUE, skip=3)
rchop2 <- read_excel(path,sheet = "RCHOP 2", col_names = TRUE, skip=3)
rchop3 <- read_excel(path,sheet = "RCHOP 3", col_names = TRUE, skip=3)
rchopRQ <- read_excel(path,sheet = "R-CHOP RQ ", col_names = TRUE, skip=3)
rchopRQREF <- read_excel(path,sheet = "R-CHOP RQ REF", col_names = TRUE, skip=3)


# ====================   ===
# IDs                   ====
# ====================   ===

cleanRCHOP<- function(df,x){
        TMA_ID=factor(with(df, paste(Row, Column, sep="-")))
        HMDS_ID= factor(gsub("\\/","\\_", as.matrix(df[,1])))
        return(
                data.frame(BATCH=factor(rep(x, length(TMA_ID))),
                           TMA_ID, HMDS_ID) %>%
                        drop_na()
        )
}

IDs <- rbind(
        cleanRCHOP(rchop1, "RCHOP-1"),
        cleanRCHOP(rchop2, "RCHOP-2"),
        cleanRCHOP(rchop3, "RCHOP-3"),
        cleanRCHOP(rchopRQ,"RCHOP-RQ"),
        cleanRCHOP(rchopRQREF, "RCHOP-RQREF")
)
saveRDS(IDs, "./outputs/IDs.RDS")

