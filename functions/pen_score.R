# ====================   ===
# Function               ===
# Penalisation           ===                
# A Behzadnia            ===
# ====================   ===

# calculate penalised score for misidentification

pen_score <- function(wdf, stain = c("stain1","stain2")){
        ls <- list()
        x <- data.frame(STAIN = stain,
                        TUMOUR = NA,
                        TOTAL = NA,
                        T.SCORE = NA,
                        SCORE = NA)
        for (i in 1:length(stain)){
                ls[[i]] <- wdf %>% dplyr::filter(STAIN %in% stain[i]) %>%
                        dplyr::select(QP_Tumour)
                
                tum <- as.numeric(as.data.frame(ls[i])[[1]])
                tum <- tum[!is.na(tum)]

                x$TUMOUR[i] <-  length(tum[tum==1]) + 
                                length(tum[tum==2]) + 
                                length(tum[tum==3])
               
                x$TOTAL[i] <-  length(tum)
                
                x$T.SCORE[i] <- round(x$TUMOUR[i]/x$TOTAL[i],3)
                x$SCORE[i] <- round(x$TUMOUR[i]/(2*x$TOTAL[i]),3)
                        
                
                }
        return(x)
}
