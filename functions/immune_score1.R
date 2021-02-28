# ====================   ===
# Function               ===
# Immune Scoring         ===                
# A Behzadnia            ===
# ====================   ===

# calculate immune scoring
# PD1 low threshold to be scored as 1
immune_score1 <- function(df, stain = "none", 
                         score = c("Tumour", "HostResponse")){
        # require(tibble)
        x = data.frame(ID = df[2])
        tumour = NA
        host = NA
        # col8 = 0
        # col9 = 0
        # col11 = 0
        # col12 = 0
        # 
        # if("Num.Negative..Immune.cells..base." %in% names(df) == FALSE){
        #         col8 <- 1
        # }
        # if("Num.Negative..Tumor..base." %in% names(df) == FALSE){
        #         col9 <- 1
        # }
        # if("Num.Positive..Immune.cells..base." %in% names(df) == FALSE){
        #         col11 <- 1
        # }
        # if("Num.Positive..Tumor..base." %in% names(df) == FALSE){
        #         col12 <- 1
        # }
        # 
        # if(col8 == 1 || col9 == 1 ||
        #     col11 == 1 || col12 == 1){
        #         if(add_col == T){
        #                 if(col8 ==1){df <- add_column(df, col8 = 0, .after =7)}
        #                 if(col9 ==1){df <- add_column(df, col9 = 0, .after =8)}
        #                 if(col11 == 1){df<- add_column(df, col11 = 0, .after = 10)}
        #                 if(col12 == 1){df<-add_column(df, col12 = 0, .after = 11)}      
        #         }else(return("colnames do not match, add_col should be TRUE"))
        # }

        if (score == "Tumour"){
                if("Num.Positive..Tumor..Positive" %in% names(cd3_qp) == FALSE){
                        tumour = 0
                        y <- tumour/(tumour+df$Num.Negative..Tumor..Negative)
                } else{
                        y <- df$Num.Positive..Tumor..Positive/(df$Num.Positive..Tumor..Positive
                                                               +df$Num.Negative..Tumor..Negative)    
                }
                
        } else if(score == "HostResponse"){
                y <- df$Num.Positive..Immune.cells..Positive/(df$Num.Positive..Immune.cells..Positive
                                                              +df$Num.Negative..Immune.cells..Negative)
        }
        x$Ratio <- y
        x$Score <- cut(x$Ratio, 
                       breaks = c(0,0.01,0.25,0.5,1), include.lowest = T,
                       labels = c(0,1,2,3))
        colnames(x)[1] <- "TMA_ID"
        colnames(x)[3] <- paste(stain, score, sep = "_" )
        return(x[c(1,3)])
}




