# ====================   ===
# Function               ===
# IHC Scoring            ===                
# A Behzadnia            ===
# ====================   ===

ihc_score <- function(df, stain = "none", 
                          score = c("Tumour", "HostResponse")){
        # require(tibble)
        x = data.frame(ID = df[2])
        tumour = NA
        host = NA
       
        if (score == "Tumour"){
                if("Tumor..Positive.." %in% names(df) == FALSE){
                        tumour = 0
                        y <- tumour/df$Num.Tumor..base.
                } else{
                        y <- df$Tumor..Positive..   
                }
                
        } else if(score == "HostResponse"){
                if("Immune.cells..Positive.." %in% names(df) == FALSE){
                        host = 0
                        y <- host/df$Num.Immune.cells..base.
                } else{
                        y <- df$Immune.cells..Positive..  
                }
                
        }
        x$Ratio <- y
        x$Score <- cut(x$Ratio, 
                       breaks = c(0,1,25,50,100), include.lowest = T,
                       labels = c(0,1,2,3))
        colnames(x)[1] <- "TMA_ID"
        colnames(x)[3] <- paste(stain, score, sep = "_" )
        return(x[c(1,3)])
}
