# ====================   ===
# Function               ===
# IHC Percentage         ===                
# A Behzadnia            ===
# ====================   ===

ihc_percentage <- function(df, stain = "none", 
                          score = c("Tumour", "HostResponse")){
        # require(tibble)
        x = data.frame(ID = df[2])
        tumour = NA
        host = NA
       
        if (score == "Tumour"){
                if("Tumor..Positive.." %in% names(df) == FALSE){
                        tumour = 0
                        y <- tumour
                } else{
                        y <- df$Tumor..Positive..   
                }
                
        } else if(score == "HostResponse"){
                if("Immune.cells..Positive.." %in% names(df) == FALSE){
                        host = 0
                        y <- host
                } else{
                        y <- df$Immune.cells..Positive..  
                }
                
        }
        x$Percentage <- round(y,2)
        colnames(x)[1] <- "TMA_ID"
        colnames(x)[2] <- paste(stain, score, sep = "_" )
        return(x[c(1,2)])
}
