# ====================   ===
# Function               ===
# Fleiss Kappa           ===                
# A Behzadnia            ===
# ====================   ===

# calculate interobservatory agreement between stain scores 
# raters = RT, MK and QP
# NB stroma scores missing in MK

immune_score_kappa <- function(df, 
                          raters = c("none","none2"), 
                          stain = c("none", "none2"),
                          test = c("Tumour", "Host Response"),
                          stroma = FALSE){
        x <- data.frame(RATER=rep(paste(raters,collapse="."),
                                  length(stain)*length(test)),
                        STAIN = rep(stain, length(test)),
                        TEST=rep(test, each=length(stain)))
        tum_key = 0
        hr_key = 0
        if("Tumour" %in% test == TRUE){
                tum_key = 1
                tumour <- list()
                for(i in 1:length(stain)){
                        tumour[[i]] <- kappam.fleiss(
                        df %>%
                                dplyr::filter(STAIN == stain[i]) %>%
                                dplyr::select(paste(raters,"Tumour",sep="_"))
                        )}
        }
        if("Host Response" %in% test == TRUE){
                hr_key = 1
                host <- list()
                for(i in 1:length(stain)){
                        host[[i]] <- kappam.fleiss(
                        df %>%
                                dplyr::filter(STAIN == stain[i]) %>%
                                dplyr::select(paste(raters,"HostResponse",sep="_"))
                )}
        } 
        if(stroma == TRUE){
                stroma = list()
                raters2 = c("RT", "QP")
                x <- rbind(x, data.frame(RATER = rep("RT.QP", length(stain)),
                                         STAIN = stain,
                                         TEST = rep("Stroma", length(stain))))
                for(i in 1:length(stain)){
                        stroma[[i]] <- kappam.fleiss(
                                df %>%
                                        dplyr::filter(STAIN == stain[i]) %>%
                                        dplyr::select(paste(raters2,"Stroma",sep="_"))
                        )
                }
                ls1 <- c(tumour, host, stroma)
                multiplier = 3
        } else if(stroma == FALSE){
                if(tum_key == 1 & hr_key == 1){
                        ls1 <- c(tumour, host)
                        multiplier = 2      
                }else if(tum_key == 1 & hr_key == 0){
                        ls1 <- c(tumour)
                        multiplier = 1 
                }else if(tum_key == 0 & hr_key ==1){
                        ls1 <- c(host)
                        multiplier = 1 
                }
                
        }
        
        for (k in 1:(length(stain)*multiplier)){
                for (j in 2:length(ls1[[1]])){
                        x[k,2+j] <- ls1[[k]][j]
                }
        }
        names(x)[5] <- "n.raters"
        names(x)[7] <- "Kappa"
        x[7] <- round(x[7],3)
        names(x)[9] <- "Z stat"
        x[9] <- round(x[9], 3)
        x[10] <- round(x[10], 4)
        return(x[-c(6,8)])
}
