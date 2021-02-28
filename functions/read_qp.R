# ====================   ===
# Function               ===
# Import QP outputs      ===                
# A Behzadnia            ===
# ====================   ===


# get IDs from master file
IDs <- readRDS("./outputs/IDs.RDS")

# import measurements from QP
read_qp <- function(stain=NULL,
                    batch = NULL){
        x <- paste("../QuPath_DLBCL_v2/Measurements/",batch,"/",
                   stain, ".qptma.data/TMA results - ",stain,".txt",sep="")
        y <- read.delim(x, header = T)
        z <- IDs %>% filter(BATCH == batch)  %>% select(TMA_ID)
        z$TMA_ID <- factor(z$TMA_ID)
        y <- filter(y, Name %in% levels(TMA_ID))
        return(y)
}

