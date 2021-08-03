coloradoMRI
================
Jadyn Park
5/25/2021

Cleaning CO and NU DTI stats

``` r
## import data files
filelist=list.files("~/Downloads/FrontalSub_MRI_coloradoMRI/DTI_Stats/", pattern = ".txt",
                      full.names=T, recursive=T)

## assuming tab separated values with a header; creating one large list    
datalist=lapply(filelist, FUN=read.table, header=TRUE) ## datalist[[1]], [[2]], [[3]], ... etc. to view individual dataset

## assuming the same header/columns for all files
datafr=do.call("cbind", datalist) 

## "lh_unc" and "lh.unc" are duplicates, going with "lh.unc"
## "lh_ifl" and "lh.ilf" are duplicates, going with "lf.ilf"
## "rh_unc" and "rh.unc" are duplicates, going with "rh.unc"
## "rh_ifl" and "rh.ilf" are duplicates, going with "rh.ilf"

require(openxlsx)
wb <- list("fmajor"=datalist[[1]], "fminor"=datalist[[2]],
           "lh.atr"=datalist[[3]], "lh.cab"=datalist[[4]],
           "lh.ccg"=datalist[[5]], "lh.cst"=datalist[[6]],
           "lh.ilf"=datalist[[7]], "lh.slfp"=datalist[[8]],
           "lh.slft"=datalist[[9]], "lh.unc"=datalist[[10]],
           "rh.atr"=datalist[[11]], "rh.cab"=datalist[[12]],
           "rh.ccg"=datalist[[13]], "rh.cst"=datalist[[14]],
           "rh.ilf"=datalist[[15]], "rh.slfp"=datalist[[16]],
           "rh.slft"=datalist[[17]],"rh.unc"=datalist[[18]])

## export as a single workbook, with one structure under each sheet
write.xlsx(wb, file = "colorado-DTI-allsub.xlsx")
             
## export as a single workbook, with all structure in a single sheet
fmajor <- datalist[[1]]
colnames(fmajor) <- paste0("fmajor_", colnames(fmajor))
fminor <- datalist[[2]]
colnames(fminor) <- paste0("fminor_", colnames(fminor))

lh.atr <- datalist[[3]]
colnames(lh.atr) <- paste0("lh.atr_", colnames(lh.atr))
lh.cab <- datalist[[4]]
colnames(lh.cab) <- paste0("lh.cab_", colnames(lh.cab))
lh.ccg <- datalist[[5]]
colnames(lh.ccg) <- paste0("lh.ccg_", colnames(lh.ccg))
lh.cst <- datalist[[6]]
colnames(lh.cst ) <- paste0("lh.cst_", colnames(lh.cst ))
lh.ilf <- datalist[[7]]
colnames(lh.ilf) <- paste0("lh.ilf_", colnames(lh.ilf))
lh.slfp <- datalist[[8]]
colnames(lh.slfp ) <- paste0("lh.slfp_", colnames(lh.slfp ))
lh.slft <- datalist[[9]]
colnames(lh.slft) <- paste0("lh.slft_", colnames(lh.slft))
lh.unc <- datalist[[10]]
colnames(lh.unc) <- paste0("lh.unc_", colnames(lh.unc))

rh.atr <- datalist[[11]]
colnames(rh.atr) <- paste0("rh.atr_", colnames(rh.atr))
rh.cab <- datalist[[12]]
colnames(rh.cab) <- paste0("rh.cab_", colnames(rh.cab))
rh.ccg <- datalist[[13]]
colnames(rh.ccg) <- paste0("rh.ccg_", colnames(rh.ccg))
rh.cst <- datalist[[14]]
colnames(rh.cst) <- paste0("rh.cst_", colnames(rh.cst))
rh.ilf <- datalist[[15]]
colnames(rh.ilf) <- paste0("rh.ilf_", colnames(rh.ilf))
rh.slfp <- datalist[[16]]
colnames(rh.slfp) <- paste0("rh.slfp_", colnames(rh.slfp))
rh.slft <- datalist[[17]]
colnames(rh.slft) <- paste0("rh.slft_", colnames(rh.slft))
rh.unc <- datalist[[18]]
colnames(rh.unc) <- paste0("rh.unc_", colnames(rh.unc))

allsub <- cbind(fmajor, fminor, lh.atr, lh.cab, lh.ccg, lh.cst, lh.ilf, lh.slfp, lh.slft,lh.unc, 
                rh.atr, rh.cab, rh.ccg, rh.cst, rh.ilf, rh.slfp, rh.slft, rh.unc)

## remove subid variables
allsub <- allsub %>% select(-fminor_fminor, -lh.atr_lh.atr, -lh.cab_lh.cab, -lh.ccg_lh.ccg,-lh.cst_lh.cst, 
                            -lh.ilf_lh.ilf, -lh.slfp_lh.slfp, -lh.slft_lh.slft,-lh.unc_lh.unc,
                            -rh.atr_rh.atr, -rh.cab_rh.cab, -rh.ccg_rh.ccg,-rh.cst_rh.cst, 
                            -rh.ilf_rh.ilf, -rh.slfp_rh.slfp, -rh.slft_rh.slft,-rh.unc_rh.unc)

write.xlsx(allsub, file = "colorado-DTI-master.xlsx")
```
