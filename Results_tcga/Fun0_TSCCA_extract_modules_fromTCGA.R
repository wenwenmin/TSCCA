library("irlba")
library("ggplot2")
library("cowplot")
library("openxlsx")

# load("Data_tcga/TCGA_Tensor_Data.RData")
load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")

# TSCCA function
source('TSCCA.R') 

# parameters of TSCCA
ku = 10
kv = 100
kw = 20
J  = 50 # 50 modules

# Apply TSCCA to tcga data across 33 cancer types
ptm = proc.time()
out = TSCCA(DiffExp_Tensor, ku, kv, kw, J)
time1 = proc.time() - ptm; print(time1)
ptm = proc.time()

# save
tscca.parameter = list(ku=ku,kv=kv,kw=kw,J=J)
tscca.out = out
rm(DiffExp_Tensor)
save.image(file = "Result_ku10_kv100_kw20_J50.RData")


