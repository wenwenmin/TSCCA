# Load data
load("Result_ku10_kv100_kw20_J50.RData") # Output of TSCCA

# Enrichment analysis data
load("Fun2_InputData_EnrichmentAnalysis.RData")
source("Fun2_EnrichmentAnalysisFunctions.R")

# ------------------------------------------------------------------------------------
get_TSCCA_modules_information = function(DiffExp_Tensor.attribute, TSCCA.out){
  # gene Names
  gene.nameList1 = gene.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
  for(i in 1:length(DiffExp_Tensor.attribute$geneNames)){
    temp = unlist(strsplit(DiffExp_Tensor.attribute$geneNames[i], "[|]"))
    gene.nameList1[i] = temp[1]
    gene.nameList2[i] = temp[2] # Genes.entrezId
  }
  DiffExp_Tensor.attribute$geneNames_C1 = gene.nameList1
  DiffExp_Tensor.attribute$geneNames_C2 = gene.nameList2 # Genes.entrezId
  
  # Idnetified modules by TSCCA
  Gene.Modules = list()
  for(i in 1:length(TSCCA.out$d)){
    Gene.Modules[[i]] = gene.nameList2[which(TSCCA.out$V[,i]!=0)]
  }
  
  TSCCA_modules_information = list(Gene.Modules=Gene.Modules, AllConsideredGenes.entrezId=gene.nameList2)
  return(TSCCA_modules_information)
}
# ------------------------------------------------------------------------------------
TSCCA_modules_information = get_TSCCA_modules_information(DiffExp_Tensor.attribute, out)
Res = GeneModules.EnrichmentAnalysis(TSCCA_modules_information,known_functional_sets)

library("openxlsx")
fileName1 = "Number of Enriched terms"
fileName2 = "Enriched KEGG Terms"
fileName3 = "Enriched GOBP terms"
fileName4 = "Enriched reactome terms"

wb <- createWorkbook()
addWorksheet(wb, fileName1)
addWorksheet(wb, fileName2)
addWorksheet(wb, fileName3)
addWorksheet(wb, fileName4)

writeData(wb, sheet = fileName1, Res$modules.StatisticalResults)
writeData(wb, sheet = fileName2, Res$KEGGTerm.Table)
writeData(wb, sheet = fileName3, Res$GOBPTerm.Table)
writeData(wb, sheet = fileName4, Res$reactomeTerm.Table)
saveWorkbook(wb,"Fun2_Tab_enriched_anlysis.xlsx",overwrite=F)