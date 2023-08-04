
# R script to prepare gene sets

library(datamisc)

if (!file.exists("tables/genesets.tsv.gz")){
  
  GS <- list(GO = getGOgenes(minsize = 3)[,c("term","gene")],
             MitoCarta = getMitoCarta(),
             MetabolicAtlas = getMetabolicAtlas(),
             MSigDB = getMSigDB(format = "df")[,c("term","gene")])
  
  GS$MetabolicAtlas <- GS$MetabolicAtlas[sapply(GS$MetabolicAtlas, length) >= 3] |> convertGeneSets(to = "df")
  
  GS$GO$term <- paste0("GO ", GS$GO$term)
  GS$MitoCarta$term <- paste0("MitoCarta ", GS$MitoCarta$term)
  GS$MetabolicAtlas$term <- paste0("MetabolicAtlas ", sub("HUMAN1_", "", GS$MetabolicAtlas$term))
  GS$MSigDB$term <- paste0(ifelse(grepl("hallmark", GS$MSigDB$term, ignore.case = TRUE), "HALLMARK", "KEGG"),
                           gsub("HALLMARK_|KEGG_|_", " ", x = tolower(GS$MSigDB$term), ignore.case = TRUE))
  
  saveRDS(GS, file = file.path("data/rnaseq/GS.rds"))
  
  genesets <- Reduce(rbind, GS) |> unique()
  write.table(genesets, file = "tables/genesets.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
  system("gzip tables/genesets.tsv")
}

