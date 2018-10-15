library(GenomicRanges)
library(rtracklayer)


pmdFolder = "D:/Data/salhia_lab/desert_ID/temp_bg_files/"

pmdFiles = list.files(pmdFolder, full.names = TRUE)
pmdFiles = pmdFiles[grepl("filter", pmdFiles)]

rangeA <- import.bed(pmdFiles[1])
rangeB <- import.bed(pmdFiles[2])

overlap <- findOverlaps(rangeA, rangeB)
