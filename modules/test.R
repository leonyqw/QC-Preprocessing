library("dplyr")
library("data.table")
library("seqinr")
library("plyranges")
library("Rsamtools")


bam <- scanBam("./results/1. aligned reads/barcode03_aligned_sorted.bam", 
index = "./results/1. aligned reads/barcode03_aligned_sorted.bam.bai")

names <- data.frame(id = bam[[1]]$qname)

path <- system.file("./results/2. extracted reads/barcode03_heavy.fasta",
package = "seqinr")
heavy <- read.fasta(file = "./results/2. extracted reads/barcode03_heavy.fasta")
head(heavy)
light <- read.fasta(file = "./results/2. extracted reads/barcode03_light.fasta")

heavy_id <- data.frame(name = attr(heavy, "name"), present = "TRUE")
light_id <- data.frame(name = attr(light, "name"), present = "TRUE")

counts <- names %>% group_by(id) %>% summarise(n = n())
counts %>% filter(n != 1)
names[names$id == "1586d83d-dd66-41fe-906b-7b50114948f7",]
join <- full_join(heavy_id, light_id, by = join_by(name == name))

head(heavy_id)
