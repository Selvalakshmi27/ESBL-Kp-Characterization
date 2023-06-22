pacman::p_load(
  tidyverse,
  ape,
  ggtree,
  rhierbaps,
  phytools,
  ggnewscale,
  ggtreeExtra,
  ggpubr,
  scales, 
  ggplot2,
  RColorBrewer
)
options("install.lock"=FALSE)

#full
fasta<- read.FASTA("core_gene_alignment_full.aln")
snp.matrix <- load_fasta(fasta)

floor(nrow(snp.matrix)/5)

hb.results <- hierBAPS(snp.matrix,  max.depth = 2, n.pops = 26)
head(hb.results$partition.df)
tree<- read.tree("core_gene_alignment_full.aln.treefile")
tree <- midpoint.root(tree)

gg <- ggtree(tree, layout = "circular") + 
  geom_treescale(fontsize = 2.5, linesize = 0.01, family = "Arial Bold")


#Extract node label
bootsv <-data.frame(tree[["node.label"]])
ntip <- length(tree[["tip.label"]])
nnode <- tree[["Nnode"]]
SHaLRT <- data.frame(lapply(bootsv, function(x) as.numeric(sub("/.*", "",x))))
UFBoot <-data.frame(lapply(bootsv, function(x) as.numeric(sub(".*/", "",x))))
support <- cbind(SHaLRT, UFBoot)
support
colnames(support) <- c("SHaLRT", "UFBoot")

row.names(support) <- c((ntip + 1) : (ntip + nnode))
#Remove values below SH-aLRT80 <80 and UFB <90
support <- support[which(support$tree...node.label... >= 95 & support$tree...node.label... >=90),]
#cut data
support$UFBoot <- cut(as.numeric(support$UFBoot), breaks = 90,right = F,include.lowest = TRUE)

gg <- tree2 %<+% hb.results$partition.df
gg <- gg + geom_tippoint(aes(color = factor(`level 1`)), size=2.8) 
gg
gg1<- gg +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_color_discrete(name = "Clonal Groups",labels=c("CG307","Rare", "CG29","CG258","CG392","CG15","CG2947")) +
  geom_hilight(node=145, fill="gold",alpha=0.1,extend=0.0033) +
  geom_hilight(node=191, fill="#90e0ef",alpha=0.2,extend=0.0018) +
  geom_hilight(node=199, fill="darkolivegreen3",alpha=0.2,extend=0.0022) +
  geom_hilight(node=217, fill="#c8b6ff",alpha=0.24,extend=0.0031)+
  geom_hilight(node=247, fill="#adb5bd",alpha=0.2,extend=0.0015) +
  geom_hilight(node=239, fill="hot pink",alpha=0.2,extend=0.0025) +
  geom_cladelab(145,"CG307",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.0001,offset.text=0.0000009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T) +
  geom_cladelab(191,"CG258",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.00009,offset.text=0.0000009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T) +
  geom_cladelab(199,"CG392",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.00009,offset.text=0.0000009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T)+
  geom_cladelab(217,"CG29",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.000009,offset.text=0.000009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T) +
  geom_cladelab(247,"CG2947",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.0005,offset.text=0.00009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T)+
  geom_cladelab(239,"CG15",barsize=0,barcolor=NA,
                hjust=.5,angle=0,offset=0.00009,offset.text=0.0000009,
                hjust=0.5,fontcolor='red', fontsize=3, fontface="bold",align=T)

gg1
meta_data1 <- read.csv("Unique_isolates_kleborate.csv")
meta_data1$strain<-gsub("_spades_scaffolds", "", meta_data1$strain, fixed=TRUE)
meta_data1$strain<-gsub("_prokka", "", meta_data1$strain, fixed=TRUE)

meta_data1
df1 <- as.data.frame(meta_data1)
df1[2] <- NULL
rownames(df1) <- df1$strain
df1






gg1
table(df1$Bla_ESBL_acquired)
df1$Bla_ESBL_acquired<- gsub("CTX-M-15*;CTX-M-15*;CTX-M-15*", "CTX-M-15", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("CTX-M-15;CTX-M-15;CTX-M-15", "CTX-M-15", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("CTX-M-15;CTX-M-15", "CTX-M-15", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("CTX-M-15^", "CTX-M-15", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("CTX-M-14^;CTX-M-15", "CTX-M-15", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("SHV-31.v1^", "SHV-31", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("sHV-7^", "SHV-7", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("SHV-7^", "SHV-7", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired<- gsub("TEM-116*", "TEM-116", df1$Bla_ESBL_acquired, fixed=TRUE)
df1$Bla_ESBL_acquired[df1$Bla_ESBL_acquired == ""] <- NA
df1$Bla_ESBL_acquired[df1$Bla_ESBL_acquired == "-"] <- NA


p1<-gheatmap(gg1,df1[, "Bla_ESBL_acquired",drop=FALSE],offset=0.0013, width=0.06, hjust = 0.9,colnames_position = "top",
             color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("ESBL"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#00bbf9", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E", "#E6AB02", "#ef233c", "#fb8500","purple","#6d2e46","#e0b1cb"),
                    name = "ESBL",
                    breaks =c("CTX-M-15", "CTX-M-27","CTX-M-3", "CTX-M-55", "CTX-M-65", "SHV-31", "SHV-7", "TEM-116","SHV-100","SHV-30","CTX-M-27"), na.value="#edede9" )

p1
df1$Bla_Carb_acquired
df1$Bla_Carb_acquired[df1$Bla_Carb_acquired == "-"] <- NA
df1$Bla_Carb_acquired[df1$Bla_Carb_acquired == ""] <- NA

p2<-p1+new_scale_fill()
p3<- gheatmap(p2, df1[, "Bla_Carb_acquired", drop=FALSE],offset=0.0016, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("Carb"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#00bbf9", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"),
                    name = "Carbapenems",
                    breaks =c( "KPC-2", "NDM-1;OXA-48", "OXA-181" , "OXA-232", "OXA-48"),na.value = "#edede9")


p3

table(df1$Bla_OXA)
df1$Bla_OXA[df1$Bla_OXA == ""] <- NA
p4<-p3+new_scale_fill()
p5<- gheatmap(p4, df1[, "Bla_OXA", drop=FALSE],offset=0.00192, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("OXA-1"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#A6D854"),
                    name = "Other",
                    breaks =c( "OXA-1"),na.value = "#edede9")

p5
df1$Bla_TEM[df1$Bla_TEM == ""] <- NA
p6<-p5+new_scale_fill()
p7<- gheatmap(p6, df1[, "Bla_TEM", drop=FALSE],offset=0.00225, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("TEM-1"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#A6D854"),
                    name = "Other",
                    breaks =c( "TEM-1"),na.value = "#edede9")

p7







table(df1$GyrA.83)
df1$GyrA.83[df1$GyrA.83 == ""] <- NA
p4<-gg1+new_scale_fill()
p5<- gheatmap(p4, df1[, "GyrA.83", drop=FALSE], offset=0.0013, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("GyrA-83"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#8338ec","#ffbe0b","#ff006e"),
                    name = "GyrA-83",
                    breaks =c("GyrA-83F", "GyrA-83I", "GyrA-83Y"),na.value = "#edede9")
p5


table(df1$GyrA.87)
df1$GyrA.87[df1$GyrA.87 == ""] <- NA

p6<-p5+new_scale_fill()
p7<- gheatmap(p6, df1[,"GyrA.87", drop=FALSE],offset=0.0016, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("GyrA-87"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#8338ec","#ffbe0b","#ff006e"),
                    name = "GyrA-87",
                    breaks =c("GyrA-87A", "GyrA-87N", "GyrA-87G"),na.value = "#edede9")
p7

table(df1$ParC)
df1$ParC[df1$ParC == ""] <- NA
df1$ParC[df1$ParC == "-"] <- NA
p8<-p7+new_scale_fill()
p9<- gheatmap(p8, df1["ParC", drop=FALSE],offset=0.0019, width=0.06, hjust = 1.0,colnames_position = "top",
              color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("ParC"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#8338ec","#ffbe0b"),
                    name = "ParC",
                    breaks =c("ParC-84K","ParC-80I"),na.value = "#edede9")

p9
table(df1$Flq_acquired)
df1$Flq_acquired[df1$Flq_acquired == "qnrA1"] <- "qnrA"
df1$Flq_acquired[df1$Flq_acquired == "qnrA1^"] <- "qnrA"
df1$Flq_acquired[df1$Flq_acquired == "qnrB1"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB1.v1"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB1.v1;qnrS1"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB1.v2^"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB1.v2^;qnrE1"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB2"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB4"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrB6^"] <- "qnrB"
df1$Flq_acquired[df1$Flq_acquired == "qnrS1"] <- "qnrS"
df1$Flq_acquired[df1$Flq_acquired == ""] <- NA
df1$Flq_acquired[df1$Flq_acquired == "-"] <- NA

p10<-p9+new_scale_fill()
p11<- gheatmap(p10, df1[, "Flq_acquired", drop=FALSE],offset=0.00398, width=0.06, hjust = 1.0,colnames_position = "top",
               color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("qnr"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = c("#ffb3c1","#a4c3b2","#b5e48c","red"),
                    name = "qnr",
                    breaks =c("qnrS","qnrB","qnrA","qepA2"),na.value = "#edede9")

gg1
table(df1$K_locus_main)
df1$K_locus[df1$K_locus == "unknown (best match = KL31)"] <- "KL31"
df1$K_locus[df1$K_locus == "unknown (best match = KL147)"] <- "KL147"
df1$K_locus[df1$K_locus == "unknown (best match = KL49)"] <- "KL49"
df1$K_locus[df1$K_locus == "unknown (best match = KL74)"] <- "KL74"
df1$K_locus[df1$K_locus == "unknown (best match = KL145)"] <- "KL145"
df1$K_locus_main[df1$K_locus_main == ""] <- NA

nb.cols1 <- 22
mycolors1 <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols1)
p10<-gg1+new_scale_fill()
p11<- gheatmap(p10, df1[, "K_locus_main", drop=FALSE],offset=0.0013, width=0.06, hjust = 1.0,colnames_position = "top",
               color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("KL"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = mycolors1,na.value = "#edede9", name = "Capsule Locus")

p11
table(df1$O_locus_main)
df1$O_locus_main[df1$O_locus_main == ""] <- NA

nb.cols1 <- 7
mycolors1 <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols1)
p12<-p11+new_scale_fill()
p13<- gheatmap(p12, df1[, "O_locus_main", drop=FALSE],offset=0.00163, width=0.06, hjust = 1.0,colnames_position = "top",
               color = "white", font.size = 2.2, family="Arial Bold", colnames_angle = -90,custom_column_labels = c("O"))+
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5)+
  scale_fill_manual(values = mycolors1,na.value = "#edede9", name = "O-Locus")















p13
ggsave("test.png",width = 15.36, height = 8.14, dpi = 900)
