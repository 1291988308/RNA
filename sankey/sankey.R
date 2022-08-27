#install.packages(c('ggalluvial', 'ggplot'))
rm(list = ls())
getwd()
dir()
library(ggalluvial)
library(xlsx)
df<- read.xlsx("ceRNAnetwork.xlsx", sheetName = "circ-mi-m", header = T,encoding = 'UTF-8')

#df <- read.csv("circRNA-miRNA-mRNA network.csv",header = T)
#Check whether the data meets requirements；
head(df)
df <- df[,c(1:4)]
is_alluvia_form(df,weight ="Freq")
# Convert to long data format；
df_lodes <- to_lodes_form(df,key ="x", value = "stratum", id = "alluvium",axes =1:3)
#Check whether the converted data meets the mapping requirements；
head(df_lodes,12)
is_lodes_form(df_lodes,key = "x",value = "stratum",id = "alluvium",weight ="Freq")
mycol3=colorRampPalette(c("#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#00abef","#64b036","#ffe743","#64b036","#00abef"))(302)
#mycol3 <-rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),3)
p<-ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                        fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0, knot.pos = 0.1) +
  geom_stratum(alpha = .9,color="grey20",width = 1/7) +
  geom_text(stat = "stratum", size =1.0,color="black") +
  scale_fill_manual(values = mycol3) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)
p
pdf(file = "ceRNAnetwork.pdf",width =8.5,height = 11)
print(p)
dev.off()


