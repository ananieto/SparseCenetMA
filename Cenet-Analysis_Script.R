#Packages: 

library(scatterplot3d)
library(ggplot2)
library(ggalt)
library(ggforce)
library(ggcorrplot)
library(ggthemes)


library(devtools)

devtools::install_github("HerveAbdi/data4PCCAR")
library(data4PCCAR)

devtools::install_version("rCUR", version = "1.3")
library(rCUR)

##################################################################################
#################                  EXAMPLE 1               #######################

#Set working directory -> dataset
#Read osiq data (example 1)

imagery.data<-read.table("osiq.txt", sep=" ")
head(imagery.data)
class(imagery.data)

#Data
dim(imagery.data)
J=dim(imagery.data)[2]
colnames(imagery.data)

#Reorder: 
imagery.data<-imagery.data[,c(1,2,3,5,6,9,11,13,14,18,20,
                              23,24,27,29,4,7,8,10,
                              12,15,16,17,19,21,22,25,26,28,30)]

a=0.5 #Alpha parameter

#HJ-Biplot: 
hj.res1<-Biplot.enet(imagery.data,Q=2,tau.u=sqrt(nrow(imagery.data)),
                           tau.v=sqrt(ncol(imagery.data)),init.transf = 4,
                           #alpha.v=0.5,
                           plot.axis=c(1,2))
hj.res1$biplot.plot
hj.res1$corrplot

#CenetHJ-Biplot: 
set.seed(2)
hj.res2<-Biplot.enet(imagery.data,Q=2,tau.u=sqrt(nrow(imagery.data)),
                     tau.v=((1-a)*sqrt(J)+a)*(3/4),init.transf = 4,
                     alpha.v=0.5,
                     plot.axis=c(1,2))
hj.res2$biplot.plot
hj.res2$corrplot


##################################################################################
#################                  EXAMPLE 2               #######################


colors<-c(rep("darkolivegreen",71),rep("firebrick2",74),rep("dodgerblue1",71))
clasification<-c(rep("Control", 71), rep("CLL", 74), rep("ALL", 71))

#Read "data216.2000" (variable selection: CUR Decomposition):  

  # cur.res<-rCUR::CUR(data216.global,c=(dim(data216.global)[2]-1),k=2,method="top.scores",weighted = T)
  # plotLeverage(cur.res)
  # variables.cur=2000
  # data216.2000<-data216.global[,cur.res@C.index[1:variables.cur]]

    ## install.packages("dCUR")
    ## library(dCUR) #Another function to execute the CUR decomposition is CUR (method=sample_cur) from the dCUR library.
  
##Set working directory -> dataset
##Read data216_2000 (example 2)
#Data
  
  data216.2000<-read.csv("data216_2000.csv", header=T)
  row.names(data216.2000)<-data216.2000[,1]
  data216.2000<-data216.2000[,-1]
  
  dim(data216.2000) #216 observations and 2000 variables
  colnames(data216.2000)

#Read info_samples
library(readxl)
info_216samples <- read_excel("info_216samples.xlsx")

#PCA
pca.res<-prcomp(data216.2000)
  (pca.res$sdev[1])^2/sum(pca.res$sdev^2)
  (pca.res$sdev[2])^2/sum(pca.res$sdev^2)

#SCORES
scores.pca<-data.frame(as.data.frame(pca.res$x[,1:2]), clasification)
ggplot(scores.pca, aes(x = PC1, y = PC2,
                     shape = clasification, 
                     color = clasification, fill=clasification)) +
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank())+
    geom_encircle(expand=0.02,alpha=0.2)+
  scale_fill_manual(values = c("#0E93DC","#F03D5D", "#5AB447"))+
  scale_color_manual(values = c("#0E93DC","#F03D5D", "#5AB447"))+
  geom_point(size=3, alpha = 0.8) 


#Loadings: 
loadings<-data.frame(var=colnames(data216.2000),as.data.frame(pca.res$rotation[,1:2]))

    ggplot(loadings, aes(x=var, y=PC1)) +
      geom_segment( aes(x=var, xend=var, y=0, yend=PC1), color="#B4B4B4") +
      geom_point( color="#9B6AD7", size=0.75) +
      xlab("Gene probes") +
      ylab("PC1 loadings")+
      geom_hline(yintercept=0, colour="gray")+
      theme_light() +
      theme(
        panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "#555555", fill = NA),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),#Para quitar las rayitas del eje x con las labels
        axis.text.x = element_blank()) 
    
    
    ggplot(loadings, aes(x=var, y=PC2)) +
      geom_segment( aes(x=var, xend=var, y=0, yend=PC2), color="#B4B4B4") +
      geom_point(color="#9B6AD7", size=0.75) +
      xlab("Gene probes") +
      ylab("PC2 loadings")+
      geom_hline(yintercept=0, colour="gray")+
      theme_light() +
      theme(
        panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "#555555", fill = NA),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 

#CenetPCA. Tau selection: 
data216.scale<-init.transformation(data216.2000,3)
#tau5<-bic.tau_enet(data216.scale,Q=2,alpha.v=0.5, ntau=10)
#tau5$tau
# tau5$BIC.plot
tauv=5.857929

#J=2000
#((1-0.5)*sqrt(J)+0.5)

#CenetPCA
set.seed(0)
pca.enet.res<-pca.enet(data216.scale,Q=2,tau.u=sqrt(216),
                    tau.v=tauv,alpha.v = 0.5)
pca.enet.res$Variance

#Scores:
  scores.enet<-as.data.frame(pca.enet.res$Scores)

  ggplot(scores.enet, aes(x = PC1, y = PC2,
                     shape = clasification, 
                     color = clasification, fill=clasification)) +
    theme_bw()+
    geom_encircle(expand=0.025,  alpha=0.2)+
    scale_fill_manual(values = c("#0E93DC","#F03D5D", "#5AB447"))+
    xlab("CenetPC1")+ylab("CenetPC2")+
    geom_point(size=3, alpha = 0.8)+
    scale_color_manual(values = c("#0E93DC","#F03D5D", "#5AB447"))


#Loading plots: 
loadings.enet<-data.frame(var=colnames(data216.2000),as.data.frame(pca.enet.res$Loadings))

ggplot(filter(loadings.enet, PC1!=0), aes(x=var, y=PC1)) +
  geom_segment( aes(x=var, xend=var, y=0, yend=PC1), color="#B4B4B4") +
  geom_point( color="#9B6AD7", size=1.5) +
  xlab("Gene probes") +
  ylab("CenetPC1 loadings")+
  geom_hline(yintercept=0, colour="gray")+
  theme_light() +
  theme(
    panel.grid.major  = element_line(color = "white"),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(color = "#555555", fill = NA),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) 


ggplot(filter(loadings.enet, PC2!=0), aes(x=var, y=PC2)) +
  geom_segment( aes(x=var, xend=var, y=0, yend=PC2), color="#B4B4B4") +
  geom_point(color="#9B6AD7", size=1.5) +
  xlab("Gene probes") +
  ylab("CenetPC2 loadings")+
  geom_hline(yintercept=0, colour="gray")+
  theme_light() +
  theme(
    panel.grid.major  = element_line(color = "white"),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(color = "#555555", fill = NA),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) 


#__________________________________________________
#Data: 216 observations, 205 genes
ind<-intersect(which(abs(loadings.enet$PC1)<0.01),which(abs(loadings.enet$PC2)<0.01))
data216.205<-data216.2000[,-ind]

#SymbolGene

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

library(affy)
library(Biobase)
library(hgu133plus2.db)
library(AnnotationDbi)

probe.gene<-select(hgu133plus2.db, keys=colnames(data216.205),
                   keytype = "PROBEID", 
                   columns=c("SYMBOL"))
names.match<-probe.gene[!duplicated(probe.gene$PROBEID),]
colnames(data216.205)<-names.match$SYMBOL
colnames(data216.205)<-make.names(colnames(data216.205), unique=TRUE)


#CenetHJ-BIPLOT
dim(data216.205) #216 observations; 205 variables
a=0.5
J=dim(data216.205)[2]


set.seed(0)
hj.res<-Biplot.enet(data216.205,Q=2,tau.u=sqrt(216),tau.v=((1-a)*sqrt(J)+a)*(1/3),init.transf = 3,alpha.v=0.5,
                    plot.axis=c(1,2))
hj.res$biplot.plot
class.tumor<-info_216samples$General_Group

i=1
j=2

df.u<-as.data.frame(hj.res$Row.coordinates)
df.v<-as.data.frame(hj.res$Col.coordinates)
row.names(df.u)<-row.names(hj.res$Row.coordinates)
colnames(df.u)<-colnames(hj.res$Row.coordinates)
row.names(df.v)<-row.names(hj.res$Col.coordinates)
colnames(df.v)<-colnames(hj.res$Col.coordinates)
df.v<-df.v[rowSums(df.v[,c(1,2)])!=0,]

axis.labels<-paste(colnames(df.u),sprintf("(%0.2f%%)", hj.res$Variance))
shape.216<-as.numeric(as.factor(clasification))
shape.216[which(shape.216==3)]<-15
shape.216[which(shape.216==2)]<-17
shape.216[which(shape.216==1)]<-16

col.216<-as.character(clasification)
col.216[which(col.216=="Control")]<-"#5AB447"
col.216[which(col.216=="CLL")]<-"#F03D5D"
col.216[which(col.216=="ALL")]<-"#0E93DC"

ggplot()+
  geom_point(data = df.u, aes(x=df.u[,i], y=df.u[,j]), color=col.216, shape=shape.216, size=2)+
  geom_segment(data = df.v, aes(x=0,y=0,xend=df.v[,i],yend=df.v[,j]),
             arrow=arrow(length=unit(2/3,"picas"), type="closed"),size=0.75)+ 
  theme_bw()+
  geom_hline(yintercept=0, color="gray", size=0.5)+
  geom_vline(xintercept=0, color="gray", size=0.5)+
  geom_text_repel(data=df.v,aes(label=row.names(df.v),x=df.v[,i],y=df.v[,j],fontface=4),hjust=2,
                   color="black",size=2.5) +#2=bold, 3=italic
  labs(x = axis.labels[i], y=axis.labels[j])




