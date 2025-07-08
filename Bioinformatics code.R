chrX$MeanHealthy<-lapply(1:nrow(chrX),function(i)
  mean(as.matrix(chrX[i,9:16]),na.rm=TRUE))

chrX$MeanCLL<-lapply(1:nrow(chrX),function(i)
  mean(as.matrix(chrX[i,17:104]),na.rm=TRUE))

chrX$Diff<-((as.numeric(chrX$MeanCLL))-(as.numeric(chrX$MeanHealthy)))

chrX$Pvalues<-lapply(1:nrow(chrX),function(i) #Mann-Whitney Test
  wilcox.test(as.matrix(chrX[i,9:16]),as.matrix(chrX[i,17:104]))$p.value)

chrX$FDR<-p.adjust(chrX$Pvalues,method="fdr")


length(chrX$Pvalues)#11120
sum(chrX$Pvalues<0.05) #3081 - plotted on volcano plot
sum(chrX$FDR<0.05) #1456  - heatmap
sum(chrX$Diff<0.02)#9,652 #gain of methylation
sum(chrX$Diff>0.02)#1,468  # loss of methylation 


chrX<-chrX[order(as.numeric(chrX$Diff), decreasing = TRUE),]
head(chrX)

#VOLCANO PLOT
library(gplots)
library(ggplot2)
chrX$Sig<-as.factor(chrX$FDR<0.05)
VolcanoPlot<-ggplot(chrX,aes(x=Diff,y=-log10(as.numeric(Pvalues)),colour=Sig))
VolcanoPlot+geom_point()
VolcanoPlot+geom_point(size=5.0)
VolcanoPlot+geom_point(size=5.0,alpha=0.5)
VolcanoPlot+geom_point(size=5.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))

VolcanoPlot+geom_point(size=5.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))+theme(legend.position = "none")

VolcanoPlot+geom_point(size=5.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))+
  theme(legend.position = "none",axis.title=element_text(size=16))

VolcanoPlot+geom_point(size=2.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))+
  theme(legend.position = "none",axis.title=element_text(size=16))

VolcanoPlot+geom_point(size=2.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))+
  theme(legend.position = "none",axis.title=element_text(size=16),
        axis.text=element_text(size=16))


VolcanoPlot+geom_point(size=2.0,alpha=0.5)+
  scale_color_manual(values=c("slategrey","firebrick2"))+
  theme(legend.position = "none",axis.title=element_text(size=16),
        axis.text=element_text(size=16))+
  labs(x="Methylation change(beta)",y="-log10(P-value)")

VolcanoPlot+geom_point(size=1.0,alpha=0.5)+
  scale_color_manual(values=c("black","red"))+
  theme(legend.position = "none",axis.title=element_text(size=16),
        axis.text=element_text(size=16))+
  labs(x="Methylation change(beta)",y="-log10(P-value)")+xlim(-1,1)




dev.off()
#HEATMAPS

library(gplots)
library(RColorBrewer)
HeatCol<-colorRampPalette(c("yellow","black","blue"))(50)

chrX<-chrX[order(as.numeric(chrX$FDR), decreasing = FALSE),]

heatmap.2(as.matrix(chrX[1:200,9:104]),main="Heatmap of X-linked CpG loci",trace="none",
          col = HeatCol,labRow = chrX$ProbeID,
          ColSideColors= c(rep("paleturquoise2",8),rep("tomato",47),
                           rep("darkred",4), rep("skyblue2",37)))
key = FALSE # Turn off default legend


legend("topright", legend = c("M-CLL", "U-CLL", "Unknown", "Control"), fill = c("tomato", "skyblue2", "darkred", "paleturquoise2"),
       title = "IGHV status")

colnames(chrX)


#In case of 'invalid graphics state'
dev.off()

NewDF<-subset(chrX,subset = chrX$FDR<0.05)

NewDF<-subset(chrX,subset=chrX$Diff>0.2)

NewDF<-subset(chrX,subset=chrX$Diff<0.2)

NewDF<-subset(chrX,subset=chrX$Diff<0.2|chrX$Diff>0.2)

NewDF<-subset(chrX,subset=chrX$FDR<0.05&chrX$Diff<0.2|chrX$Diff>0.2)              

write.csv(as.matrix(NewDF),file="~/NewDF.csv",row.names = FALSE)

#Kaplan-Meier
install.packages("survminer")
library(survminer)
library(survival)

sfit<-survfit(Surv(OS,Status)~rporcn,data=expression1)

ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("High ","Low "), legend.title=".",
           palette=c("firebrick1","gray0"),
           title="Overall Survival by PORCN expression",
           risk.table.height=.2,          
           xlab = "Time (Days)")

str(expression)

dev.off()     

#Boxplots

ggplot(data=Exp_Landau,aes(x=Group,y=BCAP31,fill=Group))+geom_boxplot()

Exp_Landau$Group=factor(Exp_Landau$Group,c("Healthy","CLL"))

ggplot(data=Exp_Landau,aes(x=Group,y=BCAP31,fill=Group))+geom_boxplot()
ggplot(data=Exp_Landau,aes(x=Group,y=BCAP31,fill=Group))+geom_boxplot()+
  labs(title="BCAP31",y="Methylation(beta)",x="Group")

ggplot(data=Exp_Landau,aes(x=Group,y=BCAP31,fill=Group))+
  labs(title="BCAP31",y="Methylation(beta)",x="Group")+
  theme(legend.position="none",
        plot.title=element_text(size=24,face="bold"),
        axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12))+
  ylim(0,1)


ggplot(data = Exp_Landau, aes(x = Group, y = BCAP31, fill = Group)) +
  labs(title = "BCAP31", y = "Methylation(beta)", x = "Group") +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12)) +
  ylim(0, 1) +
  scale_fill_manual(values = c("darkorange2", "dodgerblue3"))


ggplot(data = Exp_Landau, aes(x = Group, y = BCAP31, fill = Group)) +
  geom_boxplot() +labs(title = "BCAP31", y = "Expression", x = "Group") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12)) +
  ylim(9, 14) +
  scale_fill_manual(values = c("gray65", "indianred2"))



#Wilcox
wilcox.test(CD40LG ~ Group, data = Exp_Landau)


