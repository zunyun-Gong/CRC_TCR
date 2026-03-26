############## Fig 3a ############

setwd("D:/TCR项目/result3")

####packages
library(tidyverse)
library(aplot)
library(pheatmap)
library(openxlsx)

library(ggplot2)
library(ggpubr)
library(EnvStats)


############ 筛选出显著相关的DMR
genebodydf <- read.csv("re_DMR_genebody_RNA_Pro_sig.Annot.txt",sep = "\t")
genebodydf_sig <- genebodydf[genebodydf$Sig!="no sig",]
genebodydf_sig_dmrs <- genebodydf_sig$DMR

promoterydf <- read.csv("re_DMR_promoter_RNA_Pro_sig.Annot.txt",sep = "\t")
promoterydf_sig <- promoterydf[promoterydf$Sig!="no sig",]
promoterydf_sig_dmrs <- promoterydf_sig$DMR
total_sigDMR <- c(genebodydf_sig_dmrs,promoterydf_sig_dmrs)

allsigDMR_df <- rbind(genebodydf_sig,promoterydf_sig)


########### 找出TCR组库相关的DMRs
msssum <- read.csv("MSSsum.withheader.txt",sep = "\t")
msssum <- msssum[msssum$DMR %in% total_sigDMR, ]

msssum_new <- msssum %>% group_by(cell) %>%
  mutate(adjP = p.adjust(pvalue,method = "BH")) %>%
  ungroup()

msssum_sig <- msssum_new %>%
  filter(adjP< 0.05) %>%  # 筛选 p值 < 0.05 的行
  filter(cell == "T_cell_diversity") %>%
  arrange(desc(abs(estimate))) %>%  # 按 change 列从大到小排序
  slice(1:20)

MSSDMRs <- msssum_sig[msssum_sig$cell=="T_cell_diversity" | msssum_sig$cell=="T_cell_clonality",]$DMR %>% unique


msisum <- read.csv("MSI.sum.withheader.txt",sep = "\t")
msisum <- msisum[msisum$DMR %in% total_sigDMR, ]
# msisum_sig <- msisum[msisum$pvalue<0.05 & abs(msisum$estimate)>0.3,]
msisum_sig <- msisum %>%
  filter(pvalue< 0.05) %>%  # 筛选 p值 < 0.05 的行
  filter(cell == "T_cell_diversity") %>%
  arrange(desc(abs(estimate))) %>%  # 按 change 列从大到小排序
  slice(1:20)

MSIDMRs <- msisum_sig[msisum_sig$cell=="T_cell_diversity" | msisum_sig$cell=="T_cell_clonality",]$DMR %>% unique

######### 确定是MSS 特异还是 MSI特异调控

alldmr <- c(MSIDMRs,MSSDMRs) %>% unique
type <- c()
for(i in alldmr){
  if(i %in% setdiff(MSIDMRs,MSSDMRs)){
    type <- append(type,"MSI_specific")
  }else if(i %in% setdiff(MSSDMRs,MSIDMRs)){
    type <- append(type,"MSS_specific")
  }else{
    type <- append(type,"Both")
  }
}

dmrmytype <- data.frame(DMR=alldmr,Type=type)


###############
MSSforplot <- merge(msssum,dmrmytype,by = "DMR")
MSSforplot$MSIstatus <- "MSS"

MSIforplot <- merge(msisum,dmrmytype,by = "DMR")
MSIforplot$MSIstatus <- "MSI"

forplot <- rbind(MSSforplot,MSIforplot)
DMRanno <- read.csv("group_DMR.status.promoter.genebody.bed",sep = "\t")
# DMR的注释
# DMRanno$gene <- apply(DMRanno,1,function(x){if(x[5]=="." & x[6]=="."){
#   return("Nogene")
#   } else{if(x[5]!="."){return(x[5])}else{return(x[6])}}
#   })

DMRanno_exp <- rbind(genebodydf,promoterydf)
DMRanno_exp <- DMRanno_exp[,c("DMR","Gene")]
DMRtype <- DMRanno[,c("DMR","status")]
DMRanno_exp <- merge(DMRanno_exp,DMRtype,by = "DMR")

# DMRanno <- DMRanno[,c("DMR",'status',"gene")]
forplot <- merge(forplot,DMRanno_exp,by = "DMR")

############ 将p值变为离散的
forplot<- forplot %>% mutate(sig = case_when(pvalue<0.001 ~ "p<0.001",
                                             pvalue < 0.01 ~ "p<0.01",
                                             pvalue<0.05 ~ "p<0.05",
                                             TRUE ~ "ns"))

######### 加入一个DMR的富集类型
fisherdf <- read.csv("MSSvsMSI_specific_dmrs.txt",sep = "\t")
fisherdf$Gene <- gsub("\\.","_",fisherdf$Gene)
fisherdf$OR <- (fisherdf$MSI_1/(fisherdf$MSI_1+fisherdf$MSI_0))/(fisherdf$MSS_1/(fisherdf$MSS_1+fisherdf$MSS_0))
fisherdf <- fisherdf %>% mutate(enrichtype = case_when(
  pvalue<0.05 & OR>1 ~ "MSI_enriched",
  pvalue<0.05 & OR<1 ~ "MSS_enriched",
  TRUE ~ "no_difference"
))
fisherdf <- fisherdf[,c("Gene","enrichtype")]
colnames(fisherdf) <- c("DMR","enrichtype")
forplot <- merge(forplot,fisherdf,by = "DMR")

###### 调整顺序
forplot$MSIstatus <- factor(forplot$MSIstatus,levels = c('MSS','MSI'))
forplot$Type <- factor(forplot$Type,levels = c("MSI_specific","MSS_specific","Both"))
forplot$sig <- factor(forplot$sig,levels = c("ns","p<0.05","p<0.01","p<0.001"))
forplot$cell <- factor(forplot$cell,levels = c("T_cell_diversity","T_cell_clonality","T-cells","CD8_T_cells","CD4+_T_cell","Exhausted_CD8","Th1_cells","Treg","Total_TIL"))
### 更改DMR名称
forplot <- forplot %>%
  mutate(DMR = paste0(DMR, "(", Gene, ")"))

dmrorder <- unique(forplot[order(forplot$Type,forplot$status),]$DMR)

forplot$DMR <- factor(forplot$DMR,levels = dmrorder)



status <- forplot %>% 
  mutate(p="Status") %>%
  ggplot(aes(x=p,y=DMR,fill=status))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(fill = "Status")+
  scale_fill_manual(values = c("#A91FC7","#D8D4D9"))

group <- forplot %>% 
  mutate(p="Group") %>%
  ggplot(aes(x=p,y=DMR,fill=Type))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(margin = margin(r = -1)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(fill = "Group")+
  scale_fill_manual(values = c("#F71409","#274BF8","#F79107"))


enrich <- forplot %>% 
  mutate(p="Enrichtype") %>%
  ggplot(aes(x=p,y=DMR,fill=enrichtype))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(fill = "Status")+
  scale_fill_manual(values = c("#EE715E","#5ECAEE","#D8D4D9"))




p1 <- ggplot(forplot,aes(x=cell,y=DMR,size=sig),fill=estimate)+
  facet_wrap(~ MSIstatus, nrow = NULL, ncol = NULL, scales = "fixed")+
  theme_bw()+
  geom_point(aes(fill = estimate), shape = 21, color = "black",stroke=0.05)+
  scale_fill_gradientn(colours = c("#5E98EE",'#DFE3EA','#F92717'))+
  scale_size_manual(
    values = c(`p<0.001` = 4, `p<0.01` = 3, `p<0.05` = 2,`ns`=1),  # 手动指定每个类别的大小值
    name = "sig"  # 图例标题
  )+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p1

p2 <- p1 %>%
  insert_left(status, width = .05) %>%
  insert_left(enrich, width = .05) %>%
  insert_left(group, width = .05)
p2


