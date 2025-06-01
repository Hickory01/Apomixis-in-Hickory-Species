#!/usr/bin/Rscript
library(ggplot2)
library(stringr)
library(tidyr)
library(stringr)
args<-commandArgs(TRUE)
data=read.table(args[1],header=T,check.names = FALSE)
data<-data[str_detect(data$CHROM,args[2]),]
for (i in levels(factor(data$CHROM))){
  data$LEN[data$CHROM==i]<-max(data$RANGE[data$CHROM==i])
}
all<-gather(data,key="SPECIES",value="Het_Dep",colnames(data)[3:(length(colnames(data))-1)])
all[,c('Het',"Dep")]<-str_split_fixed(all$Het_Dep,":",2)
all$Het<-as.numeric(all$Het)
new_data<-na.omit(all)
data1<-new_data[new_data$Het<=0.001,]
data2<-new_data[new_data$Het>=0.02,]
  ggplot(new_data)+
  geom_rect(data=data1,aes(xmin=RANGE,xmax=RANGE+49999,ymin=0,ymax=1),fill="DodgerBlue")+
  geom_rect(data=data2,aes(xmin=RANGE,xmax=RANGE+49999,ymin=0,ymax=1),fill="white")+
  geom_rect(aes(xmin = 0, xmax = LEN, ymin = 0, ymax = 1),color="grey81",fill="white",alpha=0)+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
		  panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_text(size=13),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=12),
          axis.ticks.y = element_blank(),
          plot.tag.position=c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.position = "right",
          plot.margin = margin(t = 10,  # 顶部边缘距离
                               r = 15,  # 右边边缘距
                               b = 20,  # 底部边缘距离
                               l = 10), # 左边边缘距离
          strip.text.y= element_text(size=8,angle=0,vjust=0.5,hjust=0),
          strip.background = element_blank(),
          panel.spacing = unit(0.05,"cm"),
          plot.tag=element_text(hjust = 0))+
  facet_grid(SPECIES~.)+
  scale_x_continuous(expand =c(0,0),breaks = seq(0,100000000,5000000),labels=seq(0,100,5),name=paste0(args[2]," (Mb)"))
ggsave(paste0(args[2],".LOH.png"),width = 10, height = 10)
