##find out common gene..  

library(readr)
library(pheatmap)
library(ggplot2)
getwd()
setwd("C:/Users/user/gsea_home/output/kegg_invivo")

library(readr)

install.packages("tidyverse")
library(tidyverse)


#initial parameter for for loop


gsea_report_up = read_tsv("gsea_report_for_CD4+_1.tsv")
gsea_report_down = read_tsv("gsea_report_for_DP_(CD4).tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)
save(gsea_report, file='gsea_normal_report_2021.05.15.rda')

plot_data = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data)[2] = 'NES_CD4_DP' 
colnames(plot_data)[3] = 'adj_p_CD4_DP' 
plot_data['adjP_thershold'] = NA
plot_data$adjP_thershold[plot_data$adj_p < 0.05] = 'black'
plot_data = plot_data[order(plot_data$NES),]
save(plot_data, file='gsea_DP_CD4_2021.05.15.rda')

####CD8_DP
gsea_report_up = read_tsv("gsea_report_for_CD8+_1.tsv")
gsea_report_down = read_tsv("gsea_report_for_DP_3_(CD8).tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)

plot_data2 = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data2)[2] = 'NES_CD8_DP' 
colnames(plot_data2)[3] = 'adj_p_CD8_DP' 
plot_data2['adjP_thershold'] = NA
plot_data2$adjP_thershold[plot_data2$adj_p < 0.05] = 'black'
plot_data2 = plot_data2[order(plot_data2$NES),]
save(plot_data2, file='gsea_DP_CD8_2021.05.15.rda')
####DN4_DP
gsea_report_up = read_tsv("gsea_report_for_DN4_2.tsv")
gsea_report_down = read_tsv("gsea_report_for_DP_1.tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)

plot_data3 = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data3)[2] = 'NES_DN4_DP' 
colnames(plot_data3)[3] = 'adj_p_DN4_DP' 
plot_data3['adjP_thershold'] = NA
plot_data3$adjP_thershold[plot_data3$adj_p < 0.05] = 'black'
plot_data3= plot_data3[order(plot_data3$NES),]
save(plot_data3, file='gsea_DN4_DP_2021.05.15.rda')
####DN3_DN4
gsea_report_up = read_tsv("gsea_report_for_DN4_1.tsv")
gsea_report_down = read_tsv("gsea_report_for_DN3_2.tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)

plot_data4 = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data4)[2] = 'NES_DN3_DN4' 
colnames(plot_data4)[3] = 'adj_p_DN3_DN4' 
plot_data4['adjP_thershold'] = NA
plot_data4$adjP_thershold[plot_data4$adj_p < 0.05] = 'black'
plot_data4 = plot_data4[order(plot_data4$NES),]
save(plot_data4, file='gsea_DN3_DN4_2021.05.15.rda')
####DN2_DN3
gsea_report_up = read_tsv("gsea_report_for_DN3_1.tsv")
gsea_report_down = read_tsv("gsea_report_for_DN2_2.tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)

plot_data5 = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data5)[2] = 'NES_DN2_DN3' 
colnames(plot_data5)[3] = 'adj_p_DN2_DN3' 
plot_data5['adjP_thershold'] = NA
plot_data5$adjP_thershold[plot_data5$adj_p < 0.05] = 'black'
plot_data5 = plot_data5[order(plot_data5$NES),]
save(plot_data5, file='gsea_DN2_DN3_2021.05.15.rda')
####DN1_DN2
gsea_report_up = read_tsv("gsea_report_for_DN2_1.tsv")
gsea_report_down = read_tsv("gsea_report_for_DN1_1.tsv")
gsea_report = rbind(gsea_report_up, gsea_report_down)
gsea_report = gsea_report[,-c(12)]

name_plot = gsub("KEGG_","", gsea_report$NAME) #to remove 'hallmark_' in name
gsea_report = cbind(name_plot, gsea_report)
str(gsea_report) #check variable type whether numberic or not
rm(gsea_report_up, gsea_report_down, name_plot)

plot_data6 = select(gsea_report, name_plot, NES, 'FDR q-val')
colnames(plot_data6)[2] = 'NES_DN1_DN2' 
colnames(plot_data6)[3] = 'adj_p_DN1_DN2' 
plot_data6['adjP_thershold'] = NA
plot_data6$adjP_thershold[plot_data6$adj_p < 0.05] = 'black'
plot_data6 = plot_data6[order(plot_data6$NES),]
save(plot_data6, file='gsea_DN1_DN2_2021.05.15.rda')

###ffff
df=load("gsea_DN1_DN2_2021.05.15.rda")
df=plot_data6
df2=load("gsea_DN2_DN3_2021.05.15.rda")
df2=plot_data5

df3=load("gsea_DN3_DN4_2021.05.15.rda")
df3=plot_data4
df4=load("gsea_DN4_DP_2021.05.15.rda")
df4=plot_data3
df_CD4=load("gsea_DP_CD4_2021.05.15.rda")
df_CD4=plot_data
df_CD8=load("gsea_DP_CD8_2021.05.15.rda")
df_CD8=plot_data2

df=na.omit(df)
df2=na.omit(df2)
df3=na.omit(df3)
df4=na.omit(df4)
df_CD4=na.omit(df_CD4)
df_CD8=na.omit(df_CD8)
df[df$adjP_thershold=='black',]a1a2a3<-data.frame(rbind(df$nameplot,df2$nameplot,df3$nameplot,df4[1],df_CD4[1],df_CD8[1]))


#unique(c(df$name_plot,df2$name_plot,df3$name_plot,df4$name_plot,df_CD4$name_plot,df_CD8$name_plot))
merge(df$name_plot,df2$name_plot,df3$name_plot,df4$name_plot,df_CD4$name_plot,df_CD8$name_plot)


df_merge=left_join(x=df, y=df3, by=c('name_plot'="name_plot"))

df_merge=left_join(x=df_merge, y=df4, by=c('name_plot'="name_plot"))

df_merge=left_join(x=df_merge, y=df_CD4, by=c('name_plot'="name_plot"))

df_merge=left_join(x=df_merge, y=df_CD8, by=c('name_plot'="name_plot"))

df_merge$adjP_thershold.x[df_merge$adjP_thershold.x=='black']=1

df_merge$adjP_thershold.y[df_merge$adjP_thershold.y=='black']=1

df_merge$adjP_thershold.x.x[df_merge$adjP_thershold.x.x=='black']=1
df_merge$adjP_thershold.y.y[df_merge$adjP_thershold.y.y=='black']=1
df_merge$adjP_thershold[df_merge$adjP_thershold=='black']=1

save(df_merge, file='gsea_df_merge_2021.05.15.rda')
