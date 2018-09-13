#------------------------------------------------
# Author: Maria Zhivagui                        #
# Script name: SBS_Mutation_count_plotting.R    #
# Last update: September 10 2018                #
#------------------------------------------------

#Packages
library(reshape2)
library(dplyr)
library(stringr)
library(plyr)
library(reshape2)
library(reshape)
library(Rmisc)
library(scales)
library(easyGgplot2)
library(ggpubr)
library(cowplot)

#Get data
setwd("~/Documents/CMM Lab documents/Bioinformatics/Tumors assignment")

#combing the 4 files
e1<-read.table("TCGA_WES_sigProfiler_SBS_signatures_in_samples.txt", sep="\t", header=T)
#calculate WES / 51 per mg
DF1=ddply(e1,.(SBS1),transform,SBS1=SBS1/51)
DF2=ddply(DF1,.(SBS2),transform,SBS2=SBS2/51)
DF3=ddply(DF2,.(SBS3),transform,SBS3=SBS3/51)
DF4=ddply(DF3,.(SBS4),transform,SBS4=SBS4/51)
DF5=ddply(DF4,.(SBS5),transform,SBS5=SBS5/51)
DF6=ddply(DF5,.(SBS6),transform,SBS6=SBS6/51)
DF7=ddply(DF6,.(SBS7a),transform,SBS7a=SBS7a/51)
DF8=ddply(DF7,.(SBS7b),transform,SBS7b=SBS7b/51)
DF9=ddply(DF8,.(SBS7c),transform,SBS7c=SBS7c/51)
DF10=ddply(DF9,.(SBS7d),transform,SBS7d=SBS7d/51)
DF11=ddply(DF10,.(SBS8),transform,SBS8=SBS8/51)
DF12=ddply(DF11,.(SBS9),transform,SBS9=SBS9/51)
DF13=ddply(DF12,.(SBS10a),transform,SBS10a=SBS10a/51)
DF14=ddply(DF13,.(SBS10b),transform,SBS10b=SBS10b/51)
DF15=ddply(DF14,.(SBS11),transform,SBS11=SBS11/51)
DF16=ddply(DF15,.(SBS12),transform,SBS12=SBS12/51)
DF17=ddply(DF16,.(SBS13),transform,SBS13=SBS13/51)
DF18=ddply(DF17,.(SBS14),transform,SBS14=SBS14/51)
DF19=ddply(DF18,.(SBS15),transform,SBS15=SBS15/51)
DF20=ddply(DF19,.(SBS16),transform,SBS16=SBS16/51)
DF21=ddply(DF20,.(SBS17a),transform,SBS17a=SBS17a/51)
DF22=ddply(DF21,.(SBS17b),transform,SBS17b=SBS17b/51)
DF23=ddply(DF22,.(SBS18),transform,SBS18=SBS18/51)
DF24=ddply(DF23,.(SBS19),transform,SBS19=SBS19/51)
DF25=ddply(DF24,.(SBS20),transform,SBS20=SBS20/51)
DF26=ddply(DF25,.(SBS21),transform,SBS21=SBS21/51)
DF27=ddply(DF26,.(SBS22),transform,SBS22=SBS22/51)
DF28=ddply(DF27,.(SBS23),transform,SBS23=SBS23/51)
DF29=ddply(DF28,.(SBS24),transform,SBS24=SBS24/51)
DF30=ddply(DF29,.(SBS25),transform,SBS25=SBS25/51)
DF31=ddply(DF30,.(SBS26),transform,SBS26=SBS26/51)
DF32=ddply(DF31,.(SBS27),transform,SBS27=SBS27/51)
DF33=ddply(DF32,.(SBS28),transform,SBS28=SBS28/51)
DF34=ddply(DF33,.(SBS29),transform,SBS29=SBS29/51)
DF35=ddply(DF34,.(SBS30),transform,SBS30=SBS30/51)
DF36=ddply(DF35,.(SBS31),transform,SBS31=SBS31/51)
DF37=ddply(DF36,.(SBS32),transform,SBS32=SBS32/51)
DF38=ddply(DF37,.(SBS33),transform,SBS33=SBS33/51)
DF39=ddply(DF38,.(SBS34),transform,SBS34=SBS34/51)
DF40=ddply(DF39,.(SBS35),transform,SBS35=SBS35/51)
DF41=ddply(DF40,.(SBS36),transform,SBS36=SBS36/51)
DF42=ddply(DF41,.(SBS37),transform,SBS37=SBS37/51)
DF43=ddply(DF42,.(SBS38),transform,SBS38=SBS38/51)
DF44=ddply(DF43,.(SBS39),transform,SBS39=SBS39/51)
DF45=ddply(DF44,.(SBS40),transform,SBS40=SBS40/51)
DF46=ddply(DF45,.(SBS41),transform,SBS41=SBS41/51)
DF47=ddply(DF46,.(SBS42),transform,SBS42=SBS42/51)
DF48=ddply(DF47,.(SBS43),transform,SBS43=SBS43/51)
DF49=ddply(DF48,.(SBS44),transform,SBS44=SBS44/51)
DF50=ddply(DF49,.(SBS45),transform,SBS45=SBS45/51)
DF51=ddply(DF50,.(SBS46),transform,SBS46=SBS46/51)
DF52=ddply(DF51,.(SBS47),transform,SBS47=SBS47/51)
DF53=ddply(DF52,.(SBS48),transform,SBS48=SBS48/51)
DF54=ddply(DF53,.(SBS49),transform,SBS49=SBS49/51)
DF55=ddply(DF54,.(SBS50),transform,SBS50=SBS50/51)
DF56=ddply(DF55,.(SBS51),transform,SBS51=SBS51/51)
DF57=ddply(DF56,.(SBS52),transform,SBS52=SBS52/51)
DF58=ddply(DF57,.(SBS53),transform,SBS53=SBS53/51)
DF59=ddply(DF58,.(SBS54),transform,SBS54=SBS54/51)
DF60=ddply(DF59,.(SBS55),transform,SBS55=SBS55/51)
DF61=ddply(DF60,.(SBS56),transform,SBS56=SBS56/51)
DF62=ddply(DF61,.(SBS57),transform,SBS57=SBS57/51)
DF63=ddply(DF62,.(SBS58),transform,SBS58=SBS58/51)
DF64=ddply(DF63,.(SBS59),transform,SBS59=SBS59/51)
DF65=ddply(DF64,.(SBS60),transform,SBS60=SBS60/51)

e10 = DF65
#subset cosine similarity>0.9 and mutations count>100
AC1<- e10 %>% filter(Accuracy >= "0.85") #do with 0.85 before:0.9
tmp1= AC1[, c(1:2, 4:68)] #remove accuracy column 
AC11<- melt(tmp1)
AC11$value[AC11$value<0.9803922] <- 0 # 50 mut/mb cutoff
AC11

e2<-read.table("nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13 (1).txt", sep="\t", header=T)
#calculate WES / 51 per mg
DF1=ddply(e2,.(SBS1),transform,SBS1=SBS1/51)
DF2=ddply(DF1,.(SBS2),transform,SBS2=SBS2/51)
DF3=ddply(DF2,.(SBS3),transform,SBS3=SBS3/51)
DF4=ddply(DF3,.(SBS4),transform,SBS4=SBS4/51)
DF5=ddply(DF4,.(SBS5),transform,SBS5=SBS5/51)
DF6=ddply(DF5,.(SBS6),transform,SBS6=SBS6/51)
DF7=ddply(DF6,.(SBS7a),transform,SBS7a=SBS7a/51)
DF8=ddply(DF7,.(SBS7b),transform,SBS7b=SBS7b/51)
DF9=ddply(DF8,.(SBS7c),transform,SBS7c=SBS7c/51)
DF10=ddply(DF9,.(SBS7d),transform,SBS7d=SBS7d/51)
DF11=ddply(DF10,.(SBS8),transform,SBS8=SBS8/51)
DF12=ddply(DF11,.(SBS9),transform,SBS9=SBS9/51)
DF13=ddply(DF12,.(SBS10a),transform,SBS10a=SBS10a/51)
DF14=ddply(DF13,.(SBS10b),transform,SBS10b=SBS10b/51)
DF15=ddply(DF14,.(SBS11),transform,SBS11=SBS11/51)
DF16=ddply(DF15,.(SBS12),transform,SBS12=SBS12/51)
DF17=ddply(DF16,.(SBS13),transform,SBS13=SBS13/51)
DF18=ddply(DF17,.(SBS14),transform,SBS14=SBS14/51)
DF19=ddply(DF18,.(SBS15),transform,SBS15=SBS15/51)
DF20=ddply(DF19,.(SBS16),transform,SBS16=SBS16/51)
DF21=ddply(DF20,.(SBS17a),transform,SBS17a=SBS17a/51)
DF22=ddply(DF21,.(SBS17b),transform,SBS17b=SBS17b/51)
DF23=ddply(DF22,.(SBS18),transform,SBS18=SBS18/51)
DF24=ddply(DF23,.(SBS19),transform,SBS19=SBS19/51)
DF25=ddply(DF24,.(SBS20),transform,SBS20=SBS20/51)
DF26=ddply(DF25,.(SBS21),transform,SBS21=SBS21/51)
DF27=ddply(DF26,.(SBS22),transform,SBS22=SBS22/51)
DF28=ddply(DF27,.(SBS23),transform,SBS23=SBS23/51)
DF29=ddply(DF28,.(SBS24),transform,SBS24=SBS24/51)
DF30=ddply(DF29,.(SBS25),transform,SBS25=SBS25/51)
DF31=ddply(DF30,.(SBS26),transform,SBS26=SBS26/51)
DF32=ddply(DF31,.(SBS27),transform,SBS27=SBS27/51)
DF33=ddply(DF32,.(SBS28),transform,SBS28=SBS28/51)
DF34=ddply(DF33,.(SBS29),transform,SBS29=SBS29/51)
DF35=ddply(DF34,.(SBS30),transform,SBS30=SBS30/51)
DF36=ddply(DF35,.(SBS31),transform,SBS31=SBS31/51)
DF37=ddply(DF36,.(SBS32),transform,SBS32=SBS32/51)
DF38=ddply(DF37,.(SBS33),transform,SBS33=SBS33/51)
DF39=ddply(DF38,.(SBS34),transform,SBS34=SBS34/51)
DF40=ddply(DF39,.(SBS35),transform,SBS35=SBS35/51)
DF41=ddply(DF40,.(SBS36),transform,SBS36=SBS36/51)
DF42=ddply(DF41,.(SBS37),transform,SBS37=SBS37/51)
DF43=ddply(DF42,.(SBS38),transform,SBS38=SBS38/51)
DF44=ddply(DF43,.(SBS39),transform,SBS39=SBS39/51)
DF45=ddply(DF44,.(SBS40),transform,SBS40=SBS40/51)
DF46=ddply(DF45,.(SBS41),transform,SBS41=SBS41/51)
DF47=ddply(DF46,.(SBS42),transform,SBS42=SBS42/51)
DF48=ddply(DF47,.(SBS43),transform,SBS43=SBS43/51)
DF49=ddply(DF48,.(SBS44),transform,SBS44=SBS44/51)
DF50=ddply(DF49,.(SBS45),transform,SBS45=SBS45/51)
DF51=ddply(DF50,.(SBS46),transform,SBS46=SBS46/51)
DF52=ddply(DF51,.(SBS47),transform,SBS47=SBS47/51)
DF53=ddply(DF52,.(SBS48),transform,SBS48=SBS48/51)
DF54=ddply(DF53,.(SBS49),transform,SBS49=SBS49/51)
DF55=ddply(DF54,.(SBS50),transform,SBS50=SBS50/51)
DF56=ddply(DF55,.(SBS51),transform,SBS51=SBS51/51)
DF57=ddply(DF56,.(SBS52),transform,SBS52=SBS52/51)
DF58=ddply(DF57,.(SBS53),transform,SBS53=SBS53/51)
DF59=ddply(DF58,.(SBS54),transform,SBS54=SBS54/51)
DF60=ddply(DF59,.(SBS55),transform,SBS55=SBS55/51)
DF61=ddply(DF60,.(SBS56),transform,SBS56=SBS56/51)
DF62=ddply(DF61,.(SBS57),transform,SBS57=SBS57/51)
DF63=ddply(DF62,.(SBS58),transform,SBS58=SBS58/51)
DF64=ddply(DF63,.(SBS59),transform,SBS59=SBS59/51)
DF65=ddply(DF64,.(SBS60),transform,SBS60=SBS60/51)

e20 = DF65
#subset cosine similarity>0.9 and mutations count>100
AC2<- e20 %>% filter(Accuracy >= "0.85")
tmp2= AC2[, c(1:2, 4:68)]
AC22<- melt(tmp2)
AC22$value[AC22$value<0.9803922] <- 0

e3<-read.table("nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.txt", sep="\t", header=T)
DF1=ddply(e3,.(SBS1),transform,SBS1=SBS1/2781)
DF2=ddply(DF1,.(SBS2),transform,SBS2=SBS2/2781)
DF3=ddply(DF2,.(SBS3),transform,SBS3=SBS3/2781)
DF4=ddply(DF3,.(SBS4),transform,SBS4=SBS4/2781)
DF5=ddply(DF4,.(SBS5),transform,SBS5=SBS5/2781)
DF6=ddply(DF5,.(SBS6),transform,SBS6=SBS6/2781)
DF7=ddply(DF6,.(SBS7a),transform,SBS7a=SBS7a/2781)
DF8=ddply(DF7,.(SBS7b),transform,SBS7b=SBS7b/2781)
DF9=ddply(DF8,.(SBS7c),transform,SBS7c=SBS7c/2781)
DF10=ddply(DF9,.(SBS7d),transform,SBS7d=SBS7d/2781)
DF11=ddply(DF10,.(SBS8),transform,SBS8=SBS8/2781)
DF12=ddply(DF11,.(SBS9),transform,SBS9=SBS9/2781)
DF13=ddply(DF12,.(SBS10a),transform,SBS10a=SBS10a/2781)
DF14=ddply(DF13,.(SBS10b),transform,SBS10b=SBS10b/2781)
DF15=ddply(DF14,.(SBS11),transform,SBS11=SBS11/2781)
DF16=ddply(DF15,.(SBS12),transform,SBS12=SBS12/2781)
DF17=ddply(DF16,.(SBS13),transform,SBS13=SBS13/2781)
DF18=ddply(DF17,.(SBS14),transform,SBS14=SBS14/2781)
DF19=ddply(DF18,.(SBS15),transform,SBS15=SBS15/2781)
DF20=ddply(DF19,.(SBS16),transform,SBS16=SBS16/2781)
DF21=ddply(DF20,.(SBS17a),transform,SBS17a=SBS17a/2781)
DF22=ddply(DF21,.(SBS17b),transform,SBS17b=SBS17b/2781)
DF23=ddply(DF22,.(SBS18),transform,SBS18=SBS18/2781)
DF24=ddply(DF23,.(SBS19),transform,SBS19=SBS19/2781)
DF25=ddply(DF24,.(SBS20),transform,SBS20=SBS20/2781)
DF26=ddply(DF25,.(SBS21),transform,SBS21=SBS21/2781)
DF27=ddply(DF26,.(SBS22),transform,SBS22=SBS22/2781)
DF28=ddply(DF27,.(SBS23),transform,SBS23=SBS23/2781)
DF29=ddply(DF28,.(SBS24),transform,SBS24=SBS24/2781)
DF30=ddply(DF29,.(SBS25),transform,SBS25=SBS25/2781)
DF31=ddply(DF30,.(SBS26),transform,SBS26=SBS26/2781)
DF32=ddply(DF31,.(SBS27),transform,SBS27=SBS27/2781)
DF33=ddply(DF32,.(SBS28),transform,SBS28=SBS28/2781)
DF34=ddply(DF33,.(SBS29),transform,SBS29=SBS29/2781)
DF35=ddply(DF34,.(SBS30),transform,SBS30=SBS30/2781)
DF36=ddply(DF35,.(SBS31),transform,SBS31=SBS31/2781)
DF37=ddply(DF36,.(SBS32),transform,SBS32=SBS32/2781)
DF38=ddply(DF37,.(SBS33),transform,SBS33=SBS33/2781)
DF39=ddply(DF38,.(SBS34),transform,SBS34=SBS34/2781)
DF40=ddply(DF39,.(SBS35),transform,SBS35=SBS35/2781)
DF41=ddply(DF40,.(SBS36),transform,SBS36=SBS36/2781)
DF42=ddply(DF41,.(SBS37),transform,SBS37=SBS37/2781)
DF43=ddply(DF42,.(SBS38),transform,SBS38=SBS38/2781)
DF44=ddply(DF43,.(SBS39),transform,SBS39=SBS39/2781)
DF45=ddply(DF44,.(SBS40),transform,SBS40=SBS40/2781)
DF46=ddply(DF45,.(SBS41),transform,SBS41=SBS41/2781)
DF47=ddply(DF46,.(SBS42),transform,SBS42=SBS42/2781)
DF48=ddply(DF47,.(SBS43),transform,SBS43=SBS43/2781)
DF49=ddply(DF48,.(SBS44),transform,SBS44=SBS44/2781)
DF50=ddply(DF49,.(SBS45),transform,SBS45=SBS45/2781)
DF51=ddply(DF50,.(SBS46),transform,SBS46=SBS46/2781)
DF52=ddply(DF51,.(SBS47),transform,SBS47=SBS47/2781)
DF53=ddply(DF52,.(SBS48),transform,SBS48=SBS48/2781)
DF54=ddply(DF53,.(SBS49),transform,SBS49=SBS49/2781)
DF55=ddply(DF54,.(SBS50),transform,SBS50=SBS50/2781)
DF56=ddply(DF55,.(SBS51),transform,SBS51=SBS51/2781)
DF57=ddply(DF56,.(SBS52),transform,SBS52=SBS52/2781)
DF58=ddply(DF57,.(SBS53),transform,SBS53=SBS53/2781)
DF59=ddply(DF58,.(SBS54),transform,SBS54=SBS54/2781)
DF60=ddply(DF59,.(SBS55),transform,SBS55=SBS55/2781)
DF61=ddply(DF60,.(SBS56),transform,SBS56=SBS56/2781)
DF62=ddply(DF61,.(SBS57),transform,SBS57=SBS57/2781)
DF63=ddply(DF62,.(SBS58),transform,SBS58=SBS58/2781)
DF64=ddply(DF63,.(SBS59),transform,SBS59=SBS59/2781)
DF65=ddply(DF64,.(SBS60),transform,SBS60=SBS60/2781)

e30 = DF65
#subset cosine similarity>0.9 and mutations count>100
AC3<- e30 %>% filter(Accuracy >= "0.85")
tmp3= AC3[, c(1:2, 4:68)]
AC33<- melt(tmp3)
AC33$value[AC11$value<0.01797914] <- 0

e4<-read.table("PCAWG_sigProfiler_SBS_signatures_in_samples.txt", sep="\t", header=T)
DF1=ddply(e4,.(SBS1),transform,SBS1=SBS1/2781)
DF2=ddply(DF1,.(SBS2),transform,SBS2=SBS2/2781)
DF3=ddply(DF2,.(SBS3),transform,SBS3=SBS3/2781)
DF4=ddply(DF3,.(SBS4),transform,SBS4=SBS4/2781)
DF5=ddply(DF4,.(SBS5),transform,SBS5=SBS5/2781)
DF6=ddply(DF5,.(SBS6),transform,SBS6=SBS6/2781)
DF7=ddply(DF6,.(SBS7a),transform,SBS7a=SBS7a/2781)
DF8=ddply(DF7,.(SBS7b),transform,SBS7b=SBS7b/2781)
DF9=ddply(DF8,.(SBS7c),transform,SBS7c=SBS7c/2781)
DF10=ddply(DF9,.(SBS7d),transform,SBS7d=SBS7d/2781)
DF11=ddply(DF10,.(SBS8),transform,SBS8=SBS8/2781)
DF12=ddply(DF11,.(SBS9),transform,SBS9=SBS9/2781)
DF13=ddply(DF12,.(SBS10a),transform,SBS10a=SBS10a/2781)
DF14=ddply(DF13,.(SBS10b),transform,SBS10b=SBS10b/2781)
DF15=ddply(DF14,.(SBS11),transform,SBS11=SBS11/2781)
DF16=ddply(DF15,.(SBS12),transform,SBS12=SBS12/2781)
DF17=ddply(DF16,.(SBS13),transform,SBS13=SBS13/2781)
DF18=ddply(DF17,.(SBS14),transform,SBS14=SBS14/2781)
DF19=ddply(DF18,.(SBS15),transform,SBS15=SBS15/2781)
DF20=ddply(DF19,.(SBS16),transform,SBS16=SBS16/2781)
DF21=ddply(DF20,.(SBS17a),transform,SBS17a=SBS17a/2781)
DF22=ddply(DF21,.(SBS17b),transform,SBS17b=SBS17b/2781)
DF23=ddply(DF22,.(SBS18),transform,SBS18=SBS18/2781)
DF24=ddply(DF23,.(SBS19),transform,SBS19=SBS19/2781)
DF25=ddply(DF24,.(SBS20),transform,SBS20=SBS20/2781)
DF26=ddply(DF25,.(SBS21),transform,SBS21=SBS21/2781)
DF27=ddply(DF26,.(SBS22),transform,SBS22=SBS22/2781)
DF28=ddply(DF27,.(SBS23),transform,SBS23=SBS23/2781)
DF29=ddply(DF28,.(SBS24),transform,SBS24=SBS24/2781)
DF30=ddply(DF29,.(SBS25),transform,SBS25=SBS25/2781)
DF31=ddply(DF30,.(SBS26),transform,SBS26=SBS26/2781)
DF32=ddply(DF31,.(SBS27),transform,SBS27=SBS27/2781)
DF33=ddply(DF32,.(SBS28),transform,SBS28=SBS28/2781)
DF34=ddply(DF33,.(SBS29),transform,SBS29=SBS29/2781)
DF35=ddply(DF34,.(SBS30),transform,SBS30=SBS30/2781)
DF36=ddply(DF35,.(SBS31),transform,SBS31=SBS31/2781)
DF37=ddply(DF36,.(SBS32),transform,SBS32=SBS32/2781)
DF38=ddply(DF37,.(SBS33),transform,SBS33=SBS33/2781)
DF39=ddply(DF38,.(SBS34),transform,SBS34=SBS34/2781)
DF40=ddply(DF39,.(SBS35),transform,SBS35=SBS35/2781)
DF41=ddply(DF40,.(SBS36),transform,SBS36=SBS36/2781)
DF42=ddply(DF41,.(SBS37),transform,SBS37=SBS37/2781)
DF43=ddply(DF42,.(SBS38),transform,SBS38=SBS38/2781)
DF44=ddply(DF43,.(SBS39),transform,SBS39=SBS39/2781)
DF45=ddply(DF44,.(SBS40),transform,SBS40=SBS40/2781)
DF46=ddply(DF45,.(SBS41),transform,SBS41=SBS41/2781)
DF47=ddply(DF46,.(SBS42),transform,SBS42=SBS42/2781)
DF48=ddply(DF47,.(SBS43),transform,SBS43=SBS43/2781)
DF49=ddply(DF48,.(SBS44),transform,SBS44=SBS44/2781)
DF50=ddply(DF49,.(SBS45),transform,SBS45=SBS45/2781)
DF51=ddply(DF50,.(SBS46),transform,SBS46=SBS46/2781)
DF52=ddply(DF51,.(SBS47),transform,SBS47=SBS47/2781)
DF53=ddply(DF52,.(SBS48),transform,SBS48=SBS48/2781)
DF54=ddply(DF53,.(SBS49),transform,SBS49=SBS49/2781)
DF55=ddply(DF54,.(SBS50),transform,SBS50=SBS50/2781)
DF56=ddply(DF55,.(SBS51),transform,SBS51=SBS51/2781)
DF57=ddply(DF56,.(SBS52),transform,SBS52=SBS52/2781)
DF58=ddply(DF57,.(SBS53),transform,SBS53=SBS53/2781)
DF59=ddply(DF58,.(SBS54),transform,SBS54=SBS54/2781)
DF60=ddply(DF59,.(SBS55),transform,SBS55=SBS55/2781)
DF61=ddply(DF60,.(SBS56),transform,SBS56=SBS56/2781)
DF62=ddply(DF61,.(SBS57),transform,SBS57=SBS57/2781)
DF63=ddply(DF62,.(SBS58),transform,SBS58=SBS58/2781)
DF64=ddply(DF63,.(SBS59),transform,SBS59=SBS59/2781)
DF65=ddply(DF64,.(SBS60),transform,SBS60=SBS60/2781)

e40 = DF65

AC4<- e40 %>% filter(Accuracy >= "0.85")
tmp4= AC4[, c(1:2, 4:68)]
AC44<- melt(tmp4)
AC44$value[AC11$value<0.01797914] <- 0

key1<-unlist(apply(AC11, 1, function(x){ paste(x[1], x[2], x[3], x[4], x[5], sep="") } ))
key2<-unlist(apply(AC22, 1, function(x){ paste(x[1], x[2], x[3], x[4], x[5], sep="") } ))
key3<-unlist(apply(AC33, 1, function(x){ paste(x[1], x[2], x[3], x[4], x[5], sep="") } ))
key4<-unlist(apply(AC44, 1, function(x){ paste(x[1], x[2], x[3], x[4], x[5], sep="") } ))

exome1<-duplicated(c(key1,key2))
genome1<-duplicated(c(key3, key4))

tmpe<-rbind(AC11, AC22)
tmpg<-rbind(AC33, AC44)

exome<-tmpe[!exome1,]
genome<-tmpg[!genome1,]

#get total number of samples before filtering
totexome=rbind(e1, e2)
totgenome=rbind(e3, e4)

#Breast
Breast_samples_exome=dplyr::filter(totexome, grepl('Breast', Cancer.Types)) #1160
Breast_samples_genome=dplyr::filter(totgenome, grepl('Breast', Cancer.Types))#698


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#filter dara based on EXome cancer types
# adrenal <- dplyr::filter(loulou, grepl('Adrenal', Cancer.Types))
# AML <- dplyr::filter(loulou, grepl('AML', Cancer.Types))
# biliary_exome <- dplyr::filter(exome, grepl('Biliary', Cancer.Types))
breast_exome <- dplyr::filter(exome, grepl('Breast', Cancer.Types)) #27040
breast_genome <- dplyr::filter(genome, grepl('Breast', Cancer.Types)) #45370
# cervix <- dplyr::filter(loulou, grepl('Cervix', Cancer.Types))
# CNS <- dplyr::filter(loulou, grepl('CNS', Cancer.Types))
# colorect <- dplyr::filter(loulou, grepl('ColoRect', Cancer.Types))
# DLBC <- dplyr::filter(loulou, grepl('DLBC', Cancer.Types))
# eso <- dplyr::filter(loulou, grepl('Eso', Cancer.Types))
# eye <- dplyr::filter(loulou, grepl('Eye', Cancer.Types))
# head <- dplyr::filter(loulou, grepl('Head', Cancer.Types))
# kidney <- dplyr::filter(loulou, grepl('Kidney', Cancer.Types))
# liver <- dplyr::filter(loulou, grepl('Liver', Cancer.Types))
# lung <- dplyr::filter(loulou, grepl('Lung', Cancer.Types))
# lymph <- dplyr::filter(loulou, grepl('Lymph', Cancer.Types))
# mesothelium <- dplyr::filter(loulou, grepl('Mesothelium', Cancer.Types))
# ovary <- dplyr::filter(loulou, grepl('Ovary', Cancer.Types))
# panc <- dplyr::filter(loulou, grepl('Panc', Cancer.Types))
# prost <- dplyr::filter(loulou, grepl('Prost', Cancer.Types))
# sarcoma <- dplyr::filter(loulou, grepl('Sarcoma', Cancer.Types))
# skin <- dplyr::filter(loulou, grepl('Skin', Cancer.Types))
# stomach <- dplyr::filter(loulou, grepl('Stomach', Cancer.Types))
# thy <- dplyr::filter(loulou, grepl('Thy', Cancer.Types))
# thymoma <- dplyr::filter(loulou, grepl('Thymoma', Cancer.Types))
# transitionalCC <- dplyr::filter(loulou, grepl('Transitional-cell-carcinoma', Cancer.Types))
# UCS <- dplyr::filter(loulou, grepl('UCS', Cancer.Types))
# uterus <- dplyr::filter(loulou, grepl('Uterus', Cancer.Types))
# ALL <- dplyr::filter(loulou, grepl('ALL', Cancer.Types))
# bladder <- dplyr::filter(loulou, grepl('Bladder', Cancer.Types))
# Blood <- dplyr::filter(loulou, grepl('Blood', Cancer.Types))
# Ewings <- dplyr::filter(loulou, grepl('Ewings', Cancer.Types))
# Neuroblastoma <- dplyr::filter(loulou, grepl('Neuroblastoma', Cancer.Types))
# Oral <- dplyr::filter(loulou, grepl('Oral', Cancer.Types))
# Bone <- dplyr::filter(loulou, grepl('Bone', Cancer.Types))
# Myeloid <- dplyr::filter(loulou, grepl('Myeloid', Cancer.Types))
# SoftTissue <- dplyr::filter(loulou, grepl('SoftTissue', Cancer.Types))



#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

#Calculate the proportions for exome and genome
#exome
eraw=breast_exome
#remove zeros
eraw_violin <- eraw[apply(eraw!=0, 1, all),] #312 lines
eraw_tot <- cast(eraw_violin, id=c(Cancer.Types, Sample.Names)~variable, value="value") #176 samples in total after filtering

#genome
graw=breast_genome
#remove zeros
graw_violin <- graw[apply(graw!=0, 1, all),] #3549 lines
graw_tot <- cast(graw_violin, id=c(Cancer.Types, Sample.Names)~variable, value="value") #698 samples in total after filtering

#match names between eraw and graw
p= add_row(eraw_violin, Cancer.Types= "Breast_cancer", variable= "SBS17a", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS17b", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS8", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS9", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS18", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS21", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS34", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS37", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS40", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS41", value=0)
p= add_row(p, Cancer.Types= "Breast_cancer", variable= "SBS60", value=0)

#match names between eraw and graw
k= add_row(graw_violin, Cancer.Types= "Breast_cancer", variable= "SBS10a", value=0)
k= add_row(k, Cancer.Types= "Breast_cancer", variable= "SBS10b", value=0)
k= add_row(k, Cancer.Types= "Breast_cancer", variable= "SBS19", value=0)
k= add_row(k, Cancer.Types= "Breast_cancer", variable= "SBS29", value=0)
k= add_row(k, Cancer.Types= "Breast_cancer", variable= "SBS45", value=0)
k= add_row(k, Cancer.Types= "Breast_cancer", variable= "SBS59", value=0)


#---------------------------------------------------------------
#---------------------------------------------------------------
#plot violine plot (exome)

vp <- ggplot(p, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(trim=FALSE, size=0.3, scale="width")+ #remove the width for exome
  geom_boxplot(width=0.2, fill="white")+
  labs(title="Breast cancer (exome)",x="", y = "Log (Mutations count/Mb)")+ 
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), #setting the scale similar to the genome,
                labels = trans_format("log10", math_format(10^.x)),
                limits= c(0.0001,1000)) +
  # scale_x_discrete(expand=c(0, .6)) +
  theme(panel.grid.major.y = element_line(color="grey", size=.1), panel.background = element_blank(),
        axis.text.x = element_text(color="black", size=12, angle=90, vjust=0.5, hjust=0),
        axis.text.y = element_text(color="black", size=11), axis.line = element_line(colour = "black"),
        axis.title.y= element_text(color="black", size=12, vjust=3), axis.ticks = element_line(color="black"), legend.position = "none",
        panel.grid.major.x = element_line(color="grey", size=.1)) #or "top"
vpe=vp + geom_hline(yintercept=1, linetype="dashed", color = "red")+ scale_fill_manual(name= "variable", values = c("SBS1"="red2",
                                                                                                                    "SBS2"="green4",
                                                                                                                    "SBS3"="gold1",
                                                                                                                    "SBS4"="dodgerblue2",
                                                                                                                    "SBS5"="blueviolet",
                                                                                                                    "SBS6"="skyblue2",
                                                                                                                    "SBS7a"="palegreen2",
                                                                                                                    "SBS7b"="cadetblue1",
                                                                                                                    "SBS7c"="gainsboro",
                                                                                                                    "SBS7d"="gray70",
                                                                                                                    "SBS8"="maroon",
                                                                                                                    "SBS9"="darkturquoise",
                                                                                                                    "SBS10a"="darkorange4",
                                                                                                                    "SBS10b"="magenta",
                                                                                                                    "SBS11"="indianred1",
                                                                                                                    "SBS12"="mediumorchid1",
                                                                                                                    "SBS13"="tan1", #"deeppink2
                                                                                                                    "SBS14"="moccasin",
                                                                                                                    "SBS15"="darkolivegreen1",
                                                                                                                    "SBS16"="khaki2",
                                                                                                                    "SBS17a"="orchid1",
                                                                                                                    "SBS17b"="green1",
                                                                                                                    "SBS18"="brown",
                                                                                                                    "SBS19"="darkolivegreen4",
                                                                                                                    "SBS20"="tan4",
                                                                                                                    "SBS21"="firebrick4",
                                                                                                                    "SBS22"="darkgray",
                                                                                                                    "SBS23"="mediumvioletred",
                                                                                                                    "SBS24"="tan2",
                                                                                                                    "SBS25"="lightcyan",
                                                                                                                    "SBS26"="deeppink1",
                                                                                                                    "SBS27"="yellow4",
                                                                                                                    "SBS28"="grey",
                                                                                                                    "SBS29"="chocolate2",
                                                                                                                    "SBS30"="darkblue",
                                                                                                                    "SBS31"="yellowgreen",
                                                                                                                    "SBS32"="wheat4",
                                                                                                                    "SBS33"="seagreen",
                                                                                                                    "SBS34"="tomato3",
                                                                                                                    "SBS35"="lightcyan4",
                                                                                                                    "SBS36"="blue1",
                                                                                                                    "SBS37"="yellow3",
                                                                                                                    "SBS38"="forestgreen",
                                                                                                                    "SBS39"="lightcyan3",
                                                                                                                    "SBS40"="lightcyan2",
                                                                                                                    "SBS41"="lightsalmon",
                                                                                                                    "SBS42"="steelblue4",
                                                                                                                    "SBS43"="chocolate3",
                                                                                                                    "SBS44"="brown4",
                                                                                                                    "SBS45"="tan3",
                                                                                                                    "SBS46"="chartreuse",
                                                                                                                    "SBS47"="orange",
                                                                                                                    "SBS48"="chocolate1",
                                                                                                                    "SBS49"="seagreen1",
                                                                                                                    "SBS50"="chocolate",
                                                                                                                    "SBS51"="chocolate4",
                                                                                                                    "SBS52"="cornflowerblue",
                                                                                                                    "SBS53"="#6A3D9A",
                                                                                                                    "SBS54"="#FF7F00",
                                                                                                                    "SBS55"="#FB9A99",
                                                                                                                    "SBS56"="#CAB2D6",
                                                                                                                    "SBS57"="#FDBF6F",
                                                                                                                    "SBS58"="#DDAD4B",
                                                                                                                    "SBS59"="#7CE3D8",
                                                                                                                    "SBS60"="blue"))
vpe

ggsave("Breast_exome_violinplot_additional_colormatch_0.85_50mutation_noprop.png", width = 25, height = 10, units = "cm", dpi=600)



#---------------------------------------------------------------
#plot violine plot (genome)
vp2 <- ggplot(k, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(trim=FALSE, size=0.3, scale="width")+ #remove the width for exome
  geom_boxplot(width=0.2, fill="white")+
  labs(title="Breast cancer (genome)",x="", y = "Log (Mutations count/Mb)")+ 
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), #setting the scale similar to the genome,
                labels = trans_format("log10", math_format(10^.x)),
                limits= c(0.0001,1000)) +
  # scale_x_discrete(expand=c(0, .6)) +
  theme(panel.grid.major.y = element_line(color="grey", size=.1), panel.background = element_blank(),
        axis.text.x = element_text(color="black", size=12, angle=90, vjust=0.5, hjust=0),
        axis.text.y = element_text(color="black", size=11), axis.line = element_line(colour = "black"),
        axis.title.y= element_text(color="black", size=12, vjust=3), axis.ticks = element_line(color="black"), legend.position = "none",
        panel.grid.major.x = element_line(color="grey", size=.1)) #or "top"
vpg=vp2 + geom_hline(yintercept=1, linetype="dashed", color = "red")+ scale_fill_manual(name= "variable", values = c("SBS1"="red2",
                                                                                                                     "SBS2"="green4",
                                                                                                                     "SBS3"="gold1",
                                                                                                                     "SBS4"="dodgerblue2",
                                                                                                                     "SBS5"="blueviolet",
                                                                                                                     "SBS6"="skyblue2",
                                                                                                                     "SBS7a"="palegreen2",
                                                                                                                     "SBS7b"="cadetblue1",
                                                                                                                     "SBS7c"="gainsboro",
                                                                                                                     "SBS7d"="gray70",
                                                                                                                     "SBS8"="maroon",
                                                                                                                     "SBS9"="darkturquoise",
                                                                                                                     "SBS10a"="darkorange4",
                                                                                                                     "SBS10b"="magenta",
                                                                                                                     "SBS11"="indianred1",
                                                                                                                     "SBS12"="mediumorchid1",
                                                                                                                     "SBS13"="tan1", #"deeppink2
                                                                                                                     "SBS14"="moccasin",
                                                                                                                     "SBS15"="darkolivegreen1",
                                                                                                                     "SBS16"="khaki2",
                                                                                                                     "SBS17a"="orchid1",
                                                                                                                     "SBS17b"="green1",
                                                                                                                     "SBS18"="brown",
                                                                                                                     "SBS19"="darkolivegreen4",
                                                                                                                     "SBS20"="tan4",
                                                                                                                     "SBS21"="firebrick4",
                                                                                                                     "SBS22"="darkgray",
                                                                                                                     "SBS23"="mediumvioletred",
                                                                                                                     "SBS24"="tan2",
                                                                                                                     "SBS25"="lightcyan",
                                                                                                                     "SBS26"="deeppink1",
                                                                                                                     "SBS27"="yellow4",
                                                                                                                     "SBS28"="grey",
                                                                                                                     "SBS29"="chocolate2",
                                                                                                                     "SBS30"="darkblue",
                                                                                                                     "SBS31"="yellowgreen",
                                                                                                                     "SBS32"="wheat4",
                                                                                                                     "SBS33"="seagreen",
                                                                                                                     "SBS34"="tomato3",
                                                                                                                     "SBS35"="lightcyan4",
                                                                                                                     "SBS36"="blue1",
                                                                                                                     "SBS37"="yellow3",
                                                                                                                     "SBS38"="forestgreen",
                                                                                                                     "SBS39"="lightcyan3",
                                                                                                                     "SBS40"="lightcyan2",
                                                                                                                     "SBS41"="lightsalmon",
                                                                                                                     "SBS42"="steelblue4",
                                                                                                                     "SBS43"="chocolate3",
                                                                                                                     "SBS44"="brown4",
                                                                                                                     "SBS45"="tan3",
                                                                                                                     "SBS46"="chartreuse",
                                                                                                                     "SBS47"="orange",
                                                                                                                     "SBS48"="chocolate1",
                                                                                                                     "SBS49"="seagreen1",
                                                                                                                     "SBS50"="chocolate",
                                                                                                                     "SBS51"="chocolate4",
                                                                                                                     "SBS52"="cornflowerblue",
                                                                                                                     "SBS53"="#6A3D9A",
                                                                                                                     "SBS54"="#FF7F00",
                                                                                                                     "SBS55"="#FB9A99",
                                                                                                                     "SBS56"="#CAB2D6",
                                                                                                                     "SBS57"="#FDBF6F",
                                                                                                                     "SBS58"="#DDAD4B",
                                                                                                                     "SBS59"="#7CE3D8",
                                                                                                                     "SBS60"="blue"))
vpg

ggsave("Breast_genome_violinplot_additional_colormatch_0.85_50mutation_noprop.png", width = 25, height = 10, units = "cm", dpi=600)

#---------------------------------------------------------------
#merge both plots
#---------------------------------------------------------------
plot_grid(vpg, vpe, ncol=1, nrow = 2)

ggsave("Breast_genome_exome_violinplot_additional_colormatch_0.85_50mutation_noprop.png", width = 30, height = 20, units = "cm", dpi=600)
