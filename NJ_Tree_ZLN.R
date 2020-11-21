setwd("~/Dropbox/CastoeLabFolder/projects/schist_project/Manuscripts/_MS2_Shist_PopGen_Jonathan/Schisto_NJ_tree/")
library(RColorBrewer)
library(vcfR)
library(ape)
library(adegenet)
schisto.FULL.vcf <- read.vcfR("schisto_FULL.vcf")
schisto.gl <- vcfR2genlight(schisto.FULL.vcf)
schisto.gl
Schisto_distance<- dist(schisto.gl)

tree <- nj(Schisto_distance)
tree
plot(tree)
key <- read.table("sample_key_2018.txt", header=T)
tree$tip.names <- tree$tip.label
cl <- colors(distinct = TRUE)

# get a set of colors used to distinguish villages
set.seed(15885) # to set random generator seed, this ensures that the set of random colors chosen can be reproduced
#village_colors <- sample(cl, max(key$vill_col))
vill_colors = c("110107" = "#4E79A7", "110501" = "#86BCB6", "121105" = "#499894", "130408" = "#F1CE63", "130903" = "#59A14F", "210101" = "#9D7660", "210107" = "#E15759", "210302" = "#79706E", "210304" = "#D37295", "220504" = "#B07AA1", "220605" = "#F28E2B", "230703" = "#8CD17D","Other"="Grey50")
year_colors= c("2007"="blue", "2008" ="darkgray", "2010"="darkgreen", "2016F" = "goldenrod", "2016S"= "purple")
host_colors= c("1602"="#AA4B49","5101"="#C08F18","601"="#7769B7","3302"="#CD8FF2","5001"="#A14992","701"="#933E4D","2003"="#D74852","2002"="#A24E2F","702"="#9F8E3D","5501"="#4F4579","101"="#A32DDF","501"="#9F4C84","901"="#BF2A3F","1001"="#89020A","1701"="#E5762D","2801"="#C9E91E","7001"="#452F86","8301"="#7F03B0","8601"="#86436E","2201"="#B54C59","3001"="#944346","3601"="#7B6847","3602"="#54CF3E","301"="#8767DA","9201"="#9A02D5","6301"="#892E5A","2102"="#D44254","2501"="#A17819","2901"="#3FA9DD","3201"="#65439B","3401"="#7C3D8B","3502"="#96164B")

 

village_tip_colors <- vill_colors[key$vill_col[match(tree$tip.names,key$idmiracid)]]
year_tip_colors <- year_colors[key$samplewave[match(tree$tip.names,key$idmiracid)]]     
host_tip_colors <- host_colors[key$host_col[match(tree$tip.names,key$idmiracid)]]

village_edge_colors <- c(village_tip_colors, rep("gray", nrow(tree$edge)-length(village_tip_colors)+1))
year_edge_colors <- c(year_tip_colors, rep("gray", nrow(tree$edge)-length(year_tip_colors)+1))
host_edge_colors <- c(host_tip_colors, rep("gray", nrow(tree$edge)-length(host_tip_colors)+1))

# By Village
par(mar=c(0,0,0,0))
plot.phylo(tree, type="unrooted", cex = 0.3, 
           show.tip.label = F, edge.width = 2,
           edge.col=village_edge_colors[tree$edge[,2]])
tiplabels( pch=21, col=village_tip_colors, bg=village_tip_colors, cex=0.9)
legend("right", inset=.02, title="Village",
       c("110107","110501","121105","130408","130903","210101","210107","210302","210304","220504","220605","230703"), fill=c("#4E79A7","#86BCB6","#499894","#F1CE63","#59A14F","#9D7660","#E15759","#79706E","#D37295","#B07AA1","#F28E2B","#8CD17D"), horiz=FALSE, cex=1.25)

# By Year
plot.phylo(tree, type="unrooted", cex = 0.3, 
           show.tip.label = F, edge.width = 2,
           edge.col=year_edge_colors[tree$edge[,2]])
tiplabels( pch=21, col=year_tip_colors, bg=year_tip_colors, cex=0.9)
legend("right", inset=.02, title="Year",
       c("2007","2008","2010", "2016F", "2016S"), fill=c("blue", "darkgray", "darkgreen","goldenrod", "purple"), horiz=FALSE, cex=1.25)

# By Host
par(mar=c(0,0,0,0))
plot.phylo(tree, type="unrooted", cex = 0.3, 
           show.tip.label = F, edge.width = 2,
           edge.col=host_edge_colors[tree$edge[,2]])
tiplabels( pch=21, col=host_tip_colors, bg=host_tip_colors, cex=0.9)
legend("right", inset=.02, title="Host",
       c("1602","5101", "601","3302","5001", "701","2003","2002", "702","5501", "101", "501", "901","1001","1701","2801","7001","8301","8601","2201","3001","3601","3602", "301","9201","6301","2102","2501","2901","3201","3401","3502"), fill=c("#AA4B49","#C08F18","#7769B7","#CD8FF2","#A14992","#933E4D","#D74852","#A24E2F","#9F8E3D","#4F4579","#A32DDF","#9F4C84","#BF2A3F","#89020A","#E5762D","#C9E91E","#452F86","#7F03B0","#86436E","#B54C59","#944346","#7B6847","#54CF3E","#8767DA","#9A02D5","#892E5A","#D44254","#A17819","#3FA9DD","#65439B","#7C3D8B","#96164B"), horiz=FALSE, cex=1.25)

# Branch = Village Tips = Year
plot.phylo(tree, type="unrooted", cex = 0.3, 
           show.tip.label = F, edge.width = 2,
           edge.col=village_edge_colors[tree$edge[,2]])
tiplabels( pch=21, col=year_tip_colors, bg=year_tip_colors, cex=0.9)
