d = ("/media/archeogen/Elements/PROJECTS/Zalavár_finalcut/IBD/")
setwd(d)
tab1 = read.table("filtered_IBD.tsv", stringsAsFactors = F, header = T, sep = "\t")
tab2 = read.table("names_filtered_IBD_backup.tsv", stringsAsFactors = F, header = T, sep = "\t")
tab2 = tab2[!(tab2$Period3 == 250),]
cents = sort((unique(tab2$Period3_4)))
fin = as.data.frame(matrix(NA, length(cents), 6))
rownames(fin) = cents
colnames(fin) = c("AC_Eu_TD", "AC_Eu_CD", "AC_Eu", "AC_All_TD", "AC_All_CD", "AC_All")

#AC_Eu (repeat for all types, AC - cumulative kB, Eu - CB-EUR, All - all, CD - GHP, TD - Transdanubia)
# (sum(connections(ix))/(ni x nx)
subtab = NULL
for(i in 1:nrow(fin)){
  subtab = rbind(subtab, tab2[tab2$Period3_4 == rownames(fin)[i],]) #get the samples of the correspondent period
  subtab = subtab[!(subtab$f4_clust == "NO" | subtab$f4_clust == "YES"),] #exclude non european group
  subid = tab1[tab1$iid1 %in% subtab$Edges,] #get list
  subid = subid[subid$iid2 %in% subtab$Edges,]
  subid = subid[!(subid$sum_IBD.8 < 10),] #exclude no hits
  subid = subid[!(subid$sum_IBD.8 >= 1400),] #exclude close relatives
  groups = unique(subtab$Group2)
  counter = 0
  for(j in 1:length(groups)){
    subsubtab = subtab[subtab$Group2 == groups[j],]
    subsubid1 = subid[subid$iid1 %in% subsubtab$Edges,] #first col
    subsubid2 = subid[subid$iid2 %in% subsubtab$Edges,] #second col
    subsubid3 = subsubid1[subsubid1$iid2 %in% subsubtab$Edges,] #overlap, i.e. ingroup connections
    counter = counter + ((nrow(subsubid1) + nrow(subsubid2) - 2*nrow(subsubid3))/(nrow(subsubtab)*(nrow(subtab)-nrow(subsubtab)))) #(col1 + col2 - overlapx2)/(nrow1*(nrow2-nrow1))
  }
  fin$AC_Eu[i] = counter
}

#AC_Eu
# (sum(connections(ix)))/2n
subtab = NULL
for(i in 1:nrow(fin)){
  subtab = rbind(subtab, tab2[tab2$Period3 == rownames(fin)[i],]) #get the samples of the correspondent period
  subtab = subtab[!(subtab$f4_clust == "NO" | subtab$f4_clust == "YES"),] #exclude non european group
  subid = tab1[tab1$iid1 %in% subtab$Edges,] #get list
  subid = subid[subid$iid2 %in% subtab$Edges,]
  subid = subid[!(subid$sum_IBD.8 < 10),] #exclude no hits
  subid = subid[!(subid$sum_IBD.8 >= 1400),] #exclude close relatives
  groups = unique(subtab$Group2)
  counter = 0
  for(j in 1:length(groups)){
    subsubtab = subtab[subtab$Group2 == groups[j],]
    subsubid1 = subid[subid$iid1 %in% subsubtab$Edges,] #first col
    subsubid2 = subid[subid$iid2 %in% subsubtab$Edges,] #second col
    subsubid3 = subsubid1[subsubid1$iid2 %in% subsubtab$Edges,] #overlap, i.e. ingroup connections
    counter = counter + (nrow(subsubid1) + nrow(subsubid2) - 2*nrow(subsubid3)) #col1 + col2 - overlapx2
  }
  fin$AC_Eu[i] = (counter/2)/(nrow(subtab))
}