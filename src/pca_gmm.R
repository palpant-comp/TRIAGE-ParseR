library(pheatmap)

### Jaccard index heatmap
jaccard_table=function(input_table, top_no=100){
  results = data.frame(matrix(ncol=ncol(input_table), nrow=ncol(input_table)))
  rownames(results) = colnames(input_table)
  colnames(results) = colnames(input_table)
  for (col1 in colnames(input_table)){
    r1 = rownames(input_table)[order(input_table[,col1], decreasing=T)[1:top_no]]
    for (col2 in colnames(input_table)){
      r2 = rownames(input_table)[order(input_table[,col2], decreasing=T)[1:top_no]]
      value = length(intersect(r1,r2))/length(unique(c(r1,r2)))  # Jaccard index
      results[col1,col2] = value
    }
  }
  return (results)
}

temp = read.table('/Users/woojunshim/Research/extended_TRIAGE/data/sophie/211206_rrcpint0Xbc_meanpeakdisc.txt', stringsAsFactors = F, row.names=NULL)
aa = table(temp$row.names)
temp = temp[which(temp$row.names %in% names(aa)[which(aa==1)]),]
rownames(temp) = temp$row.names
temp = temp[,2:ncol(temp)]
corr_table = jaccard_table(temp)
pdf('/Users/woojunshim/Research/extended_TRIAGE/data/sophie/211206_rrcpint0Xbc_meanpeakdisc.pdf', width=10, height=10)
pheatmap(corr_table)
dev.off()

### PCA GMM GO SIGNIFICANT HEATMAP
library(viridis)
groups = c('37','11','20','21','25','26','38','42','55','62','64')
for (group in groups){
  table_ = read.table(paste0('/Users/woojunshim/Research/sc_clustering_paper/data/go/',group,'_go.txt'), stringsAsFactors = F)
  table_ = -log10(table_)
  rownames(table_) = gsub('_',' ',rownames(table_))
  idx = idx_for_top_values(matrix=table_, no_terms=10, col_order=colnames(table_))
  pdf(paste0('/Users/woojunshim/Research/sc_clustering_paper/data/gmm_figures/',group,'_go_heatmap.pdf'), width=7, height=6)
  pheatmap(table_[idx,], cluster_cols = F, cluster_rows = F, color=inferno(10), border_color = NA)
  dev.off()
}
temp = read.table('/Users/woojunshim/Research/sc_clustering_paper/data/go/64_go.txt', stringsAsFactors = F)
temp = -log10(temp)
idx = order(temp$cluster3, decreasing=T)
table_ = data.frame(temp[idx,])
rownames(table_) = rownames(temp)[idx]
colnames(table_) = colnames(temp)
pdf(paste0('/Users/woojunshim/Research/sc_clustering_paper/data/gmm_figures/64_go_heatmap.pdf'), width=5, height=6)
pheatmap(table_, cluster_cols = F, cluster_rows = F, color=inferno(10))
dev.off()