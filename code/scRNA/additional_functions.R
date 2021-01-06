### Additional functions for reproduction of analyses presented in Wilk, et al. bioRxiv (2020) ###


new_dotplot <- function(object = NULL, features = NULL, group.by = NULL, genes.on.x = TRUE, 
                        size.breaks.values = NULL, color.breaks.values = c(-3, -2, -1, 0, 1, 2, 3), shape.scale = 12, 
                        dend_x_var = "Average expression", dend_y_var = "Average expression",
                        cols.use = c("lightgrey", "blue"), scale.by = "radius", col.min = -2.5, col.max = 2.5,
                        dot.min = 0) {
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
  data.features <- FetchData(object = object, vars = features)
  object[[group.by, drop = TRUE]]
  data.features$id <- object[[group.by, drop = TRUE]]
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if(genes.on.x){
    data.final <- data.frame(features.plot = data.plot$features.plot, id = data.plot$id)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  else {
    data.final <- data.frame(id = data.plot$id, features.plot = data.plot$features.plot)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  colnames(data.final)[3:4] <- c("Percent expressed", "Average expression")
  dot_plot(data.final, size_var = "Percent expressed", "Average expression", 
           dend_y_var = dend_y_var, dend_x_var = dend_x_var,
           dist = "euclidean", hclust_method = "ward.D2", x.lab.pos = "bottom",
           display_max_sizes = FALSE, size.breaks.values = size.breaks.values,
           shape.scale = shape.scale, color.breaks.values = color.breaks.values, cols.use = cols.use)
}


filter_quality <- function(.data) {
  final <- .data[-c(grep("^RPS", .data$gene), 
                    grep("^RPL", .data$gene), 
                    grep("^MT-", .data$gene),
                    grep("^MTR", .data$gene),
                    grep("MALAT1", .data$gene),
                    grep("^RNA18S5", .data$gene),
                    grep("^RNA28S5", .data$gene)),] 
  return(final)
}

plotPosCells <- function(seu, emat, features, grep = F) {
  if (grep) {
    pos <- emat[grep(paste(features,collapse="|"), rownames(emat)),]
    pos.cells <- colnames(pos[,colSums(pos) != 0])
  }
  else {
    if (length(features) >1) {
      pos <- emat[features,]
      pos.cells <- colnames(pos[,colSums(pos) != 0])
    }
    else {
      pos <- emat[features,] != 0
      pos.cells <- names(pos)[pos==T]
    }
    
  }
  
  DimPlot(seu, cells.highlight = pos.cells) + NoLegend()  
}


#base (non-numeric) Abundance function
covid.seuAbundances <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"),
                                 meta.include = NULL, group_by = NULL, shape_by = NULL, each.pt = "orig.ident",
                                 custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                 pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                 ncol = 4, label = "p.signif", select.idents = NULL, 
                                 label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               each.pt))  
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               each.pt))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=max(uniques[c(each.pt,"seurat_clusters")])])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
    df$cell.type <- factor(df$cell.type, levels = select.idents)
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = each.pt, fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    
    
  }
  
} 


#numeric abundance function (used when the x-axis is numeric, correlation coeffs are calculated)
covid.seuAbundances.num <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"),
                                     meta.include = NULL, group_by = NULL, shape_by = NULL,
                                     custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                     pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                     ncol = 4, label = "p.signif", select.idents = NULL, 
                                     label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               "orig.ident"))  
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               "orig.ident"))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=uniques["seurat_clusters"]])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    if(is.numeric(df[,group_by])){
      
      
      p + facet_wrap("cell.type", scales = "free_y", ncol = ncol) + 
        guides(fill = FALSE) +
        geom_point(size = pt.size, position = position_jitter(width = 0.25), 
                   aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
        scale_shape_manual(values = shapes) + 
        theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
        scale_color_manual(values = custom_fill_colors) + 
        scale_fill_manual(values = custom_fill_colors) + 
        ggpubr::stat_cor(label.x = label.x)
    }  
    
  }
  
} 

#used when dropping in non-proportion information 
covid.seuAbundances.nonprop <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"), 
                                         non.prop.data = NULL, meta.include = NULL, group_by = NULL, shape_by = NULL, 
                                         each.pt = "orig.ident", custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                         pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                         ncol = 4, label = "p.signif", select.idents = NULL, 
                                         label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               each.pt))  
  }
  if (!is.null(non.prop.data)) {
    df <- aggregate(seu@meta.data[,non.prop.data], 
                    by = list(seu@meta.data[,by], seu@meta.data[,each.pt]), 
                    FUN = mean)
    colnames(df) <- c(by,each.pt,"freq")
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               each.pt))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=max(uniques[c(each.pt,"seurat_clusters")])])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = each.pt, fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    
    
  }
  
} 

#used when dropping in non-proportion information 
covid.seuAbundances.nonprop <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"), 
                                         non.prop.data = NULL, meta.include = NULL, group_by = NULL, shape_by = NULL, 
                                         each.pt = "orig.ident", custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                         pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                         ncol = 4, label = "p.signif", select.idents = NULL, 
                                         label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  
  if (!is.null(non.prop.data)) {
    if(!non.prop.data %in% colnames(seu@meta.data)) {
      #exprs.query <- log(GetAssayData(seu, assay = "RNA", slot = "data")[non.prop.data,]+1)
      exprs.query <- GetAssayData(seu, assay = "SCT", slot = "data")[non.prop.data,]
      seu <- AddMetaData(seu, metadata = exprs.query, col.name = non.prop.data)
    }
    
    df <- aggregate(seu@meta.data[,non.prop.data], 
                    by = list(seu@meta.data[,by], seu@meta.data[,each.pt]), 
                    FUN = median)
    colnames(df) <- c(by,each.pt,"freq")
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               each.pt))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=max(uniques[c(each.pt,"seurat_clusters")])])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = each.pt, fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    
    
  }
  
} 

#used when dropping in non-proportion information and using a numeric x-axis 
covid.seuAbundances.nonprop.num <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"), 
                                             non.prop.data = NULL, meta.include = NULL, group_by = NULL, shape_by = NULL, 
                                             each.pt = "orig.ident", custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                             pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                             ncol = 4, label = "p.signif", select.idents = NULL, 
                                             label.x = NA, pt.size = NA, label_size = NULL) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               each.pt))  
  }
  if (!is.null(non.prop.data)) {
    df <- aggregate(seu@meta.data[,non.prop.data], 
                    by = list(seu@meta.data[,by], seu@meta.data[,each.pt]), 
                    FUN = mean)
    colnames(df) <- c(by,each.pt,"freq")
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               each.pt))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=max(uniques[c(each.pt,"seurat_clusters")])])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = each.pt, fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    if(is.numeric(df[,group_by])){
      
      
      p + facet_wrap("cell.type", scales = "free_y", ncol = ncol) + 
        guides(fill = FALSE) +
        geom_point(size = pt.size, position = position_jitter(width = 0.25), 
                   aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
        scale_shape_manual(values = shapes) + 
        theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
        scale_color_manual(values = custom_fill_colors) + 
        scale_fill_manual(values = custom_fill_colors) + 
        ggpubr::stat_cor(label.x = label.x, aes_string(color=color_by)) + 
        geom_smooth(method = lm, aes_string(fill=color_by, color = color_by), size = label_size)
    }  
    
  }
  
} 


writeAnnData <- function(emat = NULL, nmat = NULL, seu = NULL, idents = NULL, obs.vars = c("orig.ident", "seurat_clusters", "cell.type"), name = "seu", dir = NULL) {
  subset <- Seurat:::subset.Seurat(seu, idents = idents)
  seu.emat <- emat[,which(colnames(emat) %in% colnames(subset@assays$RNA)), drop=FALSE] 
  seu.nmat <- nmat[,which(colnames(nmat) %in% colnames(subset@assays$RNA)), drop=FALSE]
  seu.emat <- seu.emat[which(rownames(seu.emat) %in% rownames(seu.nmat)),, drop = FALSE]
  seu.nmat <- seu.nmat[which(rownames(seu.nmat) %in% rownames(seu.emat)),, drop = FALSE]
  seu.nmat <- seu.nmat[match(rownames(seu.emat), rownames(seu.nmat)),]
  if (!identical(colnames(seu.emat), colnames(seu.nmat), FALSE, FALSE, FALSE, FALSE)) {
    stop("Process failed: colnames not identical")
  }
  if(!identical(rownames(seu.emat), rownames(seu.nmat), FALSE, FALSE, FALSE, FALSE)) {
    stop("Process failed: rownames not identical")
  }
  
  #seu.obs <- as.data.frame(cbind(subset@meta.data$orig.ident, subset@meta.data$seurat_clusters, subset@meta.data$cell.type))
  seu.obs <- subset@meta.data[,colnames(subset@meta.data) %in% obs.vars]
  colnames(seu.obs)[grep("seurat_clusters", colnames(seu.obs))] <- "clusters"
  #colnames(seu.obs) <- obs.vars
  rownames(seu.obs) <- rownames(subset@meta.data)
  seu.obs <- seu.obs[match(colnames(seu.emat), rownames(seu.obs)),]
  seu.obsm <- as.data.frame(Embeddings(subset, reduction = "umap"))
  seu.obsm <- seu.obsm[match(colnames(seu.emat), rownames(seu.obsm)),]
  seu.obs_names <- as.data.frame(colnames(seu.emat))
  seu.var_names <- as.data.frame(rownames(seu.emat))
  
  writeMM(t(seu.emat), file = paste0(dir,name,".emat.mtx"))
  writeMM(t(seu.nmat), file = paste0(dir,name,".nmat.mtx"))
  write.csv(seu.obs, file =  paste0(dir,name,".obs.csv"), row.names = FALSE)
  write.csv(seu.obsm, file =  paste0(dir,name,".obsm.csv"), row.names = FALSE)
  write.table(seu.obs_names, file = paste0(dir,name,".obs_names.txt"), row.names = FALSE, col.names = FALSE)
  write.table(seu.var_names, file = paste0(dir,name,".var_names.txt"), row.names = FALSE, col.names = FALSE)
  message("Finished writing data")
}


build.deg.matrix <- function(compartment = NULL, p.cutoff = 0.05, 
                             path = "/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/COVID_sc/scRNA/analyses/DEG_results_final/"){
  data.load <- paste0(path, list.files(path = path, pattern = compartment))
  markers.list <- lapply(data.load, read_csv)
  markers.list <- ldply(markers.list, data.frame)
  
  markers.list.mat <- markers.list[markers.list$p_val_adj<p.cutoff,]
  markers.list.mat <- markers.list.mat[,c("gene", "avg_logFC", "Donor")]
  
  markers.list.mat <- reshape2::dcast(markers.list.mat, formula = gene~Donor, value.var = "avg_logFC")
  markers.list.mat[is.na(markers.list.mat)] = 0
  markers.list.mat <- markers.list.mat %>% column_to_rownames(var = "gene")
  colnames(markers.list.mat) <- gsub("covid_", "", colnames(markers.list.mat))
  return(markers.list.mat)
}


covid.markers.heatmap <- function(markers.matrix = NULL, color = c("blue", "white", "red"), 
                                  paletteLength = 100, fontsize_row = 10, color.rows = T,
                                  title = "log(FC)", add.flag = T, de.cutoff = 4, top.n = NULL, 
                                  repel.degree = 0, legend = T, annotation_legend = T,
                                  cellwidth = NA, cellheight = NA,
                                  save = F, width = 7, height = 11, file = "~/Downloads/p.pdf") {
  
  paletteLength = 100
  if(sum(markers.matrix>=0)==dim(markers.matrix)[1]*dim(markers.matrix)[2]) {
    color = color[2:3]
  }
  myColor = colorRampPalette(color)(paletteLength)
  myBreaks <- unique(c(seq(min(markers.matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(max(markers.matrix)/paletteLength, max(markers.matrix),
                           length.out=floor(paletteLength/2))))
  annotation_colors = list(
    Ventilated = c(ARDS="red3", NonVent=RColorBrewer::brewer.pal(9, "Oranges")[4]))
  p <- pheatmap(markers.matrix, color = myColor, breaks = myBreaks, 
                heatmap_legend_param = list(title = title), angle_col = "90", 
                annotation_col = row_annotation, annotation_colors = annotation_colors, legend = legend,
                annotation_legend = annotation_legend, fontsize_row = fontsize_row,
                cellwidth = cellwidth, cellheight = cellheight)
  if(color.rows){
    markers.matrix <- markers.matrix[match(p$gtable$grobs[[5]]$label,rownames(markers.matrix)),]
    p$gtable$grobs[[5]]$gp=gpar(col=ifelse((rowSums(markers.matrix))>0, "red", "blue"), fontsize = fontsize_row)
  }
  
  if(add.flag){
    if(save){
      pdf(file = file, width = width, height = height)
    }
    if(!is.null(top.n)) {
      kept.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    add.flag(p, kept.labels = kept.labels,
             repel.degree = repel.degree)
    if(save){
      dev.off()
    }
  }
  else{
    p
  }
}



processNewSeurat <- function(parent.object, idents = NULL, cells = NULL) {
  if (!is.null(idents)) {
    seu <- subset(parent.object, idents = idents)
  }
  if (!is.null(cells)) {
    seu <- subset(parent.object, cells = cells)
  }
  message("Running PCA")
  seu <- RunPCA(seu, verbose = FALSE)
  message("Running UMAP")
  seu <- RunUMAP(seu, dims = 1:50, verbose = FALSE)
  message("Clustering")
  seu <- FindNeighbors(seu, dims = 1:50, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
  return(seu)
}



myFeaturePlot <- function(object = NULL, features = features, ncol = NULL, save = T, save.as = "png", height = 7, width = 7, cols = c("lightgrey", "blue"), reduction = "umap", facet.size = 20, feature.face = "bold.italic") {
  plots <- lapply(features, function(x) {FeaturePlot(object,features = x, cols = cols, reduction = reduction) +  labs(x = "UMAP1", y = "UMAP2") + 
      theme(aspect.ratio = 1, 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title = element_blank(), 
            axis.line = element_blank(),
            plot.title = element_text(face = feature.face, size = facet.size))})
  CombinePlots(plots = plots, ncol = ncol)
  if(save) {
    if(save.as=="pdf"){
      ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
    }
    if(save.as=="png"){
      ggsave("p.png", path = "~/Downloads/", height = height, width = width)
    }
  }
}

myHighlightCells <- function(object = NULL, idents = NULL, group_by = "seurat_clusters", ncol = NULL, save = T, height = 7, width = 7, col = "black") {
  cells.highlight=list()
  for(i in 1:length(idents)){
    cells = colnames(object)[grep(idents[i], object@meta.data[,group_by])]
    cells.highlight[[i]]=cells
  }
  
  plots <- lapply(cells.highlight, function(x) {
    DimPlot(object,cells.highlight = x, cols.highlight = col) + 
      NoLegend() + labs(x = "UMAP1", y = "UMAP2") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  })
  CombinePlots(plots = plots, ncol = ncol)
  if(save) {
    ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
  }
}


constructConsensus <- function(markers.matrix, donors.use = 1:8, save = T, height = 5, width = 7) {
  markers.matrix.m <- markers.matrix %>% rownames_to_column(var = "gene")
  markers.matrix.m <- reshape2::melt(markers.matrix.m)
  order.markers <- rownames(markers.matrix[order(rowSums(-markers.matrix)),])
  markers.sum <- as.data.frame(ifelse(markers.matrix != 0, 1, 0))
  markers.sum$total <- rowSums(markers.sum)
  keep <- rownames(markers.sum[markers.sum$total>3,])
  markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% keep,]
  ggplot(markers.matrix.m, aes(x = factor(gene, level = order.markers), 
                               y = value, fill = variable)) +
    geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() + 
    labs(x = "", y = "Cumulative average log(fold-change)", fill = "Sample") +
    scale_fill_manual(values = custom_fill_colors[donors.use]) + 
    ggpubr::rotate_x_text()
  if(save){
    ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
  }
}



new.markers.heatmap <- function(markers.matrix = NULL, fontsize_row = 10, color.rows = T,
                                title = "log(FC)", add.flag = T, de.cutoff = 4, percentile = 0.999,
                                top.n = NULL, color = c("blue", "white", "red"), 
                                repel.degree = 0, legend = T, annotation_legend = T,
                                cellwidth = NA, cellheight = NA, markers.label = NULL,
                                save = T, width = 7, height = 11, cluster_columns = T,
                                file = "~/Downloads/p.pdf", feature.face = "italic") {
  row_annotation <- unique(row_annotation[row_annotation$Donor.full %in% colnames(markers.matrix),])
  row_annotation <- row_annotation[match(colnames(markers.matrix),row_annotation$Donor.full),]
  ta = columnAnnotation("Age" = row_annotation$Age, 
                        "Peak Severity" = row_annotation$Severity.final.max,
                        "Current Severity" = row_annotation$Severity.final.current,
                        "Gender" = row_annotation$Sex,
                        "Acuity" = row_annotation$acuity,
                        "DPT" = row_annotation$days_post_pos_test,
                        "Days To Death" = row_annotation$time_to_death,
                        col = list("Peak Severity" = c("0" = severity.cols[1],
                                                       "1-3" = severity.cols[2],
                                                       "4-5" = severity.cols[3],
                                                       "6-7" = severity.cols[4],
                                                       "8" = severity.cols[5]),
                                   "Current Severity" = c("0" = severity.cols[1],
                                                          "1-3" = severity.cols[2],
                                                          "4-5" = severity.cols[3],
                                                          "6-7" = severity.cols[4]),
                                   "Age" = colorRamp2(c(20,100),c("white","purple")),
                                   "Gender" = c("M" = "brown", "F" = "blue"),
                                   "Acuity" = c("Acute" = "orange", "Convalescent" = "purple", "Healthy" = "steelblue"),
                                   "DPT" = colorRamp2(c(0, 80),c("white","yellowgreen")),
                                   "Days To Death" = colorRamp2(c(0, 70),c("red","white"))
                        ))
  if(add.flag){
    if(!is.null(top.n)) {
      top.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
      kept.labels <- rownames(markers.matrix)[rownames(markers.matrix) %in% top.labels]
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    if(!is.null(markers.label)) {
      kept.labels <- markers.label[order(markers.label)]
    }
    t = 1:nrow(markers.matrix)
    
    ha = rowAnnotation(foo = anno_mark(at = t[rownames(markers.matrix) %in% kept.labels], 
                                       labels = kept.labels,
                                       labels_gp = gpar(col =
                                                          ifelse((rowSums(markers.matrix[rownames(markers.matrix) %in% kept.labels,]))>0,
                                                                 "red", "blue"),
                                                        fontface = feature.face)))
  }
  
  markers.matrix <- as.matrix(markers.matrix)
  col.lim = quantile(as.matrix(abs(markers.matrix)), probs = percentile, na.rm=T)
  p <- Heatmap(markers.matrix, name = title,
               right_annotation = ha, 
               top_annotation = ta, 
               row_names_gp = gpar(fontsize = 0),
               cluster_columns = cluster_columns,
               col = colorRamp2(breaks = c(-col.lim,0,col.lim), colors = color),
               column_names_gp = gpar(fontsize=0))
  if(save){
    pdf(file = file, width = width, height = height)
    p
    #dev.off()
  }
  else {
    p
  }
}


build.cp.matrix <- function(compartment = NULL, by = "z", p.cutoff = 0.05, z.cutoff = T, top.n = 100, path = "/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/COVID_sc/scRNA/analyses/DEG_results/IPA/results/") {
  files <- list.files(path)
  files_import <- grep(paste0(compartment,".cp.txt"),files,value = T)
  donor_names <- gsub(paste0("_",compartment,".cp.txt"),"",files_import)
  datalist = lapply(files_import, function(x)read_delim(paste0(path,x), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1))
  datalist <- mapply(cbind, datalist, "Donor"=donor_names,SIMPLIFY = F)
  datalist.df <- ldply(datalist,data.frame)
  datalist.df$X6 <- NULL
  colnames(datalist.df) <- mapvalues(colnames(datalist.df), from = c("Ingenuity.Canonical.Pathways",
                                                                     "X.log.p.value.",
                                                                     "Ratio",
                                                                     "z.score",
                                                                     "Molecules",
                                                                     "Donor"),
                                     to = c("Pathway",
                                            "p",
                                            "Ratio",
                                            "z",
                                            "Molecules",
                                            "Donor"))
  datalist.final <- datalist.df
  if(!is.null(p.cutoff)){
    p.converted = -log(0.05,10)
    datalist.final <- datalist.final[datalist.final$p>p.converted,]
  }
  if(z.cutoff){
    datalist.final <- datalist.final[!is.na(datalist.final$z),]
    datalist.final <- datalist.final[!datalist.final$z==0,]
  }
  datalist.final[is.na(datalist.final$z),"z"] <- 0
  datalist.mat <- reshape2::dcast(data = datalist.final, formula = Pathway~Donor,fun.aggregate = sum, value.var = by) %>% column_to_rownames(var = "Pathway")
  if(by=="p" && z.cutoff){
    datalist.z <- as.matrix(reshape2::dcast(data = datalist.final, formula = Pathway~Donor,fun.aggregate = sum, value.var = "z") %>% column_to_rownames(var = "Pathway"))
    datalist.mat <- as.matrix(datalist.mat)
    datalist.mat[datalist.z<0] <- -datalist.mat
    datalist.mat <- as.data.frame(datalist.mat)
  }
  
  if(!is.null(top.n)) {
    top <- names((abs(rowSums(datalist.mat)) %>% sort(decreasing = T))[1:top.n])
    datalist.mat <- datalist.mat[rownames(datalist.mat) %in% top,]
  }
  
  return(datalist.mat)
}

build.ur.matrix <- function(compartment = NULL, by = "z", p.cutoff = 0.05, z.cutoff = T, top.n = 100, path = "/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/COVID_sc/scRNA/analyses/DEG_results/IPA/results/") {
  files <- list.files(path)
  files_import <- grep(paste0(compartment,".ur.txt"),files,value = T)
  donor_names <- gsub(paste0("_",compartment,".ur.txt"),"",files_import)
  datalist = lapply(files_import, function(x)read_delim(paste0(path,x), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1))
  datalist <- mapply(cbind, datalist, "Donor"=donor_names,SIMPLIFY = F)
  datalist.df <- ldply(datalist,data.frame)
  colnames(datalist.df) <- mapvalues(colnames(datalist.df), from = c("Upstream.Regulator",
                                                                     "Expr.Log.Ratio",
                                                                     "Molecule.Type",
                                                                     "Predicted.Activation.State",
                                                                     "Activation.z.score",
                                                                     "p.value.of.overlap",
                                                                     "Target.Molecules.in.Dataset",
                                                                     "Mechanistic.Network",
                                                                     "Donor",
                                                                     "Flags"),
                                     to = c("Regulator",
                                            "Expr",
                                            "Type",
                                            "State",
                                            "z",
                                            "p",
                                            "Molecules",
                                            "Network",
                                            "Donor",
                                            "Flags"))
  datalist.final <- datalist.df
  if(!is.null(p.cutoff)){
    datalist.final <- datalist.final[datalist.final$p<p.cutoff,]
  }
  if(z.cutoff){
    datalist.final <- datalist.final[!is.na(datalist.final$z),]
    datalist.final <- datalist.final[!datalist.final$z==0,]
  }
  datalist.final[is.na(datalist.final$z),"z"] <- 0
  datalist.mat <- reshape2::dcast(data = datalist.final, formula = Regulator~Donor,fun.aggregate = sum, value.var = by) %>% column_to_rownames(var = "Regulator")
  if(by=="p") {
    datalist.mat <- -log(datalist.mat,10)
  }
  if(!is.null(top.n)) {
    top <- names((abs(rowSums(datalist.mat)) %>% sort(decreasing = T))[1:top.n])
    datalist.mat <- datalist.mat[rownames(datalist.mat) %in% top,]
  }
  return(datalist.mat)
}

`%notin%` <- Negate(`%in%`)

update.metadata <- function(seu = NULL, parent = covid_combined) {
  orig_meta <- seu@meta.data
  seu@meta.data <- parent@meta.data[colnames(parent) %in% colnames(seu),]
  seu@meta.data <- cbind(seu@meta.data, orig_meta[,colnames(orig_meta) %notin% colnames(parent@meta.data)])
  return(seu)
}



saveit <- function(format = "png", height = 7, width = 7, path = "~/Downloads/", name = "p") {
  if(format=="png") {
    ggsave(paste0(name,".png"), path = path, height = height, width = width, dpi = "retina")
  }
  if(format=="pdf") {
    ggsave(paste0(name,".pdf"), path = path, height = height, width = width)
  }
}






new.markers.heatmap.batch <- function(markers.matrix = NULL, fontsize_row = 10, color.rows = T,
                                      title = "log(FC)", add.flag = T, de.cutoff = 4, percentile = 0.999,
                                      top.n = NULL, color = c("blue", "white", "red"), 
                                      repel.degree = 0, legend = T, annotation_legend = T,
                                      cellwidth = NA, cellheight = NA, markers.label = NULL,
                                      save = T, width = 7, height = 11, cluster_columns = T,
                                      file = "~/Downloads/p.pdf") {
  row_annotation <- unique(row_annotation.batch[row_annotation.batch$Donor.full %in% colnames(markers.matrix),])
  row_annotation <- row_annotation[match(colnames(markers.matrix),row_annotation$Donor.full),]
  ta = columnAnnotation(Age = row_annotation$Age, 
                        MaxSeverity = row_annotation$Severity.final.max,
                        CurrentSeverity = row_annotation$Severity.final.current,
                        Gender = row_annotation$Sex,
                        Acuity = row_annotation$acuity,
                        DPT = row_annotation$days_post_pos_test,
                        Batch = row_annotation$Batch,
                        col = list(MaxSeverity = c("0" = severity.cols[1],
                                                   "1-3" = severity.cols[2],
                                                   "4-5" = severity.cols[3],
                                                   "6-7" = severity.cols[4],
                                                   "8" = severity.cols[5]),
                                   CurrentSeverity = c("0" = severity.cols[1],
                                                       "1-3" = severity.cols[2],
                                                       "4-5" = severity.cols[3],
                                                       "6-7" = severity.cols[4]),
                                   Age = colorRamp2(c(20,100),c("white","purple")),
                                   Gender = c("M" = "brown", "F" = "blue"),
                                   Acuity = c("Acute" = "orange", "Convalescent" = "purple", "Healthy" = "steelblue"),
                                   DPT = colorRamp2(c(0, 80),c("white","yellowgreen")),
                                   Batch = c("1" = hue_pal()(6)[1],
                                             "2" = hue_pal()(6)[2],
                                             "3" = hue_pal()(6)[3],
                                             "4" = hue_pal()(6)[4],
                                             "5" = hue_pal()(6)[5],
                                             "6" = hue_pal()(6)[6])))
  if(add.flag){
    if(!is.null(top.n)) {
      top.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
      kept.labels <- rownames(markers.matrix)[rownames(markers.matrix) %in% top.labels]
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    if(!is.null(markers.label)) {
      kept.labels <- markers.label[order(markers.label)]
    }
    t = 1:nrow(markers.matrix)
    
    ha = rowAnnotation(foo = anno_mark(at = t[rownames(markers.matrix) %in% kept.labels], 
                                       labels = kept.labels,
                                       labels_gp = gpar(col =
                                                          ifelse((rowSums(markers.matrix[rownames(markers.matrix) %in% kept.labels,]))>0,
                                                                 "red", "blue"))))
  }
  
  markers.matrix <- as.matrix(markers.matrix)
  col.lim = quantile(as.matrix(abs(markers.matrix)), probs = percentile, na.rm=T)
  p <- Heatmap(markers.matrix, name = title,
               right_annotation = ha, 
               top_annotation = ta, 
               row_names_gp = gpar(fontsize = 0),
               cluster_columns = cluster_columns,
               col = colorRamp2(breaks = c(-col.lim,0,col.lim), colors = color),
               column_names_gp = gpar(fontsize=0))
  if(save){
    pdf(file = file, width = width, height = height)
    p
    #dev.off()
  }
  else {
    p
  }
}

new.markers.heatmap.limitmeta <- function(markers.matrix = NULL, fontsize_row = 10, color.rows = T,
                                          title = "log(FC)", add.flag = T, de.cutoff = 4, percentile = 0.999,
                                          top.n = NULL, color = c("blue", "white", "red"), 
                                          repel.degree = 0, legend = T, annotation_legend = T,
                                          cellwidth = NA, cellheight = NA, markers.label = NULL,
                                          save = T, width = 7, height = 11, cluster_columns = T,
                                          file = "~/Downloads/p.pdf") {
  row_annotation <- unique(row_annotation[row_annotation$Donor.full %in% colnames(markers.matrix),])
  row_annotation <- row_annotation[match(colnames(markers.matrix),row_annotation$Donor.full),]
  ta = columnAnnotation(Severity = row_annotation$Severity.final.max,
                        DaysPostTest = row_annotation$days_post_pos_test,
                        DaysToDeath = row_annotation$time_to_death,
                        col = list(Severity = c("0" = severity.cols[1],
                                                "1-3" = severity.cols[2],
                                                "4-5" = severity.cols[3],
                                                "6-7" = severity.cols[4],
                                                "8" = severity.cols[5]),
                                   DaysPostTest = colorRamp2(c(0, 80),c("white","yellowgreen")),
                                   DaysToDeath = colorRamp2(c(0, 70),c("red","white"))
                        ))
  if(add.flag){
    if(!is.null(top.n)) {
      top.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
      kept.labels <- rownames(markers.matrix)[rownames(markers.matrix) %in% top.labels]
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    if(!is.null(markers.label)) {
      kept.labels <- markers.label[order(markers.label)]
    }
    t = 1:nrow(markers.matrix)
    
    ha = rowAnnotation(foo = anno_mark(at = t[rownames(markers.matrix) %in% kept.labels], 
                                       labels = kept.labels,
                                       labels_gp = gpar(col =
                                                          ifelse((rowSums(markers.matrix[rownames(markers.matrix) %in% kept.labels,]))>0,
                                                                 "red", "blue"))))
  }
  
  markers.matrix <- as.matrix(markers.matrix)
  col.lim = quantile(as.matrix(abs(markers.matrix)), probs = percentile, na.rm=T)
  p <- Heatmap(markers.matrix, name = title,
               right_annotation = ha, 
               top_annotation = ta, 
               row_names_gp = gpar(fontsize = 0),
               cluster_columns = cluster_columns,
               col = colorRamp2(breaks = c(-col.lim,0,col.lim), colors = color),
               column_names_gp = gpar(fontsize=0))
  if(save){
    pdf(file = file, width = width, height = height)
    p
    #dev.off()
  }
  else {
    p
  }
}



FeatureBoxPlot <- function(seu, by = "features", assay = "RNA", slot = "data", exprs.function = "mean",
                           features = NULL, meta.include = NULL, group_by = "Severity.current.final", shape_by = "acuity", 
                           each.pt = "Donor", custom_fill_colors = severity.cols, color_by = "Severity.max.final", 
                           comparisons = my_comparisons.current, 
                           ncol = 4, label = "p.signif", select.idents = NULL, 
                           label.x = NA, pt.size = NA, facet.face = "bold.italic") 
{
  DefaultAssay(seu) <- assay
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by == "cell.type" & length(features)>1) {
    message("Cannot facet by two variables, faceting by features")
    by = "features"
  }
  if(is.null(color_by)){
    color_by = group_by
  }
  data.features = (FetchData(seu, vars = features, slot = slot))
  if(by=="cell.type") {
    df <- aggregate(data.features, by = list(seu@meta.data[,each.pt], seu$cell.type), FUN = exprs.function)
    colnames(df)[1] <- each.pt
    colnames(df)[2] <- "cell.type"   
  }
  else {
    df <- aggregate(data.features, by = list(seu@meta.data[,each.pt]), FUN = exprs.function)
    colnames(df)[1] <- each.pt
  }
  meta.include <- unique(c(each.pt,group_by,shape_by,color_by))
  ei = unique(seu[[]][,meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  colnames(df) <- gsub("impadt_","",colnames(df))
  if(length(features)>1){
    df <- melt(df, id = colnames(df)[!colnames(df) %in% features])
  }
  else{
    df$variable = features
    df$value = df[,features]
  }
  if(by=="cell.type"){
    df$variable = df$cell.type
  }
  else {
    df$cell.type <- NULL
    df <- unique(df)
  }
  ggplot(df, aes_string(y = "value", x = group_by)) + labs(x = "Severity score at draw", y = "Average expression") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major = element_blank(), 
                       strip.background = element_rect(fill = NA, color = NA), 
                       strip.text = element_text(face = facet.face),
                       axis.ticks.x = element_blank(), 
                       axis.text = element_text(color = "black"), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("variable", scales = "free_y", ncol = ncol) + 
    guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
    geom_point(size = 4, position = position_jitter(width = 0.25),
               aes_string(x = group_by, y = "value", color = color_by, shape = shape_by)) +
    scale_shape_manual(values = shapes) + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.25), aspect.ratio = 1,
          text = element_text(size = 20)) +
    scale_color_manual(values = custom_fill_colors) +
    scale_fill_manual(values = custom_fill_colors) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
    ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label) + 
    labs(color = "Peak severity\nscore", shape = "Acuity")
  
} 


plot_traj_single_gene <- function(gene = NULL, group.by = "Severity.max.final"){
  
  myScale <- function(x) rescale(x, to = c(0,1), from = c(quantile(x,.05),
                                                          quantile(x,.95)))
  expr <- data.frame(GetAssayData(devneut_alone.nc)[gene,])
  expr <- apply(expr,2,myScale)
  
  aframe <- data.frame((expr), pt = t, 
                       group.by = devneut_alone.nc@meta.data[,group.by])
  colnames(aframe) <- c(gene,"pt",group.by)
  #melted <- melt(data = aframe, id.vars = "pt", measure.vars = setdiff(colnames(aframe), 'pt'))
  
  ggplot(aframe, aes_string(x = "pt", y = gene, group = group.by, color = group.by)) + 
    geom_smooth(method = "loess") +
    ylab(paste("Relative expression (scaled)")) + xlab("Latent time") + theme_bw() + 
    scale_color_manual(values = severity.cols) + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = NA, color = NA), 
          strip.text = element_text(face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text = element_text(color = "black"), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(color = group.by)
}

my_comparisons.current <- list(c("0","1-3"),
                               c("0","4-5"),
                               c("0","6-7"))
severity.cols <- c("steelblue","green3","yellow3","red","black")


plot_traj_genes <- function(seu = devneut_alone.nc, genes = 'MPO'){
  
  myScale <- function(x) rescale(x, to = c(0,1), from = c(quantile(x,.05),
                                                          quantile(x,.95)))
  expr <- t(data.matrix(GetAssayData(seu,assay="SCT",slot="data")[genes,]))
  expr <- apply(expr,2,myScale)
  
  aframe <- data.frame((expr), pt = t)
  
  melted <- melt(data = aframe, id.vars = "pt", measure.vars = setdiff(colnames(aframe), 'pt'))
  
  ggplot(melted, aes(x = pt, y = value, group = variable, color = variable)) + geom_smooth(method = "loess") +
    ylab(paste("Relative expression (scaled)")) + xlab("Latent time") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                           panel.grid.major = element_blank(), 
                                                                                           strip.background = element_rect(fill = NA, color = NA), 
                                                                                           strip.text = element_text(face = "bold"),
                                                                                           axis.ticks.x = element_blank(), 
                                                                                           axis.text = element_text(color = "black"), 
                                                                                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                                                                           legend.text = element_text(face = "italic")) + labs(color = "Gene")
}



## https://stackoverflow.com/questions/12381117/add-header-to-file-created-by-write-csv
my.write <- function(x, file, header, f = write.csv, ...){
  # create and open the file connection
  datafile <- file(file, open = 'wt')
  # close on exit
  on.exit(close(datafile))
  # if a header is defined, write it to the file (@CarlWitthoft's suggestion)
  if(!missing(header)) writeLines(header,con=datafile)
  # write the file using the defined function and required addition arguments  
  f(x, datafile,...)
}



FeatureBoxPlot.stroke <- function(seu, by = "features", assay = "RNA", slot = "data", exprs.function = "mean",
                                  features = NULL, meta.include = NULL, group_by = "Severity.current.final", shape_by = "stroke", 
                                  each.pt = "Donor", color_by = "Severity.max.final", custom_fill_colors = severity.cols, 
                                  comparisons = my_comparisons.current, alpha_by = "ddim",
                                  ncol = 4, label = "p.signif", select.idents = NULL, 
                                  label.x = NA, pt.size = NA, facet.face = "bold.italic") 
{
  DefaultAssay(seu) <- assay
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by == "cell.type" & length(features)>1) {
    message("Cannot facet by two variables, faceting by features")
    by = "features"
  }
  if(is.null(color_by)){
    color_by = group_by
  }
  data.features = (FetchData(seu, vars = features, slot = slot))
  if(by=="cell.type") {
    df <- aggregate(data.features, by = list(seu@meta.data[,each.pt], seu$cell.type), FUN = exprs.function)
    colnames(df)[1] <- each.pt
    colnames(df)[2] <- "cell.type"   
  }
  else {
    df <- aggregate(data.features, by = list(seu@meta.data[,each.pt]), FUN = exprs.function)
    colnames(df)[1] <- each.pt
  }
  meta.include <- unique(c(each.pt,group_by,shape_by,color_by,alpha_by))
  ei = unique(seu[[]][,meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  colnames(df) <- gsub("impadt_","",colnames(df))
  if(length(features)>1){
    df <- melt(df, id = colnames(df)[!colnames(df) %in% features])
  }
  else{
    df$variable = features
    df$value = df[,features]
  }
  if(by=="cell.type"){
    df$variable = df$cell.type
  }
  else {
    df$cell.type <- NULL
    df <- unique(df)
  }
  if(!is.null(alpha_by)){
    df[,alpha_by] <- ifelse(is.na(df[,alpha_by]),0,df[,alpha_by])
  }
  ggplot(df, aes_string(y = "value", x = group_by)) + labs(x = "Severity score at draw", y = "Average expression") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major = element_blank(), 
                       strip.background = element_rect(fill = NA, color = NA), 
                       strip.text = element_text(face = facet.face),
                       axis.ticks.x = element_blank(), 
                       axis.text = element_text(color = "black"), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("variable", scales = "free_y", ncol = ncol) + 
    guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
    geom_point(size = 4, position = position_jitter(width = 0.25),
               aes_string(x = group_by, y = "value", color = color_by, shape = shape_by, alpha = alpha_by)) +
    scale_shape_manual(values = shapes) + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.25), aspect.ratio = 1,
          text = element_text(size = 20)) +
    scale_color_manual(values = custom_fill_colors) +
    scale_fill_manual(values = custom_fill_colors) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
    ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label) + 
    labs(color = "Peak severity\nscore", shape = "Stroke", alpha = "D-dimer")
  
} 

propPlot <- function(seu, by = "cell.type", meta.include = NULL, group_by = "Severity.current.final", shape_by = "acuity", 
                     each.pt = "Donor", custom_fill_colors = severity.cols, color_by = "Severity.max.final", 
                     comparisons = my_comparisons.current, 
                     ncol = 4, label = "p.signif", select.idents = NULL, 
                     label.x = NA, pt.size = NA) 
{
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[,shape_by]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if(is.null(color_by)){
    color_by = group_by
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu$cell.type, seu[,each.pt]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               each.pt))  
  }
  meta.include <- unique(c(each.pt,group_by,shape_by,color_by))
  ei = unique(seu[,meta.include])
  df <- merge(df, ei, by = each.pt)
  df <- cbind(df, null.group = paste("1"))
  df[,each.pt] <- as.factor(df[,each.pt])
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  
  ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = "WHO severity score at draw", y = "Average expression") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major = element_blank(), 
                       strip.background = element_rect(fill = NA, color = NA), 
                       strip.text = element_text(face = "bold"),
                       axis.ticks.x = element_blank(), 
                       axis.text = element_text(color = "black"), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("cell.type", scales = "free_y", ncol = ncol) + 
    guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
    geom_point(size = 4, position = position_jitter(width = 0.25),
               aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
    scale_shape_manual(values = shapes) + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.25), aspect.ratio = 1,
          text = element_text(size = 20)) +
    scale_color_manual(values = custom_fill_colors) +
    scale_fill_manual(values = custom_fill_colors) + 
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.1))) + 
    ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label, 
                               p.adjust.methods = "bonferroni",
                               method = "wilcox.test") + 
    labs(color = "Peak severity\nscore", shape = "Acuity")
  
} 