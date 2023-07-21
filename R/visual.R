#' 
#' @keywords internal
#' 
get_random_colours <- function(n = 10, category = "qual") {
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == category,]
  all.colours <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                               rownames(qual_col_pals)))
  if (n > length(all.colours) ) {
    return( sample(all.colours, size = n, replace = TRUE) )
  } else {
    return( sample(all.colours, size = n, replace = FALSE) )
  }
}



#' Generate multiple histograms for two values
#' 
#' @keywords internal
#' 
#' @import ggplot2
#' 
get_multi_hist <- function(df = NULL, ylab = "Frequency", 
                           binwidth = NULL, size = 15, fill = "#63707E", 
                           hjust = 0, path = NULL, font.scale = 1.0, 
                           strip.face = "plain") {
  
  # Parameters
  message ("Generating ", length(unique(df[, 1])), " histograms ...")
  
  
  # Plotting
  if (!is.null(binwidth)) {
    p <- ggplot(df, aes(x = df[, 2], fill = df[, 1])) +
      geom_histogram(binwidth = binwidth, color = "black", fill = fill)
  } else {
    p <- ggplot(df, aes(x = df[, 2], fill = df[, 1])) +
      geom_bar(color = "black", fill = fill)
  }
  p <- p + facet_wrap(~ df[, 1]) + ylab(label = ylab) + 
    theme(axis.text.x = element_text(colour = "black", family = "Arial", size = size), 
          axis.text.y = element_text(colour = "black", family = "Arial", size = size), 
          axis.title.y = element_text(colour = "black", family = "Arial", size = font.scale * size), 
          axis.title.x = element_blank(), 
          panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "none", 
          strip.text = element_text(colour = "black", family = "Arial", 
                                    size = font.scale * size, face = strip.face, hjust = hjust), 
          strip.background = element_blank(), 
          axis.line = element_line(colour = "black"))

  return(p)
}



#' Generate multiple histograms for :
#' 1. The number of enhancers linked to each gene
#' 2. The rank of enhancers linked to each gene
#' 
#' @import dplyr
#' @export
#' @rdname plot_stats_genes_and_enhs
#' 
#' @param obj The \code{Seurat} object to identify enhancer regulons (eRegulons)
#' @param en.regs The list of eRegulons
#' @return Return a \code{ggplot} object
#' 
plot_stats_genes_and_enhs <- function(obj = NULL, peak.assay = "ATAC", 
                                     en.regs = NULL, 
                                     binwidth = NULL, size = 15, fill = "#63707E", 
                           hjust = 0, font.scale = 1.0, 
                           strip.face = "plain") {
  
  # Calculate the number of enhancers linked to each gene
  n.enhs <- unlist(pbmcapply::pbmclapply(en.regs, mc.cores = parallel::detectCores(), 
                                  function(x) {
                                    table(x$links$gene)
                                  }) )
  message ("Calculated the number of enhancers linked to each of the ", length(n.enhs), 
    " genes.")
  
  
  # Load gene coordinate annotations
  gene.coords <- CollapseToLongestTranscript(ranges = Signac::Annotation(object = obj[[peak.assay]]))
  gene.coords <- gene.coords[gene.coords$gene_name %in% Reduce("union", sapply(en.regs, "[[", "genes"))]
  tss <- GenomicRanges::resize(x = gene.coords, width = 1, fix = 'start') # find the TSS location
  peaks <- Signac::StringToGRanges(rownames(obj[[peak.assay]])) # get peak coordinates
  peaks <- peaks[Signac::GRangesToString(peaks) %in% Reduce("union", sapply(en.regs, "[[", "peaks"))]
  summ <- Matrix::summary(as(gUtils::gr.dist(gr1 = peaks, gr2 = tss), "sparseMatrix") )
  summ$i <- Signac::GRangesToString(peaks)[summ$i]
  summ$j <- tss$gene_name[summ$j]
  summ <- dplyr::filter(summ, is.not.na(x))
  
  
  # Calculate the rank of enhancers among all that linked to a gene
  rank.enhs <- unlist(pbmcapply::pbmclapply(en.regs, mc.cores = max(1, parallel::detectCores() / 3),
                                     function(x) {
                                       unlist(sapply(split(x$links, x$links$gene), function(y) {
                                         summ[summ$i %in% Signac::GRangesToString(y) &
                                              summ$j == unique(y$gene), "x"] %>% order
                                       }))
                                     }))
  message ("Calculated the rank of enhancers in increasing order of distance to each gene.")
  
  
  # Generate plot
  df <- as.data.frame(
    rbind(
        cbind(
            rep("# of enhs. per gene", length(n.enhs)),
            n.enhs
        ),
        cbind(
            rep("nth. gene per enh", length(rank.enhs)),
            rank.enhs
        )
    ))
  rownames(df) <- NULL
  return( get_multi_hist(df = df, ylab = "Frequency", 
                           binwidth = binwidth, size = size, fill = fill, 
                           hjust = hjust, font.scale = font.scale, 
                           strip.face = strip.face) )
}


#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
#' 
get_dotplot_heatmap <- function(df = NULL, 
                                ht.colours = NULL, 
                                ht.fill = "EX z-score", 
                                ht.breaks = NULL,
                                dot.name = "AC level", 
                                font.size = 15,
                                text.font = "Arial", 
                                path = NULL, angle = 30, 
                                row.anno = NULL,
                                vjust = 1, hjust = 1
                                ) {
  
  # Parameters
  message ("Preparing to generate dortplot-heatmap for ", 
           length(unique(df[, 1])), " rows and ", length(unique(df[, 2])), 
           " columns.")
  
  
  # Prepare sizes for dotplots
  dot.breaks <- cut(df[, 4], breaks = quantile(df[, 4]), include.lowest = TRUE)
  levels(dot.breaks) <- 1:4
  df <- cbind(df, dot.breaks)
  message ("Calculated the dot sizes in dotplot.")
  
  
  # Generate heatmap
  if (is.null(ht.breaks)) {
    ht.breaks <- trunc(range(df[, 3]) * 10) / 10
  }
  p <- ggplot(df, aes(y = ordered(df[, 1], levels = rev(levels(df[, 1]))),  
                      x = df[, 2])) +  
    geom_tile(aes(fill = df[, 3])) + 
    scale_fill_gradientn(colours = ht.colours, breaks = ht.breaks) + 
    geom_point(aes(size = df[, 5])) + 
    scale_size_manual(values = setNames(as.numeric(levels(df[, 5])), 
                                        levels(df[, 5]))) + 
    labs(fill = ht.fill, size = dot.name) + 
    theme(legend.title = element_text(size = font.size, family = text.font), 
          legend.text = element_text(size = font.size, family = text.font), 
          legend.key = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text.y = element_text(family = text.font, colour = "black", size = font.size), 
          axis.text.x = element_text(family = text.font, colour = "black", size = font.size, 
                                     angle = angle, vjust = vjust, hjust = hjust))
  
  
  # Split rows 
  if (!is.null(row.anno)) {
    message ("Annotating rows ...")
    right.anno.df <- data.frame(columnv = seq_along(row.anno), rowv = rep(1, length(row.anno)), 
                                value = rev(names(row.anno)))
    p.right.anno <- ggplot(right.anno.df, aes(x = rowv, y = columnv, fill = value)) + 
      geom_tile() + scale_fill_manual(values = row.anno) + 
      scale_y_discrete(position = "right") + 
      theme_minimal() + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank()) + 
      xlab(NULL) + ylab(NULL) + 
      theme(legend.title = element_blank(), 
            legend.text = element_blank(), 
            legend.key = element_blank(), 
            legend.position = "none",
            panel.background = element_blank(), 
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(), 
            axis.text.x = element_blank())
    # require(aplot)
    p <- ggplotify::as.ggplot(aplot::insert_right(p, p.right.anno, width = .05) )
  }
  return(p)
}



#' Get the dotplot-heatmap for cell-type-specific enhancer regulons (eRegulons)
#' 
#' @rdname plot_dotplot_heatmap
#' @export
#' 
#' @import dplyr
#' 
#' @param obj A \code{Seurat} object, NULL by default
#' @param cts.en.regs The list of cell-type-specific enhancer regulons (eRegulons)
#' @param n.regs The number of top-ranked eRegulons for plotting, 5 by default
#' @param dot What dot sizes represent, "AC" by default
#' @param peak.assay The assay denoting chromatin accessibility, "ATAC" by default
#' @param celltypes The list of cell types arranged in an order defined by users
#' @param celltype.col The metadata column indicating cell types, "seurat_clusters" by default
#' @param ht.colours Color list used to define the gradients of the heatmap, 
#' c(
#' grDevices::colorRampPalette(c("#99D5FF", "#E8F9FD"))(3),
#' "#fffaf7",
#' grDevices::colorRampPalette(c("#ffe4e1", "#ee7261"))(3)
#' )
#' by default
#' @param ht.breaks The breaks marked on the legend of heatmap colors
#' @param font.size The font size of legends
#' @param text.font The text font of legends
#' @param angle The angle of labels of the x-axis in the heatmap, 30 by default
#' @param vjust Adjustment of positions in the vertical direction, 1 by default 
#' @param hjust Adjustment of positions in the horizontal direction, 1 by default
#' @param celltype.colours The list of colors of cell types in the row annotations
#' 
plot_dotplot_heatmap <- function(obj = NULL, cts.en.regs = NULL,
                                 n.regs = 5, celltype.col = "seurat_clusters",
                                 dot = "AC", peak.assay = "ATAC",
                                 ht.colours = c(
                                   grDevices::colorRampPalette(c("#99D5FF", "#E8F9FD"))(3),
                                   "#fffaf7",
                                   grDevices::colorRampPalette(c("#ffe4e1", "#ee7261"))(3)
                                 ),
                                 ht.breaks = NULL,
                                 font.size = 15,
                                 text.font = "Arial",
                                 angle = 30,
                                 celltypes = NULL, 
                                 vjust = 1,
                                 hjust = 1, 
                                 celltype.colours = NULL
                                 ) {
  
  # Calculate the matrix for heatmap
  heat.lst <- pbmcapply::pbmclapply(names(celltypes), 
                                                   mc.cores = parallel::detectCores(),
                                                   function(x) {
                                                     x.regs <- cts.en.regs[sapply(cts.en.regs, 
                                                                                  "[[", "celltype") == x]
                                                     x.heat <- Reduce("rbind", lapply(x.regs, 
                                                                            function(y) {
                                                       Seurat::AverageExpression(obj, assays = "RNA", 
                                                                                 features = y$genes, 
                                                                                 group.by = celltype.col)$RNA %>% 
                                                                                apply(., 1, scale) %>% t %>% 
                                                                                apply(., 2, mean)
                                                     }) )
                                                     if (is.vector(x.heat)) {
                                                       x.heat <- matrix(x.heat, nrow = 1)
                                                     }
                                                     rownames(x.heat) <- which(sapply(cts.en.regs, 
                                                                                      "[[", "celltype") == x)
                                                     x.heat
                                                   })
  message ("Calculated the matrix for heatmap plotting.")
  
  
  # Select the top-ranked eRegulons for each cell type
  heat.m <- do.call("rbind", pbmcapply::pbmclapply(seq_along(heat.lst), mc.cores = parallel::detectCores(), 
                        function(i) {
                          score.lst <- heat.lst[[i]][, i, drop = FALSE] / 
                            (apply(heat.lst[[i]][, -i, drop = FALSE], 1, mean) + 1e-03)
                          head(heat.lst[[i]][order(score.lst, decreasing = TRUE), , drop = FALSE], 
                               n = n.regs)
                        }) )
  message ("Filtered the eRegulons by retaining the top-ranked ", n.regs, " ones for each of the ", 
           length(celltypes), " cell types,\n", 
           "obtaining a ", nrow(heat.m), " x ", ncol(heat.m), " matrix.")
  
  
  # Construct matrix for dotplot
  if (dot == "AC") {
    message ("Dotplot is used to visualize chromatin accessibility ...")
    dot.m <- do.call("rbind", pbmcapply::pbmclapply(rownames(heat.m), mc.cores = parallel::detectCores(), 
                                                    function(i) {
                                                      Seurat::AverageExpression(obj, assays = peak.assay, 
                                                                                features = cts.en.regs[[as.numeric(i)]]$enhancers, 
                                                                                group.by = 
                                                                                  celltype.col)[[peak.assay]] %>% 
                                                        apply(., 1, scale) %>% t %>%
                                                        apply(., 2, mean)
                                                    }) )
  } else {
    message ("Dotplot is used to visualize TF expression ...")
    dot.m <- do.call("rbind", pbmcapply::pbmclapply(rownames(heat.m), mc.cores = parallel::detectCores(), 
                                                    function(i) {
                                                      Seurat::AverageExpression(obj, assays = "RNA", 
                                                                                features = cts.en.regs[[as.numeric(i)]]$TF, 
                                                                                group.by = 
                                                                                  celltype.col)$RNA %>% 
                                                        apply(., 1, scale) %>% t %>%
                                                        apply(., 2, mean)
                                                    }) )
  }
  rownames(dot.m) <- rownames(heat.m)
  

  heat.summ <- Matrix::summary(as(heat.m, "sparseMatrix"))
  heat.summ <- cbind(heat.summ, sapply(1:nrow(heat.summ), function(i) {
    dot.m[heat.summ[i, 1], heat.summ[i, 2]]
  }))
  heat.summ$i <- paste0("eR", rownames(heat.m),
                        "_",
                        sapply(cts.en.regs[as.numeric(rownames(heat.m))],
                               "[[", "TF")[heat.summ$i])
  heat.summ$j <- celltypes[heat.summ$j]
  heat.summ$i <- factor(heat.summ$i, levels = unique(heat.summ$i))
  heat.summ$j <- factor(heat.summ$j, levels = unique(heat.summ$j))
  colnames(heat.summ) <- c("rowv", "columnv", "EX z-score", "AC level")
  if (dot == "TF") {
    dot.name <- "TF z-score"
  } else {
    dot.name <- "AC level"
  }
  
  
  # Set row annotations
  row.anno <- celltype.colours[sapply(cts.en.regs[as.numeric(rownames(heat.m))],
                     "[[", "celltype")]
  message ("Prepared row annotation recrtangles.")
  
  
  # Generate plot
  p <- get_dotplot_heatmap(df = heat.summ, ht.colours = ht.colours, 
                           ht.fill = "EX z-score", 
                           ht.breaks = ht.breaks,
                           dot.name = dot.name, 
                           font.size = font.size,
                           text.font = text.font, 
                           angle = angle, 
                           vjust = vjust, 
                           hjust = hjust, 
                           row.anno = row.anno)
  message ("Returning the dotplot-heatmap plot ...")
  return(p)
}



#' 
#' @import dplyr
#' @import igraph
#' @import rTRM
#' @import RColorBrewer
#' 
#' @keywords internal
#' 
df_to_igraph_net <- function(df = NULL, color.df = NULL, vertex.frame.width = 0.3, 
                             vertex.frame.color = NULL, tf.shape = "circle", 
                             enh.shape = "square", gene.shape = "circle", 
                             edge.curved = 0.0, path = "./", title = "temp", 
                             format = ".png",
                             width = 4000, height = 4000, res = 300,
                             layout = "fr",
                             vertex.label.family = "Arial", vertex.label.color = "black", 
                             vertex.size = NULL, tf.cex = 6, enh.cex = 1, gene.cex = 2, 
                             ifSaveRda = FALSE,
                             max.width = 5 ) {
  
  # Parameters
  message ("The data frame contains ", length(unique(df$TF)), 
           " TFs, ", length(unique(unique(df$enhancer))), " enhancers, and ", 
           length(unique(df$gene)), " genes.")
  if (length(setdiff(unique(df$TF), names(color.df))) > 0 | 
      is.null(color.df)) {
    message ("There are ", length(unique(df$TF)), " TFs but only ",
          length(unique(names(color.df))), " colors are given.", 
          "Using random color list ...")
    color.df <- setNames(get_random_colours(n = length(unique(df$TF))), 
                        unique(df$TF))
  }
  
  
  TR.links <- distinct(df[, c(1:2, 4)]) %>% dplyr::rename(., c(source = TF, target = enhancer, 
                                                               weight = T2R))
  RG.links <- distinct(df[, c(2:3, 5)]) %>% dplyr::rename(., c(source = enhancer, target = gene, 
                                                               weight = R2G))
  links <- rbind(TR.links, RG.links)
  links$weight <- max.width * (links$weight - min(links$weight) + 1e-04 ) / 
    (max(links$weight) - min(links$weight) + 1e-04)
  nodes <- apply(df[, 1:3], 2, unique) %>% unlist %>% unique
  message ("There are ", nrow(links), " linkages and ", length(nodes), " nodes.")
  
  
  # Function for plotting an elliptical node
  myellipse <- function(coords, v = NULL, params) {
    require(plotrix)
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1 / 30 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    
    plotrix::draw.ellipse(x = coords[,1], y = coords[,2],
                 a = vertex.size, b = vertex.size / 2, col = vertex.color)
  }
  
  
  # Register the shape with igraph
  igraph::add_shape("ellipse", clip = shapes("circle")$clip,
            plot = myellipse)
  
  
  # Directly save qsave files
  network <- graph_from_data_frame(d = links, vertices = nodes, directed = FALSE )
  if (ifSaveRda) {
    message ("Save the igraph object to file: ", 
             path, title, ".qsave.")
    qs::qsave(network, paste0(path, title, ".qsave"))
  }
  
  
  # Colors
  tf.color <- color.df[unique(df$TF)]
  tf.shape <- rep(tf.shape, length(unique(df$TF)))
  enh.shape <- rep(enh.shape, length(unique(df$enhancer)))
  gene.shape <- rep(gene.shape, length(unique(df$gene)))
  vertex.shape <- c(tf.shape, enh.shape, gene.shape)
  enh.color <- setNames(rep("#928A97", length(unique(df$enhancer))), 
                        unique(df$enhancer))
  gene.color <- setNames(rep("#FDB44B", length(unique(df$gene))), 
                         unique(df$gene))
  vertex.color <- c(tf.color, enh.color, gene.color)
  vertex.label <- c(names(tf.color), rep(NA, length(enh.color)), 
                    names(gene.color))
  TR.color <- tf.color[TR.links$source]
  RG.color <- rep("#5F6769", nrow(RG.links))
  edge.color <- c(TR.color, RG.color)
  if (is.null(vertex.size)) {
    tf.size <- setNames(rep(tf.cex, length(unique(df$TF))), 
                        unique(df$TF))
    enh.size <- setNames(rep(enh.cex, length(unique(df$enhancer))), 
                         unique(df$enhancer))
    gene.size <- setNames(rep(gene.cex, length(unique(df$gene))), 
                          unique(df$gene))
    vertex.size <- c(tf.size, enh.size, gene.size)
  }
  if (is.null(vertex.frame.color)) {
    vertex.frame.color <- vertex.color
  }
  if (format == ".tif") {
    tiff(filename = paste0(path, title, format), width = width, 
        height = height, res = res)
  } else {
    png(filename = paste0(path, title, format), width = width,
        height = height, res = res)
  }
  if (layout == "concentric") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = rTRM::layout.concentric(network, concentric = list(names(tf.color), 
                                                                     names(enh.color), 
                                                                     names(gene.color))) )
  } else if (layout == "fr") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_fr(network,
                                         minx = rep(-Inf, igraph::vcount(network)),
                                         maxx = rep(Inf, igraph::vcount(network)),
                                         miny = rep(-Inf, igraph::vcount(network)),
                                         maxy = rep(Inf, igraph::vcount(network)) )
    )
  } else if (layout == "lgl") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_lgl
    )
  } else if (layout == "mds") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_mds )
  } else if (layout == "kk") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_kk )
  } else if (layout == "dh") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_dh )
  } else if (layout == "nicely") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_nicely )
  } else if (layout == "drl") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_drl )
  } else if (layout == "gem") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_gem )
  } else if (layout == "graphopt") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = igraph::layout_with_graphopt )
  } else if (layout == "qgraph.fr") {
    plot(network, vertex.color = vertex.color, vertex.frame.width = vertex.frame.width, 
         vertex.frame.color = vertex.frame.color, edge.curved = edge.curved, edge.color = edge.color, 
         vertex.label.family = vertex.label.family, vertex.label.color = vertex.label.color, 
         vertex.shape = vertex.shape, vertex.size = vertex.size, vertex.label = vertex.label,
         layout = qgraph::qgraph.layout.fruchtermanreingold(get.edgelist(network, names = FALSE), 
                                                            vcount = vcount(network),
                                                    area = 8 * (vcount(network) ^ 2), 
                                                    repulse.rad = (vcount(network) ^ 3.1)) )
  }
  dev.off()
}



#' Plot an enhancer-driven gene regulatory network (eGRN)
#' 
#' @keywords internal
#' 
#' @export
#' 
#' @import dplyr
#' 
#' @param en.grn The enhancer gene regulatory network (eGRN) composed of \code{GRanges}, cell type, and cells
#' @param peak.assay The scATAC-seq assay
#' @param obj A \code{Seurat} object
#' @param n.links The cutoff of the number of TF-enhancer-gene relations for each TF to showcase.
#' We retain the top-ranked relations in decreasing order of the number of cells, where the genes including the 
#' ones encoding the TFs are expressed, while the enhancers are accessible, 50 by default
#' @param title The title of the eGRN to be plotted, "temp" by default
#' @param path The path to save the image file, "./" by default
#' @param format The format of the image file, ".png" by default
#' @param n.TFs The maximum number of top-ranked TFs to generate eGRN plot, 10 by default
#' @param TFs The list of TFs for plotting
#' @param layout The layout for plotting networks, "fr" by default
#' @param ifSaveRda Whether save the \code{igraph} objects in a qsdave file
#' 
plot_eGRN <- function(en.grn = NULL, obj = NULL, peak.assay = "ATAC", n.links = 20, 
                      title = "temp", path = "./", format = ".png", n.TFs = 15, 
                      layout = "fr", max.links = 5000, ifSaveRda = FALSE,
                      TFs = NULL) {
  
  # Parameters
  message (length(en.grn$links), " TF-enhancer-gene relations were input, \n", 
           "covering ", length(unique(en.grn$links$TF)), " TFs, ", 
           length(unique(Signac::GRangesToString(en.grn$links))), " enhancers, and ", 
           length(unique(en.grn$links$gene)), " genes.")
  message ("The input Seurat object contains ", nrow(obj[["RNA"]]), " genes and ", 
           nrow(obj[[peak.assay]]), " peaks.")

  
  # Filtering TFs, enhancers, and genes
  ranges <- en.grn$links
  ranges <- ranges[ranges$TF %in% rownames(obj[["RNA"]])]
  if (!is.null(TFs)) {
    ranges <- ranges[ranges$TF %in% TFs]
    if (length(ranges) < 1) {
      message ("There is no TF-enhancer-gene relation in the cell-type-specific eGRN!")
      return(NULL)
    }
    n.TFs <- Inf
    message ("Only showcase the eGRN composed of the ", 
             length(TFs), " TFs provided by the user.\n", 
             "Removed the restriction of the maximum number of TF-enhancer-gene relations for each TF.")
  }
  all.genes <- union(unique(ranges$TF), 
                     unique(ranges$gene) )
  all.enhs <- unique(Signac::GRangesToString(ranges))
  rna.m <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")[all.genes, en.grn$cells] > 0
  atac.m <- Seurat::GetAssayData(obj, assay = peak.assay, slot = "counts")[all.enhs, en.grn$cells] > 0
  message ("Filtered out the TFs whose encoding genes are not incorporated in the Seurart object,\n", 
           "leading to ", length(ranges), " relations.")
  
  
  # Filter relations for each TF
  if (length(unique(en.grn$links$TF)) + length(all.genes) + length(all.enhs) > max.links) {
    n.links <- round(max.links / length(unique(en.grn$links$TF)) )
    message ("Identifying the top-ranked ", n.links,
             " relations for each TF.")
  }
  ranges.lst <- pbmcapply::pbmclapply(unique(ranges$TF), 
                                                        mc.cores = parallel::detectCores(), 
                        function(x) {
                          x.ranges <- ranges[ranges$TF == x]
                          str.lst <- Signac::GRangesToString(x.ranges)
                          x.ranges[sapply(seq_along(x.ranges), function(i) {
                            intersect(
                              which(rna.m[x.ranges[i]$gene,] > 0 & rna.m[x.ranges[i]$TF,] > 0), 
                              which(atac.m[str.lst[i],] > 0)
                            ) %>% length
                          }) %>% order(., decreasing = TRUE) %>% head(., n = n.links)]
                        })
  filtered.ranges <- do.call("c", ranges.lst[sapply(ranges.lst, length) %>% 
                                               order(., decreasing = TRUE) %>% 
                                               head(., n = n.TFs)] )
  message ("Filtered GRanges to retain top-ranked ", length(filtered.ranges), 
           " TF-enhancer-gene relations.")
  
  
  # Plotting
  if (!dir.exists(path)) {
    dir.create(path)
    message("Created directory: ", path)
  }
  message ("Plotting the eGRN image to file: ", 
           path, title, format, " ...")
  df <- data.frame(
    TF = filtered.ranges$TF,
    enhancer = Signac::GRangesToString(filtered.ranges), 
    gene = filtered.ranges$gene
  )
  df <- cbind(df, 
              T2R = pbmcapply::pbmclapply(1:nrow(df), 
                                          mc.cores = parallel::detectCores(), 
                                          function(i) {
                                            length(which(rna.m[df$TF[i], ] > 0 & 
                                                           atac.m[df$enhancer[i], ] > 0))
                                          }) %>% unlist, 
              R2G = pbmcapply::pbmclapply(1:nrow(df), 
                                          mc.cores = parallel::detectCores(), 
                                          function(i) {
                                            length(which(rna.m[df$gene[i], ] > 0 & 
                                                           atac.m[df$enhancer[i], ] > 0))
                                          }) %>% unlist)
  df_to_igraph_net(
    df = df, 
    title = title, path = path, 
    format = format, 
    ifSaveRda = ifSaveRda,
    layout = layout
    )
}



#' Prepare TF-enhancer-gene relations to generate coverage plots
#' 
#' @export
#' @rdname prepare_coverage_plot
#' 
#' @param object A \code{Seurat} object
#' @param links A list of \code{GRanges} object to save the TF-enhancer-gene relations
#' @param peak.assay The scATAC-seq assay
#' @param highlight.width The width of the regions to highlight, i.e., TF binding sites
#' @param highlight.colors The colors of regions to highlight
#' 
prepare_coverage_plot <- function(object = NULL, links = NULL, 
                                  highlight.colors = NULL,
                      peak.assay = "ATAC", highlight.width = 500) {
  
  # Identify gene locations
  gene.coords <- CollapseToLongestTranscript(object[[peak.assay]]@annotation[object[[peak.assay]]@annotation$gene_name %in% unique(links$gene)])
  link.tss <- do.call("c", lapply(seq_along(links), function(x) {
    gene.coords[which(gene.coords$gene_name == links[x]$gene)]
  }))
  link.tss <- GenomicRanges::resize(link.tss, width = 1, fix = "start")
  message ("Identified the locations of genes.")
  
  
  # Determine the region
  intervals <- IRanges::punion(GenomicRanges::resize(links, width = 1, 
                                                     fix = "center"),
                               link.tss, fill.gap = TRUE, ignore.strand = TRUE)
  GenomicRanges::mcols(intervals)$gene <- link.tss$gene_name
  GenomicRanges::mcols(intervals)$peak <- Signac::GRangesToString(links)
  region <- paste0(
    unique(GenomeInfoDb::seqnames(intervals)), 
    "-",
    min(GenomicRanges::start(intervals)), 
    "-", 
    max(GenomicRanges::end(intervals))
  )
  message ("Determined the region to generate coverage plot.")
  
  
  # Calculate the regions to highlight
  region.highlight <- GenomicRanges::resize(Signac::StringToGRanges(unique(intervals$peak) ), 
                                             width = highlight.width, 
                                             fix = "center")
  if (is.null(highlight.colors)) {
    highlight.colors <- suppressWarnings(RColorBrewer::brewer.pal(length(unique(intervals$peak)), "Set1"))
  }
  highlight.colors <- head(highlight.colors, n = length(unique(intervals$peak) ))
  GenomicRanges::mcols(region.highlight)$color <- highlight.colors
  message ("Calculated the regions to highlight.")
  
  
  # Compute the scores of enhancer-gene relations
  rna.m <- Seurat::GetAssayData(object = object, assay = "RNA", slot = "counts") > 0
  atac.m <- Seurat::GetAssayData(object = object, assay = peak.assay, slot = "counts") > 0
  GenomicRanges::mcols(intervals)$score <- as.numeric(pbmcapply::pbmclapply(seq_along(intervals), 
                                                                 mc.cores = parallel::detectCores(), 
                                                                 function(x) {
                                                                   length(which(rna.m[intervals[x]$gene,] & 
                                                                     atac.m[intervals[x]$peak,] ) ) / 
                                                                     length(which(rna.m[intervals[x]$gene,] ) )
                                                                 }) )
  Signac::Links(object) <- intervals
  message ("Computed the score of enhancer-gene relations across all cells.")
  
  
  return(list(
    object = object,
    region = region,
    region.highlight = region.highlight
  ))
}



#' Generate coverage plot within an interval across different cell types
#' 
#' @export
#' @rdname get_coverage_plot
#' @import ggplot2
#' 
#' @param object A \code{Seurat} object
#' @param region A set of genomic coordinates to show. Can be a \code{GRanges} object, 
#' a string encoding a genomic position, a gene name, or a vector of strings describing 
#' the genomic coordinates or gene names to plot. If a gene name is supplied, annotations 
#' must be present in the assay.
#' @param features A vector of features present in another assay to plot alongside 
#' accessibility tracks (for example, gene names).
#' @param links Display links
#' @param ranges.group.by Grouping variable to color \code{ranges} by. Must be a variable present 
#' in the metadata stored in the \code{ranges} genomic ranges. If NULL, do not color by any variable.
#' @param peaks Display peaks
#' @param bigwig List of bigWig file paths to plot data from. Files can be remotely hosted. 
#' The name of each element in the list will determine the y-axis label given to the track.
#' @param region.highlight Region to highlight on the plot. Should be a \code{GRanges} object 
#' containing the coordinates to highlight. By default, regions will be highlighted in grey. 
#' To change the color of the highlighting, include a metadata column in the 
#' GRanges object named "color" containing the color to use for each region.
#' @param colors A list of colors to fill different genomic ranges
#' @param size The font size
#' @param idents Which identities to include in the plot. Default is all identities.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' 
get_coverage_plot <- function(object = NULL, region = NULL, 
                              features = NULL, links = FALSE, 
                              ranges.group.by = "seurat_cluster", peaks = FALSE, 
                              bigwig = NULL, region.highlight = NULL, 
                              idents = NULL, extend.upstream = 0, 
                              extend.downstream = 0,
                              colors = NULL, size = 10) {
  
  Seurat::DefaultAssay(object) <- peak.assay
  Seurat::Idents(object) <- setNames(object@meta.data[, ranges.group.by], 
                                     colnames(object))
  p <- suppressWarnings( Signac::CoveragePlot(object = object, region = region, features = features, 
                    ranges.group.by = ranges.group.by, peaks = peaks, links = links,
                    idents = idents,
                    bigwig = bigwig, region.highlight = region.highlight
  ) )
  if (is.not.null(colors)) {
    p <- p & 
      scale_fill_manual(values = colors)
  }
  # p <- p &  
  #   theme(
  #     text  = element_text(colour = "black", family = "Arial", size = size),
  #     strip.text.y = element_text(colour = "black", family = "Arial", size = size),
  #     axis.text.x = element_text(colour = "black", family = "Arial", size = size),
  #     axis.text.y = element_text(colour = "black", family = "Arial", size = size),
  #     axis.line.y = element_line(colour = "black"),
  #     panel.border = element_blank(),  
  #     panel.background = element_blank(), 
  #     legend.position = "none",
  #     axis.line.x = element_line(colour = "black"))
  
  
  return(p)
  }
