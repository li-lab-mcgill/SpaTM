## Plot Functions


# Plot ScatterPie
plotPieSCTE <- function(spe,
                        sample_col,
                        sample_id,
                        radius = 50){
  cur_spe <- spe[,colData(spe)[,sample_col] == sample_id]
  
  coords <- spatialCoords(cur_spe) %>%
    as.data.frame()
  colnames(coords) <- c('x','y')
  cell_prop <- theta(cur_spe) %>% 
    as.data.frame()
  
  input_data <- cbind(coords,cell_prop)

  
  
  ggplot() + 
    geom_scatterpie(aes(x,y,r = radius),
                    data = input_data,
                    cols = colnames(cell_prop),
                    color = NA) + coord_equal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    labs(title = paste("Slide: ",sample_id,"\n",ncol(spe)," Spots",sep = ""),
         y = 'Y Coordinates',
         x = 'X Coordinates',
         fill = 'Celltype')
}

#Plot Weights


#Plot Heatmaps

# Theta Heatmap

theta_heatmap <- function(scte,labels = NULL){
  
  viz_theta <- theta(scte)[order(scte$Celltype),]#,order(unique(sce$Celltype))]
  
  colnames(viz_theta) <- colnames(alphaPrior(scte))
  
  col_scale <- colorRamp2(c(min(viz_theta),max(viz_theta)), c("white", "red"))
  
  
  viz_theta <- viz_theta[,levels(scte$Celltype)]
  
  cell_anno = rowAnnotation(Celltype = scte$Celltype[order(scte$Celltype)])
  pheatmap(viz_theta,cluster_cols = F,cluster_rows = F,color = col_scale,
           heatmap_legend_param = list(
             title = "Topic Weights", at = c(round(min(viz_theta,2)), round(max(viz_theta),2))),
           show_rownames = FALSE,
           left_annotation = c(cell_anno),
           cellwidth = 12,
           use_raster = FALSE
  )
}