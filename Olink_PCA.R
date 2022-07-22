# Olink PCA Plotting Function Script

Olink_PCA <- function(path = path,
                      npx_data_type = npx_data_type,
                      npx_data_files = npx_data_files,
                      target_npx_data = NULL,
                      sample_grouping_info = sample_grouping_info,
                      folder_system = folder_system,
                      output_images_file_format = output_images_file_format,
                      qc_warning = qc_warning,
                      panel_pca = panel_pca,
                      firstcol = "UniProt",
                      firstrow = "SampleID",
                      cellvalue = "NPX",
                      sample_control_name = sample_control_name,
                      filtered_samples = filtered_samples) {

  ## 0.提前导入输出信息格式化设置R包
  ### BiocManager ###
  if (!requireNamespace("BiocManager", quietly = T)) {
    install.packages("BiocManager")
  }
  library(BiocManager)

  ## 1.安装和导入R包
  analibpackages <- c("rio", "crayon", "pcaMethods",
                      "ggplot2", "plotly", "FactoMineR",
                      "magrittr", "reshape2", "scales",
                      "openxlsx", "dplyr", "data.table",
                      "tibble", "parallel", "tidyr",
                      "purrr", "tidyverse", "stringr",
                      "ggsci", "ggpubr", "plotrix",
                      "cowplot", "showtext", "gridExtra",
                      "grid", "ggrepel", "plot3D",
                      "htmlwidgets")
  #
  for (a in analibpackages) {
    if (!requireNamespace(a, quietly = T)) {
      BiocManager::install(a)
    } else {
      # print(paste0("[软件包: '", a, "'已经安装完毕！]"))
      # print(packageVersion(a))
    }
  }
  #
  sapply(analibpackages, library, character.only = T)
  # options()
  options(bitmapType='cairo')

  ## 1.5 脚本运行开始
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(cyan(bold("--------------------------------[2.开始执行Olink项目OAAT2-QC-PCA: Olink_PCA()函数]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

  ## 2.根据函数参数获取Olink的NPX数据
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") { # rawtype or longtype

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one

      ### 2.1 Active R Function Rscript
      source("/glusterfs/home/local_jing_xx/software//OAAT2/Olink_NPX_Long2Wide.R", encoding = "utf-8")
      ### 2.2 NPX Long convert to Wide format
      if (npx_data_type == "rawtype") {
        olink_long2wide(raw_npx_data_path = paste0(path, "/", npx_data_files),
                        npx_data_path = NULL,
                        target_npx_data = NULL,
                        firstcol = firstcol,
                        firstrow = firstrow,
                        cellvalue = cellvalue,
                        filtered_samples = filtered_samples, # SC/Sample Control/CONTROL_SAMPLE_US_CS_AS_2
                        output_file_path = path)
      } else if (npx_data_type == "longtype") {
        olink_long2wide(raw_npx_data_path = NULL,
                        npx_data_path = paste0(path, "/", npx_data_files),
                        target_npx_data = NULL,
                        firstcol = firstcol,
                        firstrow = firstrow,
                        cellvalue = cellvalue,
                        filtered_samples = filtered_samples, # SC/Sample Control/CONTROL_SAMPLE_US_CS_AS_2
                        output_file_path = path)
      } else {
        print("### 2.2 NPX Long convert to Wide format")
      }

      ### 2.3 获取NPX宽数据
      npx_pca_file_index <- list.files(path = path, pattern = "ALL_Panel_NPX_WideData")
      npx_pca_df <- import(paste0(path, "/",npx_pca_file_index))
      # Long NPX Data
      npx_data <- import(paste0(path, "/", npx_data_files))

      ### 2.4 数据处理
      TMM <- data.frame(t(npx_pca_df))
      colnames(TMM) <- TMM[1,]
      TMM <- TMM[-c(1),]
      TMM2 <- as.data.frame(lapply(TMM, as.numeric))
      rownames(TMM2) <- rownames(TMM)
      res.pca <- PCA(TMM2, graph=F)
      dim1_lab <- paste0("PC1 (",round(res.pca$eig[1,2],1),"%)")
      dim2_lab <- paste0("PC2 (",round(res.pca$eig[2,2],1),"%)")
      npca <- ncol(res.pca$ind$coord)
      idx <- if (npca >= 3) 3 else 2
      if (idx == 3) {
        dim3_lab <- paste0("PC3 (",round(res.pca$eig[3,2],1),"%)")
      }
      #
      sgi <- import(paste0(path, "/", sample_grouping_info))
      ### 是否需要文件夹系统
      if (folder_system == T) { # 需要Olink分析文件夹系统
        # QC_Warning: 1.Warning/Pass; 2.WARN/PASS
        if (qc_warning == TRUE) { # qc_warning = TRUE

          ### 2.5 绘制2D PCA图
          coord <- as.data.frame(res.pca$ind$coord[,1:idx])
          coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
          coord$qc_warning <- sgi[match(rownames(coord), sgi$Sample),3]
          n.sample <- nrow(coord)
          n.group <- length(levels(as.factor(coord$group)))
          #### Save PCA Data ####
          rio::export(coord, file = paste0(path,"/1.QualityControl/1.3Samples_PCA/Olink_PCA_Data.xlsx"), overwrite = T, zoom = 120, rowNames = T)
          #### PC1 vs PC2 ####
          p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
            #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
            geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
            scale_colour_manual(values = rainbow(n.group))+
            guides(shape="none") +
            xlab(dim1_lab)+
            ylab(dim2_lab)+
            ggtitle("PC1 vs PC2")+
            xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
            ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
            theme_bw() +
            theme(axis.title=element_text(size=10),
                  legend.title=element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(size = guide_legend(),
                   shape = guide_legend()) # colour = guide_colorbar()

          if (n.sample < 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
          }
          ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC1_vs_PC2_QC_warning.", switch (output_images_file_format,
                                                                          "pdf" = "pdf",
                                                                          "png" = "png",
                                                                          "tiff" = "tiff",
                                                                          "jpeg" = "jpeg",
                                                                          "svg" = "svg",
                                                                          "bmp" = "bmp"
          )), width=10, height=8, units = "in")
          if (n.sample >= 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC1_vs_PC2_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
            )), width=10, height=8, units = "in")
          }
          #### PC2 vs PC3 ####
          if ( idx == 3 ) {
            p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
              #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
              geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
              scale_colour_manual(values = rainbow(n.group))+
              guides(shape="none") +
              xlab(dim2_lab)+
              ylab(dim3_lab)+
              ggtitle("PC2 vs PC3")+
              xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
              ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
              theme_bw() +
              theme(axis.title=element_text(size=10),
                    legend.title=element_blank(),
                    plot.title = element_text(hjust = 0.5)) +
              guides(size = guide_legend(),
                     shape = guide_legend()) # colour = guide_colorbar()

            if (n.sample < 9 ) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            }
            ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC2_vs_PC3_QC_warning.", switch (output_images_file_format,
                                                                            "pdf" = "pdf",
                                                                            "png" = "png",
                                                                            "tiff" = "tiff",
                                                                            "jpeg" = "jpeg",
                                                                            "svg" = "svg",
                                                                            "bmp" = "bmp"
            )), width=10, height=8, units = "in")

            if (n.sample >= 9) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC2_vs_PC3_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                               "pdf" = "pdf",
                                                                                               "png" = "png",
                                                                                               "tiff" = "tiff",
                                                                                               "jpeg" = "jpeg",
                                                                                               "svg" = "svg",
                                                                                               "bmp" = "bmp"
              )), width=10, height=8, units = "in")
            }
          }

          ### 2.6 绘制3D PCA图
          sn <- NROW(TMM)
          color_list <- colors()
          sacol <- sample(color_list, sn)
          sl <- sn/10.0
          spos <- 1:sn
          #### 保存图片格式
          if (output_images_file_format == "pdf") { #### PDF ####
            ##### scatter3D #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.pdf"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
            dev.off()

            ##### text3D #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.pdf"), width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "png") { #### PNG ####
            ##### scatter3D #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.png"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.png"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "tiff") { #### TIFF ####
            ##### scatter3D #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.tiff"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.tiff"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "jpeg") { #### JPEG ####
            ##### scatter3D #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.jpeg"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.jpeg"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "svg") { #### SVG ####
            ##### scatter3D #####
            svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.svg"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.svg"), width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "bmp") {#### BMP ####
            ##### scatter3D #####
            bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.bmp"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.bmp"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else {
            print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
          }

          ### 2.7 绘制PCA动图
          p <- plot_ly(coord,
                       text = ~rownames(coord),
                       x = ~Dim.1,
                       y = ~Dim.2,
                       z = ~Dim.3,
                       alpha = 0.8) %>%
            add_markers(color = ~group) %>%
            layout(scene = list(
                xaxis = list(title = dim1_lab),
                yaxis = list(title = dim2_lab),
                zaxis = list(title = dim3_lab)
                )
                )

          widget_file_size <- function(p) {
            d <- tempdir()
            withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
            f <- file.path(d, "index.html")
            mb <- round(file.info(f)$size / 1e6, 3)
            message("The 3D PCA plotly HTML File is: ", mb," MB")
          }

          widget_file_size(p)
          saveWidget(p, paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

        } else { # qc_warning = FALSE
          print("qc_warning == FALSE，不需要为PCA图添加QC_Warning注释！")

          ### 2.5 绘制2D PCA图
          coord <- as.data.frame(res.pca$ind$coord[,1:idx])
          coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
          n.sample <- nrow(coord)
          n.group <- length(levels(as.factor(coord$group)))
          #### PC1 vs PC2 ####
          p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
            #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
            geom_point(size = 3, aes(colour = factor(coord$group)))+
            scale_colour_manual(values = rainbow(n.group))+
            guides(shape="none") +
            xlab(dim1_lab)+
            ylab(dim2_lab)+
            ggtitle("PC1 vs PC2")+
            xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
            ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
            theme_bw() +
            theme(axis.title=element_text(size=10),
                  legend.title=element_blank(),
                  plot.title = element_text(hjust = 0.5))

          if (n.sample < 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
          }
          ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC1_vs_PC2.", switch (output_images_file_format,
                                                               "pdf" = "pdf",
                                                               "png" = "png",
                                                               "tiff" = "tiff",
                                                               "jpeg" = "jpeg",
                                                               "svg" = "svg",
                                                               "bmp" = "bmp"
          )), width=10, height=8, units = "in")
          if (n.sample >= 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC1_vs_PC2_SampleName_label.", switch (output_images_file_format,
                                                                                  "pdf" = "pdf",
                                                                                  "png" = "png",
                                                                                  "tiff" = "tiff",
                                                                                  "jpeg" = "jpeg",
                                                                                  "svg" = "svg",
                                                                                  "bmp" = "bmp"
            )), width=10, height=8, units = "in")
          }
          #### PC2 vs PC3 ####
          if ( idx == 3 ) {
            p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
              #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
              geom_point(size = 3, aes(colour = factor(coord$group)))+
              scale_colour_manual(values = rainbow(n.group))+
              guides(shape="none") +
              xlab(dim2_lab)+
              ylab(dim3_lab)+
              ggtitle("PC2 vs PC3")+
              xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
              ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
              theme_bw() +
              theme(axis.title=element_text(size=10),
                    legend.title=element_blank(),
                    plot.title = element_text(hjust = 0.5))

            if (n.sample < 9 ) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            }
            ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC2_vs_PC3.", switch (output_images_file_format,
                                                                 "pdf" = "pdf",
                                                                 "png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "svg" = "svg",
                                                                 "bmp" = "bmp"
            )), width=10, height=8, units = "in")

            if (n.sample >= 9) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_PC2_vs_PC3_SampleName_label.", switch (output_images_file_format,
                                                                                    "pdf" = "pdf",
                                                                                    "png" = "png",
                                                                                    "tiff" = "tiff",
                                                                                    "jpeg" = "jpeg",
                                                                                    "svg" = "svg",
                                                                                    "bmp" = "bmp"
              )), width=10, height=8, units = "in")
            }
          }

          ### 2.6 绘制3D PCA图
          sn <- NROW(TMM)
          color_list <- colors()
          sacol <- sample(color_list, sn)
          sl <- sn/10.0
          spos <- 1:sn
          #### 保存图片格式
          if (output_images_file_format == "pdf") { #### PDF ####
            ##### scatter3D #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.pdf"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
            dev.off()

            ##### text3D #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.pdf"), width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "png") { #### PNG ####
            ##### scatter3D #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.png"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.png"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            png(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "tiff") { #### TIFF ####
            ##### scatter3D #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.tiff"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.tiff"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "jpeg") { #### JPEG ####
            ##### scatter3D #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.jpeg"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.jpeg"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "svg") { #### SVG ####
            ##### scatter3D #####
            svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.svg"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            svg(paste0(path, "/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.svg"), width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            svg(path, paste0("/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "bmp") {#### BMP ####
            ##### scatter3D #####
            bmp(paste0(path, "/1.QualityControl/1.3Samples_PCA/PCA_3D_Point.bmp"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Text.bmp"), width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else {
            print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
          }

          ### 2.7 绘制PCA动图
          p <- plot_ly(coord,
                       text = ~rownames(coord),
                       x = ~Dim.1,
                       y = ~Dim.2,
                       z = ~Dim.3,
                       alpha = 0.8) %>%
          add_markers(color = ~group) %>%
            layout(scene = list(
              xaxis = list(title = dim1_lab),
              yaxis = list(title = dim2_lab),
              zaxis = list(title = dim3_lab)
            )
            )

          widget_file_size <- function(p) {
            d <- tempdir()
            withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
            f <- file.path(d, "index.html")
            mb <- round(file.info(f)$size / 1e6, 3)
            message("The 3D PCA plotly HTML File is: ", mb," MB")
          }

          widget_file_size(p)
          saveWidget(p, paste0(path, "/1.QualityControl/1.3Samples_PCA/PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

        }

      } else { # 不需要Olink分析文件夹系统
        # QC_Warning: 1.Warning/Pass; 2.WARN/PASS
        if (qc_warning == TRUE) { # qc_warning = TRUE

          ### 2.5 绘制2D PCA图
          coord <- as.data.frame(res.pca$ind$coord[,1:idx])
          coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
          coord$qc_warning <- sgi[match(rownames(coord), sgi$Sample),3]
          n.sample <- nrow(coord)
          n.group <- length(levels(as.factor(coord$group)))
          #### Save PCA Data ####
          rio::export(coord, file = paste0(path,"/Olink_PCA_Data.xlsx"), overwrite = T, zoom = 120, rowNames = T)
          #### PC1 vs PC2 ####
          p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
            #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
            geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
            scale_colour_manual(values = rainbow(n.group))+
            guides(shape="none") +
            xlab(dim1_lab)+
            ylab(dim2_lab)+
            ggtitle("PC1 vs PC2")+
            xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
            ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
            theme_bw() +
            theme(axis.title=element_text(size=10),
                  legend.title=element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(size = guide_legend(),
                   shape = guide_legend()) # colour = guide_colorbar()

          if (n.sample < 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
          }
          ggsave(p, filename=paste0("PCA_PC1_vs_PC2_QC_warning.", switch (output_images_file_format,
                                                                          "pdf" = "pdf",
                                                                          "png" = "png",
                                                                          "tiff" = "tiff",
                                                                          "jpeg" = "jpeg",
                                                                          "svg" = "svg",
                                                                          "bmp" = "bmp"
          )), width=10, height=8, units = "in")
          if (n.sample >= 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            ggsave(p, filename=paste0("PCA_PC1_vs_PC2_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
            )), width=10, height=8, units = "in")
          }
          #### PC2 vs PC3 ####
          if ( idx == 3 ) {
            p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
              #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
              geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
              scale_colour_manual(values = rainbow(n.group))+
              guides(shape="none") +
              xlab(dim2_lab)+
              ylab(dim3_lab)+
              ggtitle("PC2 vs PC3")+
              xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
              ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
              theme_bw() +
              theme(axis.title=element_text(size=10),
                    legend.title=element_blank(),
                    plot.title = element_text(hjust = 0.5)) +
              guides(size = guide_legend(),
                     shape = guide_legend()) # colour = guide_colorbar()

            if (n.sample < 9 ) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            }
            ggsave(p, filename=paste0("PCA_PC2_vs_PC3_QC_warning.", switch (output_images_file_format,
                                                                            "pdf" = "pdf",
                                                                            "png" = "png",
                                                                            "tiff" = "tiff",
                                                                            "jpeg" = "jpeg",
                                                                            "svg" = "svg",
                                                                            "bmp" = "bmp"
            )), width=10, height=8, units = "in")

            if (n.sample >= 9) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              ggsave(p, filename=paste0("PCA_PC2_vs_PC3_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                               "pdf" = "pdf",
                                                                                               "png" = "png",
                                                                                               "tiff" = "tiff",
                                                                                               "jpeg" = "jpeg",
                                                                                               "svg" = "svg",
                                                                                               "bmp" = "bmp"
              )), width=10, height=8, units = "in")
            }
          }

          ### 2.6 绘制3D PCA图
          sn <- NROW(TMM)
          color_list <- colors()
          sacol <- sample(color_list, sn)
          sl <- sn/10.0
          spos <- 1:sn
          #### 保存图片格式
          if (output_images_file_format == "pdf") { #### PDF ####
            ##### scatter3D #####
            pdf("PCA_3D_Point.pdf", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
            dev.off()

            ##### text3D #####
            pdf("PCA_3D_Text.pdf", width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            pdf("PCA_3D_Point_And_Text.pdf", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "png") { #### PNG ####
            ##### scatter3D #####
            png("PCA_3D_Point.png", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            png("PCA_3D_Text.png", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            png("PCA_3D_Point_And_Text.png", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "tiff") { #### TIFF ####
            ##### scatter3D #####
            tiff("PCA_3D_Point.tiff", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            tiff("PCA_3D_Text.tiff", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            tiff("PCA_3D_Point_And_Text.tiff", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "jpeg") { #### JPEG ####
            ##### scatter3D #####
            jpeg("PCA_3D_Point.jpeg", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            jpeg("PCA_3D_Text.jpeg", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            jpeg("PCA_3D_Point_And_Text.jpeg", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "svg") { #### SVG ####
            ##### scatter3D #####
            svg("PCA_3D_Point.svg", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            svg("PCA_3D_Text.svg", width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            svg("PCA_3D_Point_And_Text.svg", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "bmp") {#### BMP ####
            ##### scatter3D #####
            bmp("PCA_3D_Point.bmp", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            bmp("PCA_3D_Text.bmp", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            bmp("PCA_3D_Point_And_Text.bmp", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else {
            print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
          }

          ### 2.7 绘制PCA动图
          p <- plot_ly(coord,
                       text = ~rownames(coord),
                       x = ~Dim.1,
                       y = ~Dim.2,
                       z = ~Dim.3,
                       alpha = 0.8) %>%
            add_markers(color = ~group) %>%
            layout(scene = list(
              xaxis = list(title = dim1_lab),
              yaxis = list(title = dim2_lab),
              zaxis = list(title = dim3_lab)
            )
            )

          widget_file_size <- function(p) {
            d <- tempdir()
            withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
            f <- file.path(d, "index.html")
            mb <- round(file.info(f)$size / 1e6, 3)
            message("The 3D PCA plotly HTML File is: ", mb," MB")
          }

          widget_file_size(p)
          saveWidget(p, "PCA_3D_Plotly.html", selfcontained = F, libdir = "lib")

        } else { # qc_warning = FALSE
          print("qc_warning == FALSE，不需要为PCA图添加QC_Warning注释！")

          ### 2.5 绘制2D PCA图
          coord <- as.data.frame(res.pca$ind$coord[,1:idx])
          coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
          n.sample <- nrow(coord)
          n.group <- length(levels(as.factor(coord$group)))
          #### PC1 vs PC2 ####
          p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
            #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
            geom_point(size = 3, aes(colour = factor(coord$group)))+
            scale_colour_manual(values = rainbow(n.group))+
            guides(shape="none") +
            xlab(dim1_lab)+
            ylab(dim2_lab)+
            ggtitle("PC1 vs PC2")+
            xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
            ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
            theme_bw() +
            theme(axis.title=element_text(size=10),
                  legend.title=element_blank(),
                  plot.title = element_text(hjust = 0.5))

          if (n.sample < 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
          }
          ggsave(p, filename=paste0("PCA_PC1_vs_PC2.", switch (output_images_file_format,
                                                               "pdf" = "pdf",
                                                               "png" = "png",
                                                               "tiff" = "tiff",
                                                               "jpeg" = "jpeg",
                                                               "svg" = "svg",
                                                               "bmp" = "bmp"
          )), width=10, height=8, units = "in")
          if (n.sample >= 9) {
            p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            ggsave(p, filename=paste0("PCA_PC1_vs_PC2_SampleName_label.", switch (output_images_file_format,
                                                                                  "pdf" = "pdf",
                                                                                  "png" = "png",
                                                                                  "tiff" = "tiff",
                                                                                  "jpeg" = "jpeg",
                                                                                  "svg" = "svg",
                                                                                  "bmp" = "bmp"
            )), width=10, height=8, units = "in")
          }
          #### PC2 vs PC3 ####
          if ( idx == 3 ) {
            p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
              #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
              geom_point(size = 3, aes(colour = factor(coord$group)))+
              scale_colour_manual(values = rainbow(n.group))+
              guides(shape="none") +
              xlab(dim2_lab)+
              ylab(dim3_lab)+
              ggtitle("PC2 vs PC3")+
              xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
              ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
              theme_bw() +
              theme(axis.title=element_text(size=10),
                    legend.title=element_blank(),
                    plot.title = element_text(hjust = 0.5))

            if (n.sample < 9 ) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
            }
            ggsave(p, filename=paste0("PCA_PC2_vs_PC3.", switch (output_images_file_format,
                                                                 "pdf" = "pdf",
                                                                 "png" = "png",
                                                                 "tiff" = "tiff",
                                                                 "jpeg" = "jpeg",
                                                                 "svg" = "svg",
                                                                 "bmp" = "bmp"
            )), width=10, height=8, units = "in")

            if (n.sample >= 9) {
              p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              ggsave(p, filename=paste0("PCA_PC2_vs_PC3_SampleName_label.", switch (output_images_file_format,
                                                                                    "pdf" = "pdf",
                                                                                    "png" = "png",
                                                                                    "tiff" = "tiff",
                                                                                    "jpeg" = "jpeg",
                                                                                    "svg" = "svg",
                                                                                    "bmp" = "bmp"
              )), width=10, height=8, units = "in")
            }
          }

          ### 2.6 绘制3D PCA图
          sn <- NROW(TMM)
          color_list <- colors()
          sacol <- sample(color_list, sn)
          sl <- sn/10.0
          spos <- 1:sn
          #### 保存图片格式
          if (output_images_file_format == "pdf") { #### PDF ####
            ##### scatter3D #####
            pdf("PCA_3D_Point.pdf", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
            dev.off()

            ##### text3D #####
            pdf("PCA_3D_Text.pdf", width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            pdf("PCA_3D_Point_And_Text.pdf", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "png") { #### PNG ####
            ##### scatter3D #####
            png("PCA_3D_Point.png", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            png("PCA_3D_Text.png", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            png("PCA_3D_Point_And_Text.png", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "tiff") { #### TIFF ####
            ##### scatter3D #####
            tiff("PCA_3D_Point.tiff", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            tiff("PCA_3D_Text.tiff", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            tiff("PCA_3D_Point_And_Text.tiff", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else if (output_images_file_format == "jpeg") { #### JPEG ####
            ##### scatter3D #####
            jpeg("PCA_3D_Point.jpeg", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            jpeg("PCA_3D_Text.jpeg", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            jpeg("PCA_3D_Point_And_Text.jpeg", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "svg") { #### SVG ####
            ##### scatter3D #####
            svg("PCA_3D_Point.svg", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            svg("PCA_3D_Text.svg", width = 12, height = 12)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            svg("PCA_3D_Point_And_Text.svg", width = 12, height = 12)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #

          } else if (output_images_file_format == "bmp") {#### BMP ####
            ##### scatter3D #####
            bmp("PCA_3D_Point.bmp", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### text3D #####
            bmp("PCA_3D_Text.bmp", width = 1200, height = 1200)
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            dev.off()

            ##### Point and Text #####
            bmp("PCA_3D_Point_And_Text.bmp", width = 1200, height = 1200)
            scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
            #
            text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
            dev.off()
            #
          } else {
            print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
          }

          ### 2.7 绘制PCA动图
          p <- plot_ly(coord,
                       text = ~rownames(coord),
                       x = ~Dim.1,
                       y = ~Dim.2,
                       z = ~Dim.3,
                       alpha = 0.8) %>%
            add_markers(color = ~group) %>%
            layout(scene = list(
              xaxis = list(title = dim1_lab),
              yaxis = list(title = dim2_lab),
              zaxis = list(title = dim3_lab)
            )
            )

          widget_file_size <- function(p) {
            d <- tempdir()
            withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
            f <- file.path(d, "index.html")
            mb <- round(file.info(f)$size / 1e6, 3)
            message("The 3D PCA plotly HTML File is: ", mb," MB")
          }

          widget_file_size(p)
          saveWidget(p, "PCA_3D_Plotly.html", selfcontained = F, libdir = "lib")

        }

      }
      ##########################################################################
      ### Panel_PCA Plotting ###
      ##########################################################################
      if (panel_pca == T) { # panel_pca == T
        #
        pn <- 1 # 计数器初始化
        for (pl in unique(npx_data$Panel)) {
          cat(blue("[ ", pn, ".当前处理的Panel为："), bold(red(pl,"]\n")))
          ### 2.3 获取Panel NPX宽数据
          panel_npx_pca_file_index <- list.files(path = path, pattern = paste0("*",pl,"_Panel_WideData*"))[1]
          panel_npx_pca_df <- import(paste0(path, "/",panel_npx_pca_file_index))

          ### 2.4 Panel数据处理
          TMM <- data.frame(t(panel_npx_pca_df))
          colnames(TMM) <- TMM[1,]
          TMM <- TMM[-c(1),]
          TMM2 <- as.data.frame(lapply(TMM, as.numeric))
          rownames(TMM2) <- rownames(TMM)
          res.pca <- PCA(TMM2, graph=F) # Panel的res.pca
          dim1_lab <- paste0("PC1 (",round(res.pca$eig[1,2],1),"%)")
          dim2_lab <- paste0("PC2 (",round(res.pca$eig[2,2],1),"%)")
          npca <- ncol(res.pca$ind$coord)
          idx <- if (npca >= 3) 3 else 2
          if (idx == 3) {
            dim3_lab <- paste0("PC3 (",round(res.pca$eig[3,2],1),"%)")
          }
          #
          sgi <- import(paste0(path, "/", sample_grouping_info))
          #
          #
          #-----------------------------------------------------------------------
          ### 是否需要文件夹系统
          if (folder_system == T) { # 需要Olink分析文件夹系统
            # QC_Warning: 1.Warning/Pass; 2.WARN/PASS
            if (qc_warning == TRUE) { # qc_warning = TRUE

              ### 2.5 绘制2D PCA图
              coord <- as.data.frame(res.pca$ind$coord[,1:idx])
              coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
              coord$qc_warning <- sgi[match(rownames(coord), sgi$Sample),3]
              n.sample <- nrow(coord)
              n.group <- length(levels(as.factor(coord$group)))
              #### Save PCA Data ####
              rio::export(coord, file = paste0(path,"/1.QualityControl/1.3Samples_PCA/Olink_",pl,"_PCA_Data.xlsx"), overwrite = T, zoom = 120, rowNames = T)
              #### PC1 vs PC2 ####
              p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
                #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
                scale_colour_manual(values = rainbow(n.group))+
                guides(shape="none") +
                xlab(dim1_lab)+
                ylab(dim2_lab)+
                ggtitle("PC1 vs PC2")+
                xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
                ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                theme_bw() +
                theme(axis.title=element_text(size=10),
                      legend.title=element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(size = guide_legend(),
                       shape = guide_legend()) # colour = guide_colorbar()

              if (n.sample < 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              }
              ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/",pl,"_PCA_PC1_vs_PC2_QC_warning.", switch (output_images_file_format,
                                                                                                                           "pdf" = "pdf",
                                                                                                                           "png" = "png",
                                                                                                                           "tiff" = "tiff",
                                                                                                                           "jpeg" = "jpeg",
                                                                                                                           "svg" = "svg",
                                                                                                                           "bmp" = "bmp"
              )), width=10, height=8, units = "in")
              if (n.sample >= 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC1_vs_PC2_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                                                                                "pdf" = "pdf",
                                                                                                                                                "png" = "png",
                                                                                                                                                "tiff" = "tiff",
                                                                                                                                                "jpeg" = "jpeg",
                                                                                                                                                "svg" = "svg",
                                                                                                                                                "bmp" = "bmp"
                )), width=10, height=8, units = "in")
              }
              #### PC2 vs PC3 ####
              if ( idx == 3 ) {
                p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
                  #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                  geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
                  scale_colour_manual(values = rainbow(n.group))+
                  guides(shape="none") +
                  xlab(dim2_lab)+
                  ylab(dim3_lab)+
                  ggtitle("PC2 vs PC3")+
                  xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                  ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
                  theme_bw() +
                  theme(axis.title=element_text(size=10),
                        legend.title=element_blank(),
                        plot.title = element_text(hjust = 0.5)) +
                  guides(size = guide_legend(),
                         shape = guide_legend()) # colour = guide_colorbar()

                if (n.sample < 9 ) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                }
                ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC2_vs_PC3_QC_warning.", switch (output_images_file_format,
                                                                                                                               "pdf" = "pdf",
                                                                                                                               "png" = "png",
                                                                                                                               "tiff" = "tiff",
                                                                                                                               "jpeg" = "jpeg",
                                                                                                                               "svg" = "svg",
                                                                                                                               "bmp" = "bmp"
                )), width=10, height=8, units = "in")

                if (n.sample >= 9) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                  ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC2_vs_PC3_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                                                                                  "pdf" = "pdf",
                                                                                                                                                  "png" = "png",
                                                                                                                                                  "tiff" = "tiff",
                                                                                                                                                  "jpeg" = "jpeg",
                                                                                                                                                  "svg" = "svg",
                                                                                                                                                  "bmp" = "bmp"
                  )), width=10, height=8, units = "in")
                }
              }

              ### 2.6 绘制3D PCA图
              sn <- NROW(TMM)
              color_list <- colors()
              sacol <- sample(color_list, sn)
              sl <- sn/10.0
              spos <- 1:sn
              #### 保存图片格式
              if (output_images_file_format == "pdf") { #### PDF ####
                ##### scatter3D #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
                dev.off()

                ##### text3D #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.pdf"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "png") { #### PNG ####
                ##### scatter3D #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.png"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "tiff") { #### TIFF ####
                ##### scatter3D #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.tiff"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "jpeg") { #### JPEG ####
                ##### scatter3D #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.jpeg"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "svg") { #### SVG ####
                ##### scatter3D #####
                svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.svg"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "bmp") {#### BMP ####
                ##### scatter3D #####
                bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.bmp"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else {
                print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
              }

              ### 2.7 绘制PCA动图
              p <- plot_ly(coord,
                           text = ~rownames(coord),
                           x = ~Dim.1,
                           y = ~Dim.2,
                           z = ~Dim.3,
                           alpha = 0.8) %>%
                add_markers(color = ~group) %>%
                layout(scene = list(
                  xaxis = list(title = dim1_lab),
                  yaxis = list(title = dim2_lab),
                  zaxis = list(title = dim3_lab)
                )
                )

              widget_file_size <- function(p) {
                d <- tempdir()
                withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
                f <- file.path(d, "index.html")
                mb <- round(file.info(f)$size / 1e6, 3)
                message("The 3D PCA plotly HTML File is: ", mb," MB")
              }

              widget_file_size(p)
              saveWidget(p, paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

            } else { # qc_warning = FALSE
              print("qc_warning == FALSE，不需要为PCA图添加QC_Warning注释！")

              ### 2.5 绘制2D PCA图
              coord <- as.data.frame(res.pca$ind$coord[,1:idx])
              coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
              n.sample <- nrow(coord)
              n.group <- length(levels(as.factor(coord$group)))
              #### PC1 vs PC2 ####
              p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
                #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                geom_point(size = 3, aes(colour = factor(coord$group)))+
                scale_colour_manual(values = rainbow(n.group))+
                guides(shape="none") +
                xlab(dim1_lab)+
                ylab(dim2_lab)+
                ggtitle("PC1 vs PC2")+
                xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
                ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                theme_bw() +
                theme(axis.title=element_text(size=10),
                      legend.title=element_blank(),
                      plot.title = element_text(hjust = 0.5))

              if (n.sample < 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              }
              ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC1_vs_PC2.", switch (output_images_file_format,
                                                                                                                  "pdf" = "pdf",
                                                                                                                  "png" = "png",
                                                                                                                  "tiff" = "tiff",
                                                                                                                  "jpeg" = "jpeg",
                                                                                                                  "svg" = "svg",
                                                                                                                  "bmp" = "bmp"
              )), width=10, height=8, units = "in")
              if (n.sample >= 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC1_vs_PC2_SampleName_label.", switch (output_images_file_format,
                                                                                                                                     "pdf" = "pdf",
                                                                                                                                     "png" = "png",
                                                                                                                                     "tiff" = "tiff",
                                                                                                                                     "jpeg" = "jpeg",
                                                                                                                                     "svg" = "svg",
                                                                                                                                     "bmp" = "bmp"
                )), width=10, height=8, units = "in")
              }
              #### PC2 vs PC3 ####
              if ( idx == 3 ) {
                p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
                  #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                  geom_point(size = 3, aes(colour = factor(coord$group)))+
                  scale_colour_manual(values = rainbow(n.group))+
                  guides(shape="none") +
                  xlab(dim2_lab)+
                  ylab(dim3_lab)+
                  ggtitle("PC2 vs PC3")+
                  xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                  ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
                  theme_bw() +
                  theme(axis.title=element_text(size=10),
                        legend.title=element_blank(),
                        plot.title = element_text(hjust = 0.5))

                if (n.sample < 9 ) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                }
                ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC2_vs_PC3.", switch (output_images_file_format,
                                                                                                                    "pdf" = "pdf",
                                                                                                                    "png" = "png",
                                                                                                                    "tiff" = "tiff",
                                                                                                                    "jpeg" = "jpeg",
                                                                                                                    "svg" = "svg",
                                                                                                                    "bmp" = "bmp"
                )), width=10, height=8, units = "in")

                if (n.sample >= 9) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                  ggsave(p, filename=paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_PC2_vs_PC3_SampleName_label.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                  )), width=10, height=8, units = "in")
                }
              }

              ### 2.6 绘制3D PCA图
              sn <- NROW(TMM)
              color_list <- colors()
              sacol <- sample(color_list, sn)
              sl <- sn/10.0
              spos <- 1:sn
              #### 保存图片格式
              if (output_images_file_format == "pdf") { #### PDF ####
                ##### scatter3D #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
                dev.off()

                ##### text3D #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.pdf"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                pdf(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "png") { #### PNG ####
                ##### scatter3D #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.png"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                png(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "tiff") { #### TIFF ####
                ##### scatter3D #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.tiff"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                tiff(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "jpeg") { #### JPEG ####
                ##### scatter3D #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.jpeg"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                jpeg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "svg") { #### SVG ####
                ##### scatter3D #####
                svg(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                svg(paste0(path, "/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.svg"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                svg(path, paste0("/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "bmp") {#### BMP ####
                ##### scatter3D #####
                bmp(paste0(path, "/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Text.bmp"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                bmp(paste0(path,"/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else {
                print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
              }

              ### 2.7 绘制PCA动图
              p <- plot_ly(coord,
                           text = ~rownames(coord),
                           x = ~Dim.1,
                           y = ~Dim.2,
                           z = ~Dim.3,
                           alpha = 0.8) %>%
                add_markers(color = ~group) %>%
                layout(scene = list(
                  xaxis = list(title = dim1_lab),
                  yaxis = list(title = dim2_lab),
                  zaxis = list(title = dim3_lab)
                )
                )

              widget_file_size <- function(p) {
                d <- tempdir()
                withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
                f <- file.path(d, "index.html")
                mb <- round(file.info(f)$size / 1e6, 3)
                message("The 3D PCA plotly HTML File is: ", mb," MB")
              }

              widget_file_size(p)
              saveWidget(p, paste0(path, "/1.QualityControl/1.3Samples_PCA/", pl, "_PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

            }

          } else { # 不需要Olink分析文件夹系统
            # QC_Warning: 1.Warning/Pass; 2.WARN/PASS
            if (qc_warning == TRUE) { # qc_warning = TRUE

              ### 2.5 绘制2D PCA图
              coord <- as.data.frame(res.pca$ind$coord[,1:idx])
              coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
              coord$qc_warning <- sgi[match(rownames(coord), sgi$Sample),3]
              n.sample <- nrow(coord)
              n.group <- length(levels(as.factor(coord$group)))
              #### Save PCA Data ####
              rio::export(coord, file = paste0(path,"/Olink_", pl, "_PCA_Data.xlsx"), overwrite = T, zoom = 120, rowNames = T)
              #### PC1 vs PC2 ####
              p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
                #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
                scale_colour_manual(values = rainbow(n.group))+
                guides(shape="none") +
                xlab(dim1_lab)+
                ylab(dim2_lab)+
                ggtitle("PC1 vs PC2")+
                xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
                ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                theme_bw() +
                theme(axis.title=element_text(size=10),
                      legend.title=element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(size = guide_legend(),
                       shape = guide_legend()) # colour = guide_colorbar()

              if (n.sample < 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              }
              ggsave(p, filename=paste0(pl, "_PCA_PC1_vs_PC2_QC_warning.", switch (output_images_file_format,
                                                                                   "pdf" = "pdf",
                                                                                   "png" = "png",
                                                                                   "tiff" = "tiff",
                                                                                   "jpeg" = "jpeg",
                                                                                   "svg" = "svg",
                                                                                   "bmp" = "bmp"
              )), width=10, height=8, units = "in")
              if (n.sample >= 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                ggsave(p, filename=paste0(pl, "_PCA_PC1_vs_PC2_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                                      "pdf" = "pdf",
                                                                                                      "png" = "png",
                                                                                                      "tiff" = "tiff",
                                                                                                      "jpeg" = "jpeg",
                                                                                                      "svg" = "svg",
                                                                                                      "bmp" = "bmp"
                )), width=10, height=8, units = "in")
              }
              #### PC2 vs PC3 ####
              if ( idx == 3 ) {
                p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
                  #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                  geom_point(size = 3, aes(colour = factor(coord$group), shape = qc_warning))+
                  scale_colour_manual(values = rainbow(n.group))+
                  guides(shape="none") +
                  xlab(dim2_lab)+
                  ylab(dim3_lab)+
                  ggtitle("PC2 vs PC3")+
                  xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                  ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
                  theme_bw() +
                  theme(axis.title=element_text(size=10),
                        legend.title=element_blank(),
                        plot.title = element_text(hjust = 0.5)) +
                  guides(size = guide_legend(),
                         shape = guide_legend()) # colour = guide_colorbar()

                if (n.sample < 9 ) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                }
                ggsave(p, filename=paste0(pl, "_PCA_PC2_vs_PC3_QC_warning.", switch (output_images_file_format,
                                                                                     "pdf" = "pdf",
                                                                                     "png" = "png",
                                                                                     "tiff" = "tiff",
                                                                                     "jpeg" = "jpeg",
                                                                                     "svg" = "svg",
                                                                                     "bmp" = "bmp"
                )), width=10, height=8, units = "in")

                if (n.sample >= 9) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                  ggsave(p, filename=paste0(pl, "_PCA_PC2_vs_PC3_QC_warning_SampleName_label.", switch (output_images_file_format,
                                                                                                        "pdf" = "pdf",
                                                                                                        "png" = "png",
                                                                                                        "tiff" = "tiff",
                                                                                                        "jpeg" = "jpeg",
                                                                                                        "svg" = "svg",
                                                                                                        "bmp" = "bmp"
                  )), width=10, height=8, units = "in")
                }
              }

              ### 2.6 绘制3D PCA图
              sn <- NROW(TMM)
              color_list <- colors()
              sacol <- sample(color_list, sn)
              sl <- sn/10.0
              spos <- 1:sn
              #### 保存图片格式
              if (output_images_file_format == "pdf") { #### PDF ####
                ##### scatter3D #####
                pdf(paste0(pl, "_PCA_3D_Point.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
                dev.off()

                ##### text3D #####
                pdf(paste0(pl, "_PCA_3D_Text.pdf"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                pdf(paste0(pl, "_PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "png") { #### PNG ####
                ##### scatter3D #####
                png(paste0(pl, "_PCA_3D_Point.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                png(paste0(pl, "_PCA_3D_Text.png"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                png(paste0(pl, "_PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "tiff") { #### TIFF ####
                ##### scatter3D #####
                tiff(paste0(pl, "_PCA_3D_Point.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                tiff(paste0(pl, "_PCA_3D_Text.tiff"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                tiff(paste0(pl, "_PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "jpeg") { #### JPEG ####
                ##### scatter3D #####
                jpeg(paste0(pl, "_PCA_3D_Point.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                jpeg(paste0(pl, "_PCA_3D_Text.jpeg"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                jpeg(paste0(pl, "_PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "svg") { #### SVG ####
                ##### scatter3D #####
                svg(paste0(pl, "_PCA_3D_Point.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                svg(paste0(pl, "_PCA_3D_Text.svg"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                svg(paste0(pl, "_PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "bmp") {#### BMP ####
                ##### scatter3D #####
                bmp(paste0(pl, "_PCA_3D_Point.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                bmp(paste0(pl, "_PCA_3D_Text.bmp"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                bmp(paste0(pl, "_PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else {
                print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
              }

              ### 2.7 绘制PCA动图
              p <- plot_ly(coord,
                           text = ~rownames(coord),
                           x = ~Dim.1,
                           y = ~Dim.2,
                           z = ~Dim.3,
                           alpha = 0.8) %>%
                add_markers(color = ~group) %>%
                layout(scene = list(
                  xaxis = list(title = dim1_lab),
                  yaxis = list(title = dim2_lab),
                  zaxis = list(title = dim3_lab)
                )
                )

              widget_file_size <- function(p) {
                d <- tempdir()
                withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
                f <- file.path(d, "index.html")
                mb <- round(file.info(f)$size / 1e6, 3)
                message("The 3D PCA plotly HTML File is: ", mb," MB")
              }

              widget_file_size(p)
              saveWidget(p, paste0(pl, "_PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

            } else { # qc_warning = FALSE
              print("qc_warning == FALSE，不需要为PCA图添加QC_Warning注释！")

              ### 2.5 绘制2D PCA图
              coord <- as.data.frame(res.pca$ind$coord[,1:idx])
              coord$group <- sgi[match(rownames(coord), sgi$Sample),2]
              n.sample <- nrow(coord)
              n.group <- length(levels(as.factor(coord$group)))
              #### PC1 vs PC2 ####
              p <- ggplot(coord, aes(x = Dim.1, y = Dim.2)) +
                #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                geom_point(size = 3, aes(colour = factor(coord$group)))+
                scale_colour_manual(values = rainbow(n.group))+
                guides(shape="none") +
                xlab(dim1_lab)+
                ylab(dim2_lab)+
                ggtitle("PC1 vs PC2")+
                xlim(c(min(coord$Dim.1)*1.2, max(coord$Dim.1)*1.2))+
                ylim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                theme_bw() +
                theme(axis.title=element_text(size=10),
                      legend.title=element_blank(),
                      plot.title = element_text(hjust = 0.5))

              if (n.sample < 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
              }
              ggsave(p, filename=paste0(pl, "_PCA_PC1_vs_PC2.", switch (output_images_file_format,
                                                                        "pdf" = "pdf",
                                                                        "png" = "png",
                                                                        "tiff" = "tiff",
                                                                        "jpeg" = "jpeg",
                                                                        "svg" = "svg",
                                                                        "bmp" = "bmp"
              )), width=10, height=8, units = "in")
              if (n.sample >= 9) {
                p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                ggsave(p, filename=paste0(pl, "_PCA_PC1_vs_PC2_SampleName_label.", switch (output_images_file_format,
                                                                                           "pdf" = "pdf",
                                                                                           "png" = "png",
                                                                                           "tiff" = "tiff",
                                                                                           "jpeg" = "jpeg",
                                                                                           "svg" = "svg",
                                                                                           "bmp" = "bmp"
                )), width=10, height=8, units = "in")
              }
              #### PC2 vs PC3 ####
              if ( idx == 3 ) {
                p <- ggplot(coord, aes(x=Dim.2, y=Dim.3)) +
                  #geom_point(size=3, aes(colour=factor(coord$group), shape=group))+
                  geom_point(size = 3, aes(colour = factor(coord$group)))+
                  scale_colour_manual(values = rainbow(n.group))+
                  guides(shape="none") +
                  xlab(dim2_lab)+
                  ylab(dim3_lab)+
                  ggtitle("PC2 vs PC3")+
                  xlim(c(min(coord$Dim.2)*1.2, max(coord$Dim.2)*1.2))+
                  ylim(c(min(coord$Dim.3)*1.2, max(coord$Dim.3)*1.2))+
                  theme_bw() +
                  theme(axis.title=element_text(size=10),
                        legend.title=element_blank(),
                        plot.title = element_text(hjust = 0.5))

                if (n.sample < 9 ) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                }
                ggsave(p, filename=paste0(pl, "_PCA_PC2_vs_PC3.", switch (output_images_file_format,
                                                                          "pdf" = "pdf",
                                                                          "png" = "png",
                                                                          "tiff" = "tiff",
                                                                          "jpeg" = "jpeg",
                                                                          "svg" = "svg",
                                                                          "bmp" = "bmp"
                )), width=10, height=8, units = "in")

                if (n.sample >= 9) {
                  p <- p + geom_text_repel(size=3, aes(label=rownames(coord)), hjust=0, vjust=0)
                  ggsave(p, filename=paste0(pl, "_PCA_PC2_vs_PC3_SampleName_label.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                  )), width=10, height=8, units = "in")
                }
              }

              ### 2.6 绘制3D PCA图
              sn <- NROW(TMM)
              color_list <- colors()
              sacol <- sample(color_list, sn)
              sl <- sn/10.0
              spos <- 1:sn
              #### 保存图片格式
              if (output_images_file_format == "pdf") { #### PDF ####
                ##### scatter3D #####
                pdf(paste0(pl, "_PCA_3D_Point.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                # colkey = list(at = spos, side = 4, addlines = FALSE, length = sl, width = 0.6, labels = rownames(coord), cex.clab = 1.5), clab = "Sample", colvar = as.integer(rownames(coord))
                dev.off()

                ##### text3D #####
                pdf(paste0(pl, "_PCA_3D_Text.pdf"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                pdf(paste0(pl, "_PCA_3D_Point_And_Text.pdf"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "png") { #### PNG ####
                ##### scatter3D #####
                png(paste0(pl, "_PCA_3D_Point.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                png(paste0(pl, "_PCA_3D_Text.png"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                png(paste0(pl, "_PCA_3D_Point_And_Text.png"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "tiff") { #### TIFF ####
                ##### scatter3D #####
                tiff(paste0(pl, "_PCA_3D_Point.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                tiff(paste0(pl, "_PCA_3D_Text.tiff"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                tiff(paste0(pl, "_PCA_3D_Point_And_Text.tiff"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else if (output_images_file_format == "jpeg") { #### JPEG ####
                ##### scatter3D #####
                jpeg(paste0(pl, "_PCA_3D_Point.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                jpeg(paste0(pl, "_PCA_3D_Text.jpeg"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                jpeg(paste0(pl, "_PCA_3D_Point_And_Text.jpeg"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "svg") { #### SVG ####
                ##### scatter3D #####
                svg(paste0(pl, "_PCA_3D_Point.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                svg(paste0(pl, "_PCA_3D_Text.svg"), width = 12, height = 12)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                svg(paste0(pl, "_PCA_3D_Point_And_Text.svg"), width = 12, height = 12)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #

              } else if (output_images_file_format == "bmp") {#### BMP ####
                ##### scatter3D #####
                bmp(paste0(pl, "_PCA_3D_Point.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### text3D #####
                bmp(paste0(pl, "_PCA_3D_Text.bmp"), width = 1200, height = 1200)
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                dev.off()

                ##### Point and Text #####
                bmp(paste0(pl, "_PCA_3D_Point_And_Text.bmp"), width = 1200, height = 1200)
                scatter3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 5, ticktype = "detailed", plot = T, adj = 0.5, font = 2, colkey = F)
                #
                text3D(x = coord$Dim.1, y = coord$Dim.2, z = coord$Dim.3, xlab = dim1_lab, ylab = dim2_lab, zlab = dim3_lab, labels = rownames(coord), col = sacol, phi = 5, theta = 32, bty = 'b2', pch = 20, cex = 1, ticktype = "detailed", plot = T, adj = 0.5, font = 2, add = T)
                dev.off()
                #
              } else {
                print("注意：output_images_file_format的选项不在设定范围内，请检查确认后重试！")
              }

              ### 2.7 绘制PCA动图
              p <- plot_ly(coord,
                           text = ~rownames(coord),
                           x = ~Dim.1,
                           y = ~Dim.2,
                           z = ~Dim.3,
                           alpha = 0.8) %>%
                add_markers(color = ~group) %>%
                layout(scene = list(
                  xaxis = list(title = dim1_lab),
                  yaxis = list(title = dim2_lab),
                  zaxis = list(title = dim3_lab)
                )
                )

              widget_file_size <- function(p) {
                d <- tempdir()
                withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
                f <- file.path(d, "index.html")
                mb <- round(file.info(f)$size / 1e6, 3)
                message("The 3D PCA plotly HTML File is: ", mb," MB")
              }

              widget_file_size(p)
              saveWidget(p, paste0(pl, "_PCA_3D_Plotly.html"), selfcontained = F, libdir = "lib")

            }

          }
          #
          #
          pn = pn + 1 # 计数器计数
        }
        ##----------------------------------------------------------------------
        ########################################################################

      } else { # panel_pca == F
        print("注意：Olink_PCA函数的参数panel_pca = F，所以不进行Panel的PCA分析！")
      }

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two
      print("暂时不支持！")

    } else {
      print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
    }

  } else if (npx_data_type == "widetype") { # widetype
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two

    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }

  ## 脚本运行结束
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(cyan(bold("--------------------------------[2.Olink项目OAAT2-QC-PCA: Olink_PCA()函数运行完成！]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))


}

# ## 执行PCA分析绘图函数
# Olink_PCA(path = path,
#           npx_data_type = npx_data_type,
#           npx_data_files = npx_data_files,
#           target_npx_data = NULL,
#           sample_grouping_info = sample_grouping_info,
#           folder_system = folder_system,
#           output_images_file_format = output_images_file_format,
#           qc_warning = qc_warning,
#           panel_pca = panel_pca,
#           firstcol = "UniProt",
#           firstrow = "SampleID",
#           cellvalue = "NPX",
#           sample_control_name = sample_control_name,
#           filtered_samples = filtered_samples)
