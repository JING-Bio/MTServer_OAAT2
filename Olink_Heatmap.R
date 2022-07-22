# Olink Heatmap Analysis and Plotting Function Script

Olink_Heatmap <- function(path = path,
                          npx_data_type = npx_data_type,
                          npx_data_files = npx_data_files,
                          # sample_grouping_info = sample_grouping_info,
                          folder_system = folder_system,
                          output_images_file_format = output_images_file_format,
                          qc_warning = qc_warning,
                          # firstcol = "UniProt",
                          # firstrow = "SampleID",
                          # cellvalue = "NPX",
                          target_npx_data = NULL,
                          heatmap_color_style = "light", # dark/light
                          panel_heatmap = F, # Logical parameters, Whether Need Heatmap for every Panel
                          olink_project_type = "T96"# Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
                          # sample_control_name = sample_control_name,
                          # filtered_samples = filtered_samples
                          ) {
  ## 0.提前导入输出信息格式化设置R包
  ### BiocManager ###
  if (!requireNamespace("BiocManager", quietly = T)) {
    install.packages("BiocManager")
  }
  library(BiocManager)

  ### devtools ###
  if (!requireNamespace("devtools", quietly = T)) {
    install.packages("devtools")
  }
  library(devtools)
  ### crayon ###
  if (!requireNamespace("crayon", quietly = T)) {
    BiocManager::install("crayon")
  }
  library(crayon) #

  ## 1.脚本运行开始
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(red(bold("--------------------------------[1.开始执行Olink项目OAAT2-QC-Heatmap: Olink_Heatmap()函数]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

  ## 2.安装和导入R包
  analibpackages <- c("rio", "visdat", "htmlwidgets",
                      "pheatmap", "ggplot2", "corrplot",
                      "magrittr", "reshape2", "FactoMineR",
                      "openxlsx", "dplyr", "data.table",
                      "tibble", "parallel", "tidyr",
                      "purrr", "tidyverse", "stringr",
                      "ggsci", "ggpubr", "plotrix",
                      "igraph", "cowplot", "showtext",
                      "gridExtra", "grid", "scales",
                      "ggrepel", "plot3D", "readxl",
                      "pcaMethods", "plotly") # "crayon"

  for (a in analibpackages) {
    if (!requireNamespace(a, quietly = T)) {
      BiocManager::install(a)
    } else {
      # print(paste0("[软件包: '", a, "'已经安装完毕！]"))
      # print(packageVersion(a))
    }
  }

  sapply(analibpackages, library, character.only = T)

  # OlinkAnalyze #
  if (!requireNamespace("OlinkAnalyze", quietly = T)) {
    devtools::install_github(repo = 'Olink-Proteomics/OlinkRPackage/OlinkAnalyze',
                             build_vignettes = TRUE)
    library("OlinkAnalyze")
  } else {
    # print("[软件包: OlinkAnalyze 已经安装完毕！]")
    library("OlinkAnalyze")
  }

  ## Sample-Panel-Protein_Heatmap ##
  if (npx_data_type == "rawtype") { # rawtype

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1 && is.null(target_npx_data)) { # one and need to read npx_data_files
      npx_data <- read_NPX(paste0(path, "/", npx_data_files))
      olink_npx_panel_data <- npx_data
      samples_name_list <- npx_data$SampleID %>% unique() # 获取样本名称列表
      samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
      colnames(samples_name_df) <- "SampleID"
      uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
      if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
        uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
      } else {
        uniprotid_noNA_data <- uniprotid_no_NA_data
      }
      #
      protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
      protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
      colnames(protein_uniprotid_df) <- "UniProt"
      npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
      sn <- 0
      for (s in samples_name_list) {
        # print(s)
        sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
        # colnames(sample_npx_heatmap_data)[2] <- s
        sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
        npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
        sn = sn + 1
        colnames(npx_heatmap_data)[3*sn + 2] <- s
      }
      # QC_Warning：1.PASS/WARN；2.Pass/Warning
      if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 1.PASS/WARN
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "Warning"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "Pass"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 2.Pass/Warning
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "WARN"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "PASS"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else {
        print("QC_Warning的数据内容【1.PASS/WARN；2.Pass/Warning】不在预设范围内，请检查！")
      }

      #
      # print(dim(npx_heatmap_data))
      npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
      # print(head(npx_heatmap_data))
      rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
      npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]

      if (folder_system == T) {
        ## export heatmap data
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      } else { # folder_system = F
        # export heatmap data to the path
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 5,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      }

      ## Every Panel Heatmap ## ------------------------------------------------
      if (panel_heatmap == TRUE) { # 为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）
        print("panel_heatmap == TRUE，为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）")
        ### data processing
        olink_panel_list <- npx_data$Panel %>% unique()
        olink_npx_data <- npx_data
        # olink_npx_data <- npx_data %>%
        #   filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) # 不对对照样本进行筛选
        for (p in olink_panel_list) { # 循环各个panel绘制
          cat(bold(red(green("当前循环处理的Panel为：", red(p),"\n"))))
          #
          olink_npx_panel_data <- olink_npx_data[which(olink_npx_data$Panel == p),]
          samples_name_list <- olink_npx_data$SampleID %>% unique() # 获取样本名称列表
          samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
          colnames(samples_name_df) <- "SampleID"
          uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
          if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
            uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
          } else {
            uniprotid_noNA_data <- uniprotid_no_NA_data
          }

          protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
          protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
          colnames(protein_uniprotid_df) <- "UniProt"
          #
          npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
          sn <- 0
          for (s in samples_name_list) {
            # print(s)
            sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
            sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
            npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
            sn = sn + 1
            colnames(npx_heatmap_data)[3*sn + 2] <- s
          }
          #
          if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # QC_Warning: Warning/Pass
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "Warning"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "Pass"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "WARN"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "PASS"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else {
            print("注意：QC_Warning的数据选项不在预设范围内！")
          }

          #
          # print(dim(npx_heatmap_data))
          npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
          # print(head(npx_heatmap_data))
          rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
          npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]
          #
          if (folder_system == T) { # 有Olink分析结果文件系统，分析结果保存至对应的文件夹中
            rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap as images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            }
          } else { # 没有Olink分析结果文件系统，分析结果保存至：path

            rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            }

          }

        }

      } else { # 为假则不绘制每个Panel的质控热图
        print("panel_heatmap == FALSE，为假则不绘制每个Panel的质控热图！")
      }

    } else if ( length(strsplit(npx_data_files, split = ",")[[1]]) == 1 && !is.null(target_npx_data)) { # one and use environment of npx_data
      olink_npx_panel_data <- npx_data
      samples_name_list <- npx_data$SampleID %>% unique()
      samples_name_list <- npx_data$SampleID %>% unique() # 获取样本名称列表
      samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
      colnames(samples_name_df) <- "SampleID"
      uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
      if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
        uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
      } else {
        uniprotid_noNA_data <- uniprotid_no_NA_data
      }
      #
      protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
      protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
      colnames(protein_uniprotid_df) <- "UniProt"
      npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
      sn <- 0
      for (s in samples_name_list) {
        # print(s)
        sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
        # colnames(sample_npx_heatmap_data)[2] <- s
        sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
        npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
        sn = sn + 1
        colnames(npx_heatmap_data)[3*sn + 2] <- s
      }
      # QC_Warning：1.PASS/WARN；2.Pass/Warning
      if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 1.PASS/WARN
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "Warning"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "Pass"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 2.Pass/Warning
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "WARN"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "PASS"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else {
        print("QC_Warning的数据内容【1.PASS/WARN；2.Pass/Warning】不在预设范围内，请检查！")
      }

      #
      # print(dim(npx_heatmap_data))
      npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
      # print(head(npx_heatmap_data))
      rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
      npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]

      if (folder_system == T) {
        ## export heatmap data
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      } else { # folder_system = F
        # export heatmap data to the path
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 5,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      }

      ## Every Panel Heatmap ## ------------------------------------------------
      if (panel_heatmap == TRUE) { # 为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）
        print("panel_heatmap == TRUE，为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）")
        ### data processing
        olink_panel_list <- npx_data$Panel %>% unique()
        olink_npx_data <- npx_data
        # olink_npx_data <- npx_data %>%
        #   filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) # 不对对照样本进行筛选
        for (p in olink_panel_list) { # 循环各个panel绘制
          cat(bold(red(green("当前循环处理的Panel为：", red(p),"\n"))))
          #
          olink_npx_panel_data <- olink_npx_data[which(olink_npx_data$Panel == p),]
          samples_name_list <- olink_npx_data$SampleID %>% unique() # 获取样本名称列表
          samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
          colnames(samples_name_df) <- "SampleID"
          uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
          if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
            uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
          } else {
            uniprotid_noNA_data <- uniprotid_no_NA_data
          }

          protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
          protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
          colnames(protein_uniprotid_df) <- "UniProt"
          #
          npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
          sn <- 0
          for (s in samples_name_list) {
            # print(s)
            sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
            sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
            npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
            sn = sn + 1
            colnames(npx_heatmap_data)[3*sn + 2] <- s
          }
          #
          if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # QC_Warning: Warning/Pass
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "Warning"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "Pass"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "WARN"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "PASS"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else {
            print("注意：QC_Warning的数据选项不在预设范围内！")
          }

          #
          # print(dim(npx_heatmap_data))
          npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
          # print(head(npx_heatmap_data))
          rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
          npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]
          #
          if (folder_system == T) { # 有Olink分析结果文件系统，分析结果保存至对应的文件夹中
            rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap as images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            }
          } else { # 没有Olink分析结果文件系统，分析结果保存至：path

            rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            }

          }

        }

      } else { # 为假则不绘制每个Panel的质控热图
        print("panel_heatmap == FALSE，为假则不绘制每个Panel的质控热图！")
      }

    } else if ( length(strsplit(npx_data_files, split = ",")[[1]]) == 2 && is.null(target_npx_data)) { # two and need to read npx_data_files
      print("Sample-Panel-Protein_Heatmap-rawtype-two and need to read npx_data_files")
    } else if ( length(strsplit(npx_data_files, split = ",")[[1]]) == 2 && !is.null(target_npx_data)) { # two and use environment of npx_data
      print("Sample-Panel-Protein_Heatmap-rawtype-two and use environment of npx_data")
    } else { # more than 2
      print("多于2个NPX Raw Data文件！暂时不能处理！")
    }

  } else if (npx_data_type == "longtype") { # longtype

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1 && is.null(target_npx_data)) { # one and need to read npx_data_files

      ## All Panel ## ----------------------------------------------------------
      ### plotting data process
      npx_data <- read.xlsx(paste0(path, "/", npx_data_files))
      olink_npx_panel_data <- npx_data
      samples_name_list <- npx_data$SampleID %>% unique() # 获取样本名称列表
      samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
      colnames(samples_name_df) <- "SampleID"
      uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
      if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
        uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
      } else {
        uniprotid_noNA_data <- uniprotid_no_NA_data
      }
      #
      protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
      protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
      colnames(protein_uniprotid_df) <- "UniProt"
      npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
      sn <- 0
      for (s in samples_name_list) {
        # print(s)
        sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
        # colnames(sample_npx_heatmap_data)[2] <- s
        sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
        npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
        sn = sn + 1
        colnames(npx_heatmap_data)[3*sn + 2] <- s
      }
      # QC_Warning：1.PASS/WARN；2.Pass/Warning
      if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 1.PASS/WARN
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "Warning"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "Pass"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 2.Pass/Warning
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "WARN"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "PASS"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else {
        print("QC_Warning的数据内容【1.PASS/WARN；2.Pass/Warning】不在预设范围内，请检查！")
      }

      #
      # print(dim(npx_heatmap_data))
      npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
      # print(head(npx_heatmap_data))
      rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
      npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]

      if (folder_system == T) {
        ## export heatmap data
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      } else { # folder_system = F
        # export heatmap data to the path
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 5,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      }

      ## Every Panel Heatmap ## ------------------------------------------------
      if (panel_heatmap == TRUE) { # 为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）
        print("panel_heatmap == TRUE，为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）")
        ### data processing
        olink_panel_list <- npx_data$Panel %>% unique()
        olink_npx_data <- npx_data
        # olink_npx_data <- npx_data %>%
        #   filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) # 不对对照样本进行筛选
        for (p in olink_panel_list) { # 循环各个panel绘制
          cat(bold(red(green("当前循环处理的Panel为：", red(p),"\n"))))
          #
          olink_npx_panel_data <- olink_npx_data[which(olink_npx_data$Panel == p),]
          samples_name_list <- olink_npx_data$SampleID %>% unique() # 获取样本名称列表
          samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
          colnames(samples_name_df) <- "SampleID"
          uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
          if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
            uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
          } else {
            uniprotid_noNA_data <- uniprotid_no_NA_data
          }

          protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
          protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
          colnames(protein_uniprotid_df) <- "UniProt"
          #
          npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
          sn <- 0
          for (s in samples_name_list) {
            # print(s)
            sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
            sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
            npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
            sn = sn + 1
            colnames(npx_heatmap_data)[3*sn + 2] <- s
          }
          #
          if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # QC_Warning: Warning/Pass
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "Warning"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "Pass"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "WARN"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "PASS"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else {
            print("注意：QC_Warning的数据选项不在预设范围内！")
          }

          #
          # print(dim(npx_heatmap_data))
          npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
          # print(head(npx_heatmap_data))
          rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
          npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]
          #
          if (folder_system == T) { # 有Olink分析结果文件系统，分析结果保存至对应的文件夹中
            rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap as images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            }
          } else { # 没有Olink分析结果文件系统，分析结果保存至：path

            rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            }

          }

        }

      } else { # 为假则不绘制每个Panel的质控热图
        print("panel_heatmap == FALSE，为假则不绘制每个Panel的质控热图！")
      }
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1 && !is.null(target_npx_data)) { # one and use environment npx_data
      ## All Panel ## ----------------------------------------------------------
      ### plotting data process
      # npx_data <- read.xlsx(paste0(path, "/", npx_data_files))
      olink_npx_panel_data <- npx_data
      samples_name_list <- npx_data$SampleID %>% unique() # 获取样本名称列表
      samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
      colnames(samples_name_df) <- "SampleID"
      uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
      if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
        uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
      } else {
        uniprotid_noNA_data <- uniprotid_no_NA_data
      }
      #
      protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
      protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
      colnames(protein_uniprotid_df) <- "UniProt"
      npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
      sn <- 0
      for (s in samples_name_list) {
        # print(s)
        sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
        # colnames(sample_npx_heatmap_data)[2] <- s
        sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
        npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
        sn = sn + 1
        colnames(npx_heatmap_data)[3*sn + 2] <- s
      }
      # QC_Warning：1.PASS/WARN；2.Pass/Warning
      if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 1.PASS/WARN
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "Warning"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "Pass"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) { # 2.Pass/Warning
        sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
        warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
        colnames(warning_samples_df) <- "SampleID"
        warning_samples_df$QC_Warning <- "WARN"
        annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
        annotation_col_df[is.na(annotation_col_df)] <- "PASS"
        rownames(annotation_col_df) <- annotation_col_df$SampleID
        annotation_col <- subset(annotation_col_df, select = -SampleID)
      } else {
        print("QC_Warning的数据内容【1.PASS/WARN；2.Pass/Warning】不在预设范围内，请检查！")
      }

      #
      # print(dim(npx_heatmap_data))
      npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
      # print(head(npx_heatmap_data))
      rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
      npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]

      if (folder_system == T) {
        ## export heatmap data
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                       "pdf" = "pdf",
                                                                                                                                       "png" = "png",
                                                                                                                                       "tiff" = "tiff",
                                                                                                                                       "jpeg" = "jpeg",
                                                                                                                                       "svg" = "svg",
                                                                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      } else { # folder_system = F
        # export heatmap data to the path
        rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_All_Panel_Protein_Heatmap.xlsx"), zoom = 120, overwrite = T, rowNames = T)

        ## plotting
        if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
        } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
          ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
        } else {
          print("QC_Warning的数据内容不在预设范围内，请检查！")
        }
        #
        if (heatmap_color_style == "light") { # heatmap_color_style = "light"
          #### Save plotting heatmap image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else if (heatmap_color_style == "dark") { # heatmap_color_style = "dark"
          #### Save plotting heatmap as image file format
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 1,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        } else { # 默认配色绘制
          #### PDF plotting heatmap
          pheatmap(npx_heatmap_data_plot,
                   cluster_rows = T,
                   display_numbers = FALSE,
                   number_format = "%.3f",
                   show_colnames = T,
                   show_rownames = F,
                   fontsize = 18,
                   color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellwidth = 12,
                   cellheight = 5,
                   scale = "column",
                   fontsize_col = 10,
                   # fontsize_row = 8,
                   angle_col = 90,
                   filename = paste0(path,"/Olink_Panel_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                       "pdf" = "pdf",
                                                                                       "png" = "png",
                                                                                       "tiff" = "tiff",
                                                                                       "jpeg" = "jpeg",
                                                                                       "svg" = "svg",
                                                                                       "bmp" = "bmp"
                   )))
          cat(bgWhite(red("Olink_Panel_UniProt_NPX_heatmap Plotting Finished!"),"\n"))
        }

      }

      ## Every Panel Heatmap ## ------------------------------------------------
      if (panel_heatmap == TRUE) { # 为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）
        print("panel_heatmap == TRUE，为真则绘制每个Panel的质控热图（所谓质控热图，就是不进行数据筛选，直接用所有的样本数据，包括对照样本的NPX数据进行热图的绘制）")
        ### data processing
        olink_panel_list <- npx_data$Panel %>% unique()
        olink_npx_data <- npx_data
        # olink_npx_data <- npx_data %>%
        #   filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) # 不对对照样本进行筛选
        for (p in olink_panel_list) { # 循环各个panel绘制
          cat(bold(red(green("当前循环处理的Panel为：", red(p),"\n"))))
          #
          olink_npx_panel_data <- olink_npx_data[which(olink_npx_data$Panel == p),]
          samples_name_list <- olink_npx_data$SampleID %>% unique() # 获取样本名称列表
          samples_name_df <- as.data.frame(samples_name_list, stringsAsFactors = F)
          colnames(samples_name_df) <- "SampleID"
          uniprotid_no_NA_data <- olink_npx_panel_data[!is.na(olink_npx_panel_data$UniProt),] # 删除UniProt列含缺失值的行
          if ("NA" %in% uniprotid_no_NA_data$UniProt == T) {
            uniprotid_noNA_data <- uniprotid_no_NA_data[-which(uniprotid_no_NA_data$UniProt == "NA"),] # 删除UniProt列含NA的行
          } else {
            uniprotid_noNA_data <- uniprotid_no_NA_data
          }

          protein_uniprotid_list <- uniprotid_noNA_data$UniProt %>% unique()
          protein_uniprotid_df <- as.data.frame(protein_uniprotid_list, stringsAsFactors = F)
          colnames(protein_uniprotid_df) <- "UniProt"
          #
          npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == samples_name_list[1]), c("UniProt", "NPX", "QC_Warning")]
          sn <- 0
          for (s in samples_name_list) {
            # print(s)
            sample_npx_heatmap_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$SampleID == s), c("UniProt", "NPX", "QC_Warning")]
            sample_npx_heatmap_data <- merge.data.frame(protein_uniprotid_df, sample_npx_heatmap_data, by.x = "UniProt", by.y = "UniProt", all.x = T)
            npx_heatmap_data <- cbind(npx_heatmap_data, sample_npx_heatmap_data)
            sn = sn + 1
            colnames(npx_heatmap_data)[3*sn + 2] <- s
          }
          #
          if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) { # QC_Warning: Warning/Pass
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "Warning"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "Warning"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "Pass"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
            sample_npx_warning_data <- uniprotid_noNA_data[which(uniprotid_noNA_data$QC_Warning == "WARN"), c("SampleID", "UniProt", "Panel", "QC_Warning")]
            warning_samples_df <-  as.data.frame(sample_npx_warning_data$SampleID %>% unique(), stringsAsFactors = F)
            colnames(warning_samples_df) <- "SampleID"
            warning_samples_df$QC_Warning <- "WARN"
            annotation_col_df <- merge.data.frame(samples_name_df, warning_samples_df, by.x = "SampleID", by.y = "SampleID", all = T)
            annotation_col_df[is.na(annotation_col_df)] <- "PASS"
            rownames(annotation_col_df) <- annotation_col_df$SampleID
            annotation_col <- subset(annotation_col_df, select = -SampleID)
          } else {
            print("注意：QC_Warning的数据选项不在预设范围内！")
          }

          #
          # print(dim(npx_heatmap_data))
          npx_heatmap_data <- npx_heatmap_data[!duplicated(npx_heatmap_data$UniProt),]
          # print(head(npx_heatmap_data))
          rownames(npx_heatmap_data) <- npx_heatmap_data$UniProt
          npx_heatmap_data_plot <- npx_heatmap_data[,c(samples_name_list)]
          #
          if (folder_system == T) { # 有Olink分析结果文件系统，分析结果保存至对应的文件夹中
            rio::export(npx_heatmap_data_plot, file = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap as images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                                                                             "pdf" = "pdf",
                                                                                                                                             "png" = "png",
                                                                                                                                             "tiff" = "tiff",
                                                                                                                                             "jpeg" = "jpeg",
                                                                                                                                             "svg" = "svg",
                                                                                                                                             "bmp" = "bmp"
                       )))
              cat(bgWhite(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n")))
            }
          } else { # 没有Olink分析结果文件系统，分析结果保存至：path

            rio::export(npx_heatmap_data_plot, file = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.xlsx"), zoom = 120, rowNames = T, colNames = T, overwrite = T)
            #
            if ("Warning" %in% unique(uniprotid_noNA_data$QC_Warning) || "Pass" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("Warning" = "red","Pass" = "lightseagreen"))
            } else if ("WARN" %in% unique(uniprotid_noNA_data$QC_Warning) || "PASS" %in% unique(uniprotid_noNA_data$QC_Warning)) {
              ann_colors = list(QC_Warning = c("WARN" = "red","PASS" = "lightseagreen"))
            } else {
              print("QC_Warning的数据内容不在预设范围内，请检查！")
            }
            #
            if (heatmap_color_style == "light") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else if (heatmap_color_style == "dark") {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("blue","black","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            } else {
              #### Save plotting heatmap images format
              pheatmap(npx_heatmap_data_plot,
                       cluster_rows = T,
                       display_numbers = FALSE,
                       number_format = "%.3f",
                       show_colnames = T,
                       show_rownames = T,
                       annotation_colors = ann_colors,
                       # color = colorRampPalette(c("blue","white","red"), bias = 1)(1000),
                       border_color = NA,
                       annotation_col = annotation_col,
                       cellwidth = 12,
                       cellheight = 10,
                       scale = "row",
                       fontsize_col = 10,
                       fontsize_row = 8,
                       angle_col = 90,
                       fontsize = 18,
                       filename = paste0(path,"/Olink_", p, "_UniProt_NPX_heatmap.", switch (output_images_file_format,
                                                                                             "pdf" = "pdf",
                                                                                             "png" = "png",
                                                                                             "tiff" = "tiff",
                                                                                             "jpeg" = "jpeg",
                                                                                             "svg" = "svg",
                                                                                             "bmp" = "bmp"
                       )))
              cat((bgGreen(blue(paste0("Olink_",p,"_UniProt_NPX_heatmap Plotting Finished!"),"\n"))))
            }

          }

        }

      } else { # 为假则不绘制每个Panel的质控热图
        print("panel_heatmap == FALSE，为假则不绘制每个Panel的质控热图！")
      }
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2 && is.null(target_npx_data)) { # two and need to read npx_data_files
      print("Sample-Panel-Protein_Heatmap-[longtype]-two and need to read npx_data_files")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2 && !is.null(target_npx_data)) { # two and use environment npx_data
      print("Sample-Panel-Protein_Heatmap-[longtype]-two and use environment of npx_data")
    } else {
      print("多于2个NPX Long Data文件！暂时不能处理！")
    }

  } else if (npx_data_type == "widetype") {
    print("暂不支持widetype的热图绘制！")
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }

  ## 脚本运行结束
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(red(bold("--------------------------------[1.Olink项目OAAT2-QC-Heatmap: Olink_Heatmap()函数运行完成！]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

}

## 执行Olink_Heatmap
# Olink_Heatmap(path = "/Users/mingma_10000455/P/OLINK项目/医创云康",
#               npx_data_type = "rawtype",
#               npx_data_files = "医创云康88血浆_Raw_NPX_Data.xlsx",
#               # sample_grouping_info = sample_grouping_info,
#               folder_system = F,
#               output_images_file_format = "pdf",
#               # qc_warning = qc_warning,
#               # firstcol = "UniProt",
#               # firstrow = "SampleID",
#               # cellvalue = "NPX",
#               target_npx_data = NULL,
#               heatmap_color_style = "light", # dark/light
#               panel_heatmap = T, # Logical parameters, Whether Need Heatmap for every Panel
#               olink_project_type = "T96" # Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
#               # sample_control_name = sample_control_name,
#               # filtered_samples = filtered_samples
#               )
##
