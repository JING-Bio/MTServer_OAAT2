# Olink Corrplot Plotting Function Script

Olink_Corrplot <- function(path = path,
                           npx_data_type = npx_data_type,
                           npx_data_files = npx_data_files,
                           folder_system = folder_system,
                           output_images_file_format = output_images_file_format,
                           qc_warning = qc_warning,
                           panel_corrplot = panel_corrplot,
                           firstcol = "UniProt",
                           firstrow = "SampleID",
                           cellvalue = "NPX",
                           sample_control_name = sample_control_name,
                           filtered_samples = filtered_samples,
                           corrplot_method = corrplot_method,
                           corrplot_type = corrplot_type,
                           corrplot_style = corrplot_style) {
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
                      "htmlwidgets", "corrplot")
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

  ## 1.5 开始执行脚本
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(green(bold("--------------------------------[3.开始执行Olink项目OAAT2-QC-Corrplot: Olink_Corrplot()函数]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

  ## 2.根据函数参数获取Olink的NPX数据
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") {

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {
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
      npx_corrplot_file_index <- list.files(path = path, pattern = "ALL_Panel_NPX_WideData")
      npx_corrplot_df <- import(paste0(path, "/",npx_corrplot_file_index))
      # Long NPX Data
      npx_data <- import(paste0(path, "/", npx_data_files))
      ### 2.4 根据不同参数绘制相关性系数热
      if (folder_system == T) { # Logical parameters，需要Olink分析结果文件系统

        if (corrplot_style == T) { # Logical parameters，直接绘制设定好的相关性系数热图风格：lower: number + upper: circle, color: default，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type不生效
          rownames(npx_corrplot_df) <- npx_corrplot_df$UniProt
          npx_corrplot_df$UniProt <- NULL
          npx_corrplot_df_noNA <- npx_corrplot_df[complete.cases(npx_corrplot_df),]

          npx_corrplot_df_noNA_cor <- npx_corrplot_df_noNA %>% cor()
          rio::export(npx_corrplot_df_noNA_cor, file = paste0(path, "/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)

          #### Images File Format ####
          if (output_images_file_format == "pdf" ) { # PDF
            pdf(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "png" ) { # PNG
            png(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "tiff" ) { # TIFF
            tiff(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "jpeg" ) { # JPEG
            jpeg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "bmp" ) { # JPEG
            bmp(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "svg" ) { # JPEG
            svg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else {
            print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
          }

        } else { # 不需要设定好的绘图风格，要自定义，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type生效
          rownames(npx_corrplot_df) <- npx_corrplot_df$UniProt
          npx_corrplot_df$UniProt <- NULL
          npx_corrplot_df_noNA <- npx_corrplot_df[complete.cases(npx_corrplot_df),]

          npx_corrplot_df_noNA_cor <- npx_corrplot_df_noNA %>% cor()
          rio::export(npx_corrplot_df_noNA_cor, file = paste0(path, "/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
          #### Images File Format ####
          if (output_images_file_format == "pdf" ) { # PDF
            pdf(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            # corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()
          } else if (output_images_file_format == "png" ) { # PNG
            png(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "tiff" ) { # TIFF
            tiff(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "jpeg" ) { # JPEG
            jpeg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "bmp" ) { # BMP
            bmp(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "svg" ) { # SVG
            svg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else {
            print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
          }
        }
      } else { # 不需要Olink分析结果文件系统
        if (corrplot_style == T) { # Logical parameters，直接绘制设定好的相关性系数热图风格：lower: number + upper: circle, color: default，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type不生效
          rownames(npx_corrplot_df) <- npx_corrplot_df$UniProt
          npx_corrplot_df$UniProt <- NULL
          npx_corrplot_df_noNA <- npx_corrplot_df[complete.cases(npx_corrplot_df),]

          npx_corrplot_df_noNA_cor <- npx_corrplot_df_noNA %>% cor()
          rio::export(npx_corrplot_df_noNA_cor, file = paste0(path, "/Olink_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
          #### Images File Format ####
          if (output_images_file_format == "pdf" ) { # PDF
            pdf(paste0(path,"/Olink_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "png" ) { # PNG
            png(paste0(path,"/Olink_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "tiff" ) { # TIFF
            tiff(paste0(path,"/Olink_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "jpeg" ) { # JPEG
            jpeg(paste0(path,"/Olink_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "bmp" ) { # BMP
            bmp(paste0(path,"/Olink_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else if (output_images_file_format == "svg" ) { # SVG
            svg(paste0(path,"/Olink_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
            dev.off()

          } else {
            print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
          }
        } else { # 不需要设定好的绘图风格，要自定义，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type生效
          rownames(npx_corrplot_df) <- npx_corrplot_df$UniProt
          npx_corrplot_df$UniProt <- NULL
          npx_corrplot_df_noNA <- npx_corrplot_df[complete.cases(npx_corrplot_df),]

          npx_corrplot_df_noNA_cor <- npx_corrplot_df_noNA %>% cor()
          rio::export(npx_corrplot_df_noNA_cor, file = paste0(path, "/Olink_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
          #### Images File Format ####
          if (output_images_file_format == "pdf" ) { # PDF
            pdf(paste0(path,"/Olink_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "png" ) { # PNG
            png(paste0(path,"/Olink_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "tiff" ) { # TIFF
            tiff(paste0(path,"/Olink_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "jpeg" ) { # JPEG
            jpeg(paste0(path,"/Olink_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "bmp" ) { # BMP
            bmp(paste0(path,"/Olink_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else if (output_images_file_format == "svg" ) { # SVG
            svg(paste0(path,"/Olink_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
            M <- corrplot(npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
            dev.off()

          } else {
            print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
          }
        }
      }
      ##########################################################################
      ### Panel_Corrplot_Plotting ###
      ##########################################################################
      if (panel_corrplot == T) { # panel_corrplot == T
        pn <- 1 # 计数器初始化
        for (pl in unique(npx_data$Panel)) {
          cat(blue("[ ", pn, ".当前处理的Panel为："), bold(red(pl,"]\n")))
          #---------------------------------------------------------------------
          ### 2.3 获取Panel NPX宽数据
          panel_npx_corrplot_file_index <- list.files(path = path, pattern = paste0("*", pl, "_Panel_WideData*"))[1]
          panel_npx_corrplot_df <- import(paste0(path, "/", panel_npx_corrplot_file_index))
          ### 2.4 根据不同参数绘制相关性系数热
          if (folder_system == T) { # Logical parameters，需要Olink分析结果文件系统

            if (corrplot_style == T) { # Logical parameters，直接绘制设定好的相关性系数热图风格：lower: number + upper: circle, color: default，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type不生效
              rownames(panel_npx_corrplot_df) <- panel_npx_corrplot_df$UniProt
              panel_npx_corrplot_df$UniProt <- NULL
              panel_npx_corrplot_df_noNA <- panel_npx_corrplot_df[complete.cases(panel_npx_corrplot_df),]

              panel_npx_corrplot_df_noNA_cor <- panel_npx_corrplot_df_noNA %>% cor()
              rio::export(panel_npx_corrplot_df_noNA_cor, file = paste0(path, "/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)

              #### Images File Format ####
              if (output_images_file_format == "pdf" ) { # PDF
                pdf(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "png" ) { # PNG
                png(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "tiff" ) { # TIFF
                tiff(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "jpeg" ) { # JPEG
                jpeg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "bmp" ) { # JPEG
                bmp(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "svg" ) { # JPEG
                svg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else {
                print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
              }

            } else { # 不需要设定好的绘图风格，要自定义，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type生效
              rownames(panel_npx_corrplot_df) <- panel_npx_corrplot_df$UniProt
              panel_npx_corrplot_df$UniProt <- NULL
              panel_npx_corrplot_df_noNA <- panel_npx_corrplot_df[complete.cases(panel_npx_corrplot_df),]

              panel_npx_corrplot_df_noNA_cor <- panel_npx_corrplot_df_noNA %>% cor()
              rio::export(panel_npx_corrplot_df_noNA_cor, file = paste0(path, "/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
              #### Images File Format ####
              if (output_images_file_format == "pdf" ) { # PDF
                pdf(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                # corrplot(npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()
              } else if (output_images_file_format == "png" ) { # PNG
                png(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "tiff" ) { # TIFF
                tiff(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "jpeg" ) { # JPEG
                jpeg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "bmp" ) { # BMP
                bmp(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "svg" ) { # SVG
                svg(paste0(path,"/1.QualityControl/1.4Samples_Corrplot/Olink_", pl, "_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else {
                print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
              }
            }
          } else { # 不需要Olink分析结果文件系统
            if (corrplot_style == T) { # Logical parameters，直接绘制设定好的相关性系数热图风格：lower: number + upper: circle, color: default，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type不生效
              rownames(panel_npx_corrplot_df) <- panel_npx_corrplot_df$UniProt
              panel_npx_corrplot_df$UniProt <- NULL
              panel_npx_corrplot_df_noNA <- panel_npx_corrplot_df[complete.cases(panel_npx_corrplot_df),]

              panel_npx_corrplot_df_noNA_cor <- panel_npx_corrplot_df_noNA %>% cor()
              rio::export(panel_npx_corrplot_df_noNA_cor, file = paste0(path, "/Olink_", pl, "_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
              #### Images File Format ####
              if (output_images_file_format == "pdf" ) { # PDF
                pdf(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "png" ) { # PNG
                png(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "tiff" ) { # TIFF
                tiff(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "jpeg" ) { # JPEG
                jpeg(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "bmp" ) { # BMP
                bmp(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else if (output_images_file_format == "svg" ) { # SVG
                svg(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = "circle", type = "upper", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                corrplot(panel_npx_corrplot_df_noNA_cor, method = "number", type = "lower", tl.col = "black", tl.cex = 0.5, order = "hclust", diag = F, add = T, tl.pos = "n", cl.cex = 0.1, col = 'black', cl.pos = "n", number.cex = 0.5)
                dev.off()

              } else {
                print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
              }
            } else { # 不需要设定好的绘图风格，要自定义，此时Olink_Corrplot函数参数：corrplot_method和corrplot_type生效
              rownames(panel_npx_corrplot_df) <- panel_npx_corrplot_df$UniProt
              panel_npx_corrplot_df$UniProt <- NULL
              panel_npx_corrplot_df_noNA <- panel_npx_corrplot_df[complete.cases(panel_npx_corrplot_df),]

              panel_npx_corrplot_df_noNA_cor <- panel_npx_corrplot_df_noNA %>% cor()
              rio::export(panel_npx_corrplot_df_noNA_cor, file = paste0(path, "/Olink_", pl, "_NPX_Samples_Corrplot_Data.xlsx"), rowNames = T, zoom = 120, overwrite = T)
              #### Images File Format ####
              if (output_images_file_format == "pdf" ) { # PDF
                pdf(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.pdf"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "png" ) { # PNG
                png(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.png"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "tiff" ) { # TIFF
                tiff(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.tiff"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "jpeg" ) { # JPEG
                jpeg(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.jpeg"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "bmp" ) { # BMP
                bmp(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.bmp"), width = 1200, height = 1200)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else if (output_images_file_format == "svg" ) { # SVG
                svg(paste0(path,"/Olink_", pl, "_NPX_Samples_Corrplot.svg"), width = 12, height = 12)
                M <- corrplot(panel_npx_corrplot_df_noNA_cor, method = corrplot_method, type = corrplot_type, tl.col = "black", tl.cex = 0.5, order = "hclust", diag = T, tl.pos = "lt")
                dev.off()

              } else {
                print("注意：函数参数output_images_file_format的选项不在预设的范围内，请检查后重试！")
              }
            }
          }
          #---------------------------------------------------------------------
          pn = pn + 1
        }
        ########################################################################
        ########################################################################
      } else { # panel_corrplot == F
        print("注意：Olink_Corrplot函数的参数panel_corrplot = F，所以不进行Panel的PCA分析！")
      }

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
      print("暂时不能处理！")
    } else {
      print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
    }
  } else if (npx_data_type == "widetype") {
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {
      print("暂时不能处理！")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
      print("暂时不能处理！")
    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }

  ## 脚本运行结束
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(green(bold("--------------------------------[3.Olink项目OAAT2-QC-Corrplot: Olink_Corrplot()函数运行完成！]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))


}

## 执行函数
# Olink_Corrplot(path = path,
#                npx_data_type = npx_data_type,
#                folder_system = folder_system,
#                output_images_file_format = output_images_file_format,
#                qc_warning = qc_warning,
#                panel_corrplot = panel_corrplot,
#                firstcol = "UniProt",
#                firstrow = "SampleID",
#                cellvalue = "NPX",
#                sample_control_name = sample_control_name,
#                filtered_samples = filtered_samples,
#                corrplot_method = corrplot_method,
#                corrplot_type = corrplot_type,
#                corrplot_style = corrplot_style)
