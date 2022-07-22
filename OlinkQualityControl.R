################################################################################
################################################################################
## Rscript Function: Olink Data QuanlityControl
## Main function: olink_quality_control()
## Author: jingxinxing
## Date: 2022-03-02-
##       2022-07-21
## Version: 2.1.0
## Usage:
################################################################################
################################################################################

olink_quality_control <- function(path = "/Users/mingma_10000455/P/OLINK项目/伯豪Olink_data_20220627", # Olink NPX Data QC Path
                                  npx_data_type = "longtype", # has three types: rawtype, longtype, widetype
                                  npx_data_files = "伯豪_Expansion_Long_NPX_Data.xlsx", # NPX Data Object or in Path of NPX Data(raw/long) Excel Files(or other files format, like csv/tsv/xls/txt, and so on, all support)
                                  sample_grouping_info = "sample_grouping_information.xlsx", # Samples grouping info or other additional info, like QC_warning
                                  sample_control_name = "SC", # Sample Control/SAMPLE CONTROL/SC/sc/Sample_Control/SAMPLE_CONTROL
                                  output_images_file_format = "pdf", # the output images file format.
                                  # output_data_file_format = "xlsx", # the output data file format.
                                  olink_bridgeselect = F, # Logical parameters, Olink different batch sample plant data bridgeselect function, the arguments is logical.
                                  bridge_samples_number = 8, # If olink_bridgeselect = T, this arguments is work, otherwise is invalid.
                                  olink_normalization = F, # Logical parameters, Whether Normalizing To Olink NPX Data
                                  olink_project_type = "T96", # Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
                                  folder_system = F, # Logical parameters, Whether Need a Folder System
                                  heatmap_color_style = "light", # dark/light
                                  panel_heatmap = F, # Logical parameters, Whether Need Heatmap for every Panel
                                  panel_pca = F, # Logical parameters, Whether Plotting Panel PCA
                                  firstcol = "UniProt", # NPX Wide Data 1st Column
                                  firstrow = "SampleID", # NPX Wide Data 1st Row
                                  cellvalue = "NPX", # NPX Wide Data Cell Value
                                  filtered_samples = "SC", # Control Samples Regular Name
                                  qc_warning = T, # Whether Plotting QC_warning Samples in PCA as diff Shape
                                  panel_corrplot = F, # Logical parameters, Whether Plotting Panel Corrplot
                                  corrplot_method = "circle", # String parameter, Optional range:
                                  corrplot_type = "upper", # String parameter, Optional range: "full", "upper", "lower"
                                  corrplot_style = F # Logical parameters, Whether Use Default Corrplot Style
                                  ) {
  ## 设置当前工作路径
  setwd(path)
  print(paste0("[Olink项目NPX数据质控路径为：", path, " ]"))

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

  ## 1.脚本运行计时开始
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  pt <- proc.time()
  cat(bgBlue(bold(paste0("当前日期时间为：", date(), "\n"))))
  cat(bgWhite(red(bold("--------------------------------[0.开始运行Olink项目分析流程：OAAT2-QC]--------------------------------\n"))))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
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

  ## 3.Olink分析结果文件夹系统是否构建
  if (folder_system == T) { # 需要构建Olink分析结果文件夹
    #############
    ## RawData ##
    #############
    if (dir.exists(paste0(path,"/0.RawData"))) {
      # cat("\n")
    } else {
      dir.create(paste0(path,"/0.RawData"))
    }
    ####################
    ## QualityControl ##
    ####################
    if (dir.exists(paste0(path,"/1.QualityControl"))) {
      # cat("\n")
      ### 1.1 样本NPX分布箱线图
      if (dir.exists(paste0(path,"/1.QualityControl/1.1NPX-Samples_Boxplot"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/1.QualityControl/1.1NPX-Samples_Boxplot"))
      }
      ### 1.2 所有蛋白表达聚类热图
      if (dir.exists(paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap"))
      }
      ### 1.3 样本主成分分析PCA
      if (dir.exists(paste0(path,"/1.QualityControl/1.3Samples_PCA"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/1.QualityControl/1.3Samples_PCA"))
      }
      ### 1.4 样本相关性系数热图
      if (dir.exists(paste0(path,"/1.QualityControl/1.4Samples_Corrplot"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/1.QualityControl/1.4Samples_Corrplot"))
      }
      ### 1.5 OLINK NPX数据标准化处理
      if (dir.exists(paste0(path,"/1.QualityControl/1.5Normalization_NPX_Data"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/1.QualityControl/1.5Normalization_NPX_Data"))
      }
      #
    } else {
      dir.create(paste0(path,"/1.QualityControl"))
      ### 1.1 样本NPX分布箱线图
      dir.create(paste0(path,"/1.QualityControl/1.1NPX-Samples_Boxplot"))
      ### 1.2 所有蛋白表达聚类热图
      dir.create(paste0(path,"/1.QualityControl/1.2Proteins_Expression_Heatmap"))
      ### 1.3 样本主成分分析PCA
      dir.create(paste0(path,"/1.QualityControl/1.3Samples_PCA"))
      ### 1.4 样本相关性系数热图
      dir.create(paste0(path,"/1.QualityControl/1.4Samples_Corrplot"))
      ### 1.5 OLINK NPX数据标准化处理
      dir.create(paste0(path,"/1.QualityControl/1.5Normalization_NPX_Data"))
    }
    #######################
    ## ProteinExpression ##
    #######################
    if (dir.exists(paste0(path,"/2.ProteinExpression"))) {
      # cat("\n")
      ### 2.1 Olink Panel 蛋白表达聚类热图
      if (dir.exists(paste0(path,"/2.ProteinExpression/2.1Panel_Proteins_Expression_Heatmap"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/2.ProteinExpression/2.1Panel_Proteins_Expression_Heatmap"))
      }
      ### 2.2 Olink Panel 相关性系数热图
      if (dir.exists(paste0(path,"/2.ProteinExpression/2.2Panel_Proteins_Correlation_Coefficient_Heatmap"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/2.ProteinExpression/2.2Panel_Proteins_Correlation_Coefficient_Heatmap"))
      }
      ### 2.3 目标蛋白表达箱线图
      if (dir.exists(paste0(path,"/2.ProteinExpression/2.3Target_Protein_NPX_Boxplot"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/2.ProteinExpression/2.3Target_Protein_NPX_Boxplot"))
      }
    } else {
      dir.create(paste0(path,"/2.ProteinExpression"))
      ### 2.1 Olink Panel 蛋白表达聚类热图
      dir.create(paste0(path,"/2.ProteinExpression/2.1Panel_Proteins_Expression_Heatmap"))
      ### 2.2 Olink Panel 相关性系数热图
      dir.create(paste0(path,"/2.ProteinExpression/2.2Panel_Proteins_Correlation_Coefficient_Heatmap"))
      ### 2.3 目标蛋白表达箱线图
      dir.create(paste0(path,"/2.ProteinExpression/2.3Target_Protein_NPX_Boxplot"))
    }
    ###################################
    ## DifferentialExpressionProtein ##
    ###################################
    if (dir.exists(paste0(path,"/3.DifferentialExpressionProtein"))) {
      # cat("\n")
      ### 3.1 差异蛋白占比统计饼图
      if (dir.exists(paste0(path,"/3.DifferentialExpressionProtein/3.1Differential_Protein_Statistical_Pie"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.1Differential_Protein_Statistical_Pie"))
      }
      ### 3.2 Panel差异蛋白聚类热图
      if (dir.exists(paste0(path,"/3.DifferentialExpressionProtein/3.2Differential_Protein_Heatmap"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.2Differential_Protein_Heatmap"))
      }
      ### 3.3 差异蛋白火山图
      if (dir.exists(paste0(path,"/3.DifferentialExpressionProtein/3.3Differential_Protein_Volcano"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.3Differential_Protein_Volcano"))
      }
    } else {
      dir.create(paste0(path,"/3.DifferentialExpressionProtein"))
      ### 3.1 差异蛋白占比统计饼图
      dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.1Differential_Protein_Statistical_Pie"))
      ### 3.2 Panel差异蛋白聚类热图
      dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.2Differential_Protein_Heatmap"))
      ### 3.3 差异蛋白火山图
      dir.create(paste0(path,"/3.DifferentialExpressionProtein/3.3Differential_Protein_Volcano"))
    }
    ########################
    ## FunctionEnrichment ##
    ########################
    if (dir.exists(paste0(path,"/4.FunctionEnrichment"))) {
      ####
      # cat("\n")
      ### FunctionEnrichment/GO_Enrichment
      if (dir.exists(paste0(path,"/4.FunctionEnrichment/4.1GO_Enrichment"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/4.FunctionEnrichment/4.1GO_Enrichment"))
      }
      ### /FunctionEnrichment/KEGG_Enrichment
      if (dir.exists(paste0(path,"/4.FunctionEnrichment/4.2KEGG_Enrichment"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/4.FunctionEnrichment/4.2KEGG_Enrichment"))
      }
      ### /FunctionEnrichment/Reactome_Enrichment
      if (dir.exists(paste0(path,"/4.FunctionEnrichment/4.3Reactome_Enrichment"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/4.FunctionEnrichment/4.3Reactome_Enrichment"))
      }
      ####
    } else {
      ### FunctionEnrichment
      dir.create(paste0(path,"/4.FunctionEnrichment"))
      ### FunctionEnrichment/GO_Enrichment
      dir.create(paste0(path,"/4.FunctionEnrichment/4.1GO_Enrichment"))
      ### /FunctionEnrichment/KEGG_Enrichment
      dir.create(paste0(path,"/4.FunctionEnrichment/4.2KEGG_Enrichment"))
      ### /FunctionEnrichment/Reactome_Enrichment
      dir.create(paste0(path,"/4.FunctionEnrichment/4.3Reactome_Enrichment"))
    }
    ######################################
    ## BioinformaticsDatabaseAnnotation ##
    ######################################
    if (dir.exists(paste0(path,"/5.BioinformaticsDatabaseAnnotation"))) {
      # cat("\n")
      ### 5.1 蛋白注释表
      if (dir.exists(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.1Annotation"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.1Annotation"))
      }
      ### 5.2 蛋白互作分析（PPI）
      if (dir.exists(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.2PPI_Analysis"))) {
        # cat("\n")
      } else {
        dir.create(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.2PPI_Analysis"))
      }
    } else {
      dir.create(paste0(path,"/5.BioinformaticsDatabaseAnnotation"))
      ### 5.1 蛋白注释表
      dir.create(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.1Annotation"))
      ### 5.2 蛋白互作分析（PPI）
      dir.create(paste0(path,"/5.BioinformaticsDatabaseAnnotation/5.2PPI_Analysis"))
    }
    #############################################
    ## Olink Project Analysis Result Directory ##
    #############################################
    if (dir.exists(paste0(path,"/data_deliver"))) {
      # cat("\n")
    } else {
      dir.create(paste0(path,"/data_deliver"))
    }
    ## Ending
    cat(blue("======================================================================\n"))
    cat(bold(blue("OAAT Project Analyze Directory Has Been Generated!\n")))
    cat(blue("======================================================================\n"))
    print(dir(path, recursive = T))
  } else {
    cat(red("======================================================================\n"))
    print("请注意：不构建Olink分析结果文件夹系统！")
    cat(red("======================================================================\n"))
  }

  ## 4.读取Olink项目的NPX数据
  if (npx_data_type == "rawtype") { # npx_data_type has three types: rawtype, longtype, widetype

    ### 读取npx_data_files
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # 只有一个npx_data_files
      npx_data <- read_NPX(paste0(path, "/", npx_data_files))
      cat("============================================================================================================\n")
      cat(bold(red("[ rawtype npx_data: ]\n")))
      # print(head(npx_data))
      cat("------------------------------------------------------------------------------------------------------------\n")
      # print(summary(npx_data))
      cat("============================================================================================================\n")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # 两个npx_data_files
      #### 遍历所有npx数据文件
      for (n in 1:length(strsplit(npx_data_files, split = ",")[[1]])) {
        # print(paste0("读取Olink项目第", n, "个NPX原始数据..."))
        assign(paste("npx_data", n, sep = ""), read_NPX(paste0(path,"/", strsplit(npx_data_files, split = ",")[[1]][n])))
        cat("============================================================================================================\n")
        cat(bold(red(paste("[ rawtype npx_data", n, ": ]", sep = ""),"\n")))

        # print(head(paste("npx_data", n, sep = "")))
        cat("------------------------------------------------------------------------------------------------------------\n")
        # print(summary(paste("npx_data", n, sep = "")))
        cat("============================================================================================================\n")
      }
      # #### 合并常见的两批次NPX数据
      # P1 <- unique(npx_data1$Project)[1]
      # P2 <- unique(npx_data2$Project)[1]
      # ## Normalizing NPX Data
      # # Find overlapping samples
      # overlap_samples <- intersect(npx_data1$SampleID, npx_data2$SampleID) %>%
      #   data.frame() %>%
      #   filter(!str_detect(., 'CONTROL_SAMPLE')) %>% #Remove control samples
      #   pull(.)
      # # Perform Bridging normalization
      # olink_normalization_data <- olink_normalization(df1 = npx_data1,
      #                                                 df2 = npx_data2,
      #                                                 overlapping_samples_df1 = overlap_samples,
      #                                                 df1_project_nr = P1,
      #                                                 df2_project_nr = P2,
      #                                                 reference_project = P1)
      #
      # ## Saving Data
      # export(olink_normalization_data, file = paste0(path,"/Olink_NPX_Normalization_Data.xlsx"), overwrite = T)
      # npx_data <- olink_normalization_data # 分析流程中一直使用默认NPX数据变量名称：npx_data

    } else {
      print("多于2个NPX raw Data文件！暂时不能处理！")
    }

    ### 读取样本分组信息表
    sgi <- import(paste0(path, "/", sample_grouping_info))
    cat("**********************************************************************\n")
    cat(bold(green("[ sample_grouping_info: ]\n")))
    # print(head(sgi))
    cat("----------------------------------------------------------------------\n")
    # print(summary(sgi))
    cat("**********************************************************************\n")

  } else if (npx_data_type == "longtype") { # longtype
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # only one NPX Long Data file
      npx_data <- read.xlsx(paste0(path, "/", npx_data_files))
      cat("============================================================================================================\n")
      cat(bold(blue("[ longtype npx_data: ]\n")))
      # print(head(npx_data))
      cat("------------------------------------------------------------------------------------------------------------\n")
      # print(summary(npx_data))
      cat("============================================================================================================\n")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two NPX Long Data files
      #### 遍历所有npx数据文件
      for (n in 1:length(strsplit(npx_data_files, split = ",")[[1]])) {
        # print(paste0("读取Olink项目第", n, "个NPX Long Data..."))
        assign(paste("npx_data", n, sep = ""), read.xlsx(paste0(path,"/", strsplit(npx_data_files, split = ",")[[1]][n])))
        cat("============================================================================================================\n")
        cat(bold(blue(paste("[ longtype npx_data", n, ": ]", sep = ""),"\n")))
        # print(head(paste("npx_data", n, sep = "")))
        cat("------------------------------------------------------------------------------------------------------------\n")
        # print(summary(paste("npx_data", n, sep = "")))
        cat("============================================================================================================\n")
      }
    } else { # now, not support more than two NPX Long Data files
      print("多于2个NPX Long Data文件！暂时不能处理！")
    }

    # reading sample_grouping_information.xlsx file
    sgi <- import(paste0(path, "/", sample_grouping_info))
    cat("**********************************************************************\n")
    cat(bold(green("[ sample_grouping_info: ]\n")))
    # print(head(sgi))
    cat("----------------------------------------------------------------------\n")
    # print(summary(sgi))
    cat("**********************************************************************\n")

  } else if (npx_data_type == "widetype") { # widetype

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {
      npx_data <- read.xlsx(paste0(path, "/", npx_data_files))
      cat("============================================================================================================\n")
      cat(bold(cyan("[ widetype npx_data: ]\n")))
      # print(head(npx_data))
      cat("------------------------------------------------------------------------------------------------------------\n")
      # print(summary(npx_data))
      cat("============================================================================================================\n")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
      #### 遍历所有npx数据文件
      for (n in 1:length(strsplit(npx_data_files, split = ",")[[1]])) {
        # print(paste0("读取Olink项目第", n, "个NPX Long Data..."))
        assign(paste("npx_data", n, sep = ""), read.xlsx(paste0(path,"/", strsplit(npx_data_files, split = ",")[[1]][n])))
        cat("============================================================================================================\n")
        cat(bold(cyan(paste("[ widetype npx_data", n, ": ]", sep = ""),"\n")))
        # print(head(paste("npx_data", n, sep = "")))
        cat("------------------------------------------------------------------------------------------------------------\n")
        # print(summary(paste("npx_data", n, sep = "")))
        cat("============================================================================================================\n")
      }
    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }

    sgi <- import(paste0(path, "/", sample_grouping_info))
    cat("**********************************************************************\n")
    cat(bold(green("[ sample_grouping_info: ]\n")))
    # print(head(sgi))
    cat("----------------------------------------------------------------------\n")
    # print(summary(sgi))
    cat("**********************************************************************\n")

  } else { # exit the loop: break/next
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }

  ## 5.样本、Panel、蛋白箱线图Boxplot
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") { # longtype and rawtype
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one

      ## plotting boxplot
      olink_panel <- npx_data$Panel %>% unique()
      for (panel in olink_panel) { # 循环处理Panel绘制样本分布图
        npx_data %>%
          filter(Panel == panel) %>% # For this example only plotting one panel.
          olink_dist_plot() +
          theme(axis.text.x = element_text(size = (300/npx_data$SampleID %>% unique() %>% length() + npx_data$SampleID %>% unique() %>% length()/50 ) , angle = 90, colour = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 14)) +
          theme(legend.title = element_text(size = 18), legend.text = element_text(size = 16))# 横坐标的样本名称值是可以调节的，这里使用(300/npx_data1$SampleID %>% unique() %>% length() + npx_data1$SampleID %>% unique() %>% length()/50 )进行自适应调节
        # theme(axis.text.x = element_blank()) # Due to the number of samples one can remove the text or rotate it

        if (folder_system == T) { # folder_system == T
        ggsave(filename = paste0(path,"/1.QualityControl/1.1NPX-Samples_Boxplot/Olink_", panel, "_NPX_Samples_Distribution_Boxplot.", switch (output_images_file_format,
          "pdf" = "pdf",
          "png" = "png",
          "jpeg" = "jpeg",
          "tiff" = "tiff",
          "svg" = "svg",
          "bmp" = "bmp"
        )), width = npx_data$SampleID %>% unique() %>% length()/5, height = 8, units = "in")
        # ggsave(filename = paste0(path,"/1.QualityControl/1.1NPX-Samples_Boxplot/Olink_", panel, "_NPX_Samples_Distribution_Boxplot.png"), width = npx_data$SampleID %>% unique() %>% length()/10, height = 8, units = "in")
        } else { # folder_system == F
          ggsave(filename = paste0(path,"/Olink_", panel, "_Samples_NPX_Distribution_Boxplot.", switch (output_images_file_format,
                                                                                                        "pdf" = "pdf",
                                                                                                        "png" = "png",
                                                                                                        "jpeg" = "jpeg",
                                                                                                        "tiff" = "tiff",
                                                                                                        "svg" = "svg",
                                                                                                        "bmp" = "bmp"
          )), width = npx_data$SampleID %>% unique() %>% length()/5, height = 8, units = "in")
          # ggsave(filename = paste0(path,"/Olink_", panel, "_Samples_NPX_Distribution_Boxplot.png"), width = npx_data$SampleID %>% unique() %>% length()/10, height = 8, units = "in")
        }

        cat(bgYellow(paste0("Olink_", red(panel), "_NPX_Samples_Distribution_Boxplot Plotting Finished!"),"\n"))
      }

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two
      print("longtype and rawtype two npx data.")

    } else {
      print("多于2个NPX Raw/Long Data文件！暂时不能处理！")
    }
  } else if (npx_data_type == "widetype") { # widetype
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one
      print("widetype one npx data.")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two
      print("widetype two npx data.")
    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }

  ## 6.样本-Panel-蛋白热图Heatmap
  source("/glusterfs/home/local_jing_xx/software/OAAT2/Olink_Heatmap.R", encoding = "utf-8")
  Olink_Heatmap(path = path,
                npx_data_type = npx_data_type,
                npx_data_files = npx_data_files,
                # sample_grouping_info = sample_grouping_info,
                folder_system = F,
                output_images_file_format = output_images_file_format,
                # qc_warning = qc_warning,
                # firstcol = "UniProt",
                # firstrow = "SampleID",
                # cellvalue = "NPX",
                target_npx_data = NULL,
                heatmap_color_style = heatmap_color_style, # dark/light
                panel_heatmap = panel_heatmap, # Logical parameters, Whether Need Heatmap for every Panel
                olink_project_type = olink_project_type # Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
                # sample_control_name = sample_control_name,
                # filtered_samples = filtered_samples
  )
  #######################
  ## 7.样本主成分分析PCA ---------------------------------------------------- 7.
  #######################
  ### 7.1 2D PCA
  ### 7.2 3D PCA
  ### 7.3 PCA html
  source("/glusterfs/home/local_jing_xx/software//OAAT2/Olink_PCA.R", encoding = "utf-8")
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") { # rawtype or longtype

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one
      #### 执行Olink_PCA()分析绘图函数
      Olink_PCA(path = path,
                npx_data_type = npx_data_type,
                npx_data_files = npx_data_files,
                sample_grouping_info = sample_grouping_info,
                folder_system = folder_system,
                output_images_file_format = output_images_file_format,
                qc_warning = qc_warning,
                panel_pca = panel_pca,
                firstcol = firstcol,
                firstrow = firstrow,
                cellvalue = cellvalue,
                sample_control_name = sample_control_name,
                filtered_samples = filtered_samples)

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two
      print("暂时不支持！")

    } else {
      print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
    }
  } else if (npx_data_type == "widetype") { # widetype
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) { # one
      print("暂时不支持！")
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) { # two
      print("暂时不支持！")
    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }
  ############################
  ## 8.样本相关性分析Corrplot ----------------------------------------------- 8.
  ############################
  source(file = "/glusterfs/home/local_jing_xx/software/OAAT2/Olink_Corrplot.R", encoding = "utf-8")
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") {

    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {
      Olink_Corrplot(path = path,
                     npx_data_type = npx_data_type,
                     npx_data_files = npx_data_files,
                     folder_system = folder_system,
                     output_images_file_format = output_images_file_format,
                     qc_warning = qc_warning,
                     panel_corrplot = panel_corrplot,
                     firstcol = firstcol,
                     firstrow = firstrow,
                     cellvalue = cellvalue,
                     sample_control_name = sample_control_name,
                     filtered_samples = filtered_samples,
                     corrplot_method = corrplot_method,
                     corrplot_type = corrplot_type,
                     corrplot_style = corrplot_style)
    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
      print("暂时不能处理！")
    } else {
      print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
    }

  } else if (npx_data_type == "widetype") {
    if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {

    } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
      print("暂时不能处理！")
    } else {
      print("多于2个NPX Wide Data文件！暂时不能处理！")
    }
  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }
  ## 9.Normalization_NPX_Data
  if (npx_data_type == "rawtype" || npx_data_type == "longtype") {

  } else if (npx_data_type == "widetype") {

  } else {
    print("暂时只支持NPX数据类型：rawtype, longtype, widetype！")
  }
  #
  if (length(strsplit(npx_data_files, split = ",")[[1]]) >= 2) {
    ## Saving Data
    rio::export(olink_normalization_data, file = paste0(path,"/Olink_NPX_Normalization_Data.xlsx"), overwrite = T)
  }
  ## 10.桥接样本
  if (olink_bridgeselect == T) { # 需要获取NPX数据的桥接样本
    cat(bold(red("[需要获取NPX数据的桥接样本：]\n")))

    if (folder_system == T) { # 需要Olink分析文件夹系统
      if (npx_data_type == "rawtype" || npx_data_type == "longtype") { # RawData/LongData
        if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {

          if (!is.null(bridge_samples_number) && bridge_samples_number <= length(unique(npx_data$SampleID))) {
            sample_num <- length(unique(npx_data$SampleID))
            olink_bridgeselector_sample_df <- olink_bridgeselector(df = npx_data, sampleMissingFreq = 0.3, n = bridge_samples_number)
            rio::export(olink_bridgeselector_sample_df, file = paste0(path,"/1.QualityControl/1.5Normalization_NPX_Data/Olink_Bridgeselector_Samples_Data.xlsx"), overwrite = T, zoom = 120)
          } else {
            sample_num <- length(unique(npx_data$SampleID))
            olink_bridgeselector_sample_df <- olink_bridgeselector(df = npx_data, sampleMissingFreq = 0.3, n = sample_num/10)
            rio::export(olink_bridgeselector_sample_df, file = paste0(path,"/1.QualityControl/1.5Normalization_NPX_Data/Olink_Bridgeselector_Samples_Data.xlsx"), overwrite = T, zoom = 120)
          }

        } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
          print("暂时不能处理！")
        } else {
          print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
        }
      } else if (npx_data_type == "widetype") { # WideData
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
    } else { # 不需要Olink文件夹系统
      if (npx_data_type == "rawtype" || npx_data_type == "longtype") { # RawData/LongData
        if (length(strsplit(npx_data_files, split = ",")[[1]]) == 1) {
          if (!is.null(bridge_samples_number) && bridge_samples_number <= length(unique(npx_data$SampleID))) {
            sample_num <- length(unique(npx_data$SampleID))
            olink_bridgeselector_sample_df <- olink_bridgeselector(df = npx_data, sampleMissingFreq = 0.3, n = bridge_samples_number)
            rio::export(olink_bridgeselector_sample_df, file = paste0(path,"/Olink_Bridgeselector_Samples_Data.xlsx"), overwrite = T, zoom = 120)
          } else {
            sample_num <- length(unique(npx_data$SampleID))
            olink_bridgeselector_sample_df <- olink_bridgeselector(df = npx_data, sampleMissingFreq = 0.3, n = sample_num/10)
            rio::export(olink_bridgeselector_sample_df, file = paste0(path,"/Olink_Bridgeselector_Samples_Data.xlsx"), overwrite = T, zoom = 120)
          }
        } else if (length(strsplit(npx_data_files, split = ",")[[1]]) == 2) {
          print("暂时不能处理！")
        } else {
          print("多于2个NPX Long/Raw Data文件！暂时不能处理！")
        }
      } else if (npx_data_type == "widetype") { # WideData
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
    }

  } else {
    print("注意：不需要获取NPX数据的桥接样本。")
  }
  ## 11.脚本耗时 ##
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(blue("[脚本运行时间为]："))
  print(proc.time() - pt)
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
  cat(bgYellow("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
}

###########################################################
## 执行olink_quality_control()函数 ##
###########################################################
# olink_quality_control(path = "/Users/mingma_10000455/P/OLINK项目/医创云康", # Olink NPX Data QC Path
#                     npx_data_type = "rawtype", # has three types: rawtype, longtype, widetype
#                     npx_data_files = "医创云康88血浆_Raw_NPX_Data.xlsx", # NPX Data Object or in Path of NPX Data(raw/long) Excel Files(or other files format, like csv/tsv/xls/txt, and so on, all support)
#                     sample_grouping_info = "sample_grouping_information.xlsx", # Samples grouping info or other additional info, like QC_warning
#                     sample_control_name = "Sample Control", # Sample Control/SAMPLE CONTROL/SC/sc/Sample_Control/SAMPLE_CONTROL
#                     output_images_file_format = "pdf", # the output images file format.
#                     # output_data_file_format = "xlsx", # the output data file format.
#                     olink_bridgeselect = F, # Olink different batch sample plant data bridgeselect function, the arguments is logical.
#                     bridge_samples_number = 8, # If olink_bridgeselect = T, this arguments is work, otherwise is invalid.
#                     olink_normalization = F,
#                     olink_project_type = "T96", # Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
#                     folder_system = F, # Whether Need a Folder System
#                     heatmap_color_style = "light", # dark/light
#                     panel_pca = T, # Whether Need PCA for every Panel
#                     panel_heatmap = T, # Whether Need Heatmap for every Panel
#                     firstcol = "UniProt", # NPX Wide Data 1st Column
#                     firstrow = "SampleID", # NPX Wide Data 1st Row
#                     cellvalue = "NPX", # NPX Wide Data Cell Value
#                     filtered_samples = "Sample Control", # Control Samples Regular Name
#                     qc_warning = T, # Whether Plotting QC_warning Samples in PCA as diff Shape
#                     panel_corrplot = T, #
#                     corrplot_method = "circle", #
#                     corrplot_type = "upper", #
#                     corrplot_style = T # Corrplot Style
#                     )
