## R长数据和宽数据相互转换（数据框）自定义工具

library("rio")
library("openxlsx")
library("data.table")
library("visdat")

setwd("/Users/mingma_10000455/R_Project/R_Custom_Function/长宽数据相互转换")

npx_data <- import("/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/东方肝胆医院_Expansion_NPX.xlsx")

npx_data_new <- npx_data[,c("SampleID", "UniProt", "NPX")]

head(npx_data_new)
# SampleID UniProt     NPX
# 1       SC  Q8IZC4  0.9054
# 2       SC  P78524  0.2475
# 3       SC  Q9H2M3  1.0249
# 4       SC  P55769 -0.3033
# 5       SC  Q9Y2W1 -0.0610
# 6       SC  O43734 -0.1381

npx_data_new2 <- npx_data_new[which(npx_data_new$SampleID != "SC"),] # Sample Control/SC

npx_clean_longdata <- npx_data_new2[which(npx_data_new2$SampleID != "SC-2"),]

head(npx_clean_longdata)
# SampleID UniProt     NPX
# 735 C21002368  Q8IZC4  0.3482
# 736 C21002368  P78524 -0.3717
# 737 C21002368  Q9H2M3 -0.0604
# 738 C21002368  P55769 -0.1202
# 739 C21002368  Q9Y2W1 -0.0728
# 740 C21002368  O43734  0.0268
dim(npx_clean_longdata)
# [1] 129448      3
summary(npx_clean_longdata)
# SampleID           UniProt               NPX
# Length:129448      Length:129448      Min.   :-9.2970
# Class :character   Class :character   1st Qu.:-0.2419
# Mode  :character   Mode  :character   Median : 0.2425
# Mean   : 0.2938
# 3rd Qu.: 0.7697
# Max.   :10.4329

samplelist <- unique(npx_clean_longdata$SampleID)
length(samplelist)
# [1] 88

proteinlist <- unique(npx_clean_longdata$UniProt)
length(proteinlist)
# [1] 1462

################################################################################
################################################################################
## 执行函数：olink_long2wide()
olink_long2wide(raw_npx_data_path = "/Users/mingma_10000455/P/OLINK项目/药明康德Olink蛋白检测数据_20220308/药明康德-细胞裂解物_NPX.xlsx",
                npx_data_path = NULL,
                target_npx_data = NULL,
                firstcol = "UniProt",
                firstrow = "SampleID",
                cellvalue = "NPX",
                filtered_samples = "Sample Control", # SC/Sample Control/CONTROL_SAMPLE_US_CS_AS_2
                output_file_path = "/Users/mingma_10000455/P/OLINK项目/药明康德Olink蛋白检测数据_20220308/")

## olink_long2wide()
olink_long2wide <- function(raw_npx_data_path = NULL,
                            npx_data_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/东方肝胆医院_Expansion_NPX.xlsx",
                            target_npx_data = NULL,
                            firstcol = "UniProt",
                            firstrow = "SampleID",
                            cellvalue = "NPX",
                            filtered_samples = "SC", # SC/Sample Control
                            output_file_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/") {
  ## 安装和导入R包
  analibpackages <- c("rio", "openxlsx", "crayon",
                      "magrittr", "reshape2",
                      "dplyr", "data.table", "tibble",
                      "parallel", "tidyr", "purrr",
                      "tidyverse", "stringr")
  for (a in analibpackages) {
    if (!requireNamespace(a, quietly = T)) {
      BiocManager::install(a)
    } else {
      #print(paste0("[软件包: '", a, "'已经安装完毕！]"))
      #print(packageVersion(a))
    }
  }

  sapply(analibpackages, library, character.only = T, quietly = T)

  ### devtools ###
  if (!requireNamespace("devtools", quietly = T)) {
    install.packages("devtools")
  } else {
    library("devtools", quietly = T)
  }


  ### OlinkAnalyze ###
  if (!requireNamespace("OlinkAnalyze", quietly = T)) {
    devtools::install_github(repo = 'Olink-Proteomics/OlinkRPackage/OlinkAnalyze',
                             build_vignettes = TRUE)
  } else {
    library("OlinkAnalyze", quietly = T)
  }
  cat(blue("*****************************************************************\n"))
  cat(bgYellow("[开始执行Olink项目长-宽数据转换任务：]\n"))
  cat(red("*****************************************************************\n"))
  ## 根据条件选择要读取的数据（是成熟的还是原始的NPX数据）
  if (is.null(target_npx_data) && is.null(raw_npx_data_path)  && !is.null(npx_data_path)) { # 需要处理的是成熟的NPX数据，NPX长数据Long Data
    ### rio::import()
    npx_data <- import(npx_data_path)
    ### data.frame subset
    npx_data_new <- npx_data[,c(firstrow, firstcol, cellvalue)]
    ### head df
    print(head(npx_data_new))
    ### filtered samples
    filtered_samples_list <- unique(npx_data_new[[firstrow]][grepl(filtered_samples, npx_data_new[[firstrow]], ignore.case = F, perl = F, fixed = F, useBytes = F)])

    ### 处理对照样本
    if (length(filtered_samples_list) == 1) { # Sample Control
      npx_clean_longdata <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list),]
    } else if (length(filtered_samples_list) == 2) { # [1] "SC"   "SC-2"
      npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list[1]),]
      npx_clean_longdata <- npx_data_new2[which(npx_data_new2[[firstrow]] != filtered_samples_list[2]),]
    } else if (length(filtered_samples_list) >= 3) { # 1, 2, 3, 4, ...
      df_temp <- data.frame()
      for (f in filtered_samples_list) {
        npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != f),]
        df_temp <- rbind(df_temp, npx_data_new2)
      }
      npx_clean_longdata <- df_temp
    } else {
      break
    }

    ###
  } else  if (is.null(target_npx_data) && !is.null(raw_npx_data_path)  && is.null(npx_data_path)) { # 需要处理原始的NPX数据
    ### OlinkAnalyze::read_NPX
    npx_data <- read_NPX(raw_npx_data_path)
    ### data.frame subset
    npx_data_new <- npx_data[,c(firstrow, firstcol, cellvalue)]
    ### head df
    print(head(npx_data_new))
    ### filtered samples
    filtered_samples_list <- unique(npx_data_new[[firstrow]][grepl(filtered_samples, npx_data_new[[firstrow]], ignore.case = F, perl = F, fixed = F, useBytes = F)])
    ### deal with control sample
    if (length(filtered_samples_list) == 1) { # Sample Control
      npx_clean_longdata <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list),]
    } else if (length(filtered_samples_list) == 2) { # [1] "SC"   "SC-2"
      npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list[1]),]
      npx_clean_longdata <- npx_data_new2[which(npx_data_new2[[firstrow]] != filtered_samples_list[2]),]
    } else if (length(filtered_samples_list) >= 3) { # 1, 2, 3, 4, ...
      df_temp <- data.frame()
      for (f in filtered_samples_list) {
        npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != f),]
        df_temp <- rbind(df_temp, npx_data_new2)
      }
      npx_clean_longdata <- df_temp
    } else {
      break
    }

  } else { # 直接用环境中的NPX对象数据：
    # target_npx_data
    ### data.frame subset
    npx_data_new <- target_npx_data[,c(firstrow, firstcol, cellvalue)]
    ### head df
    print(head(npx_data_new))
    ### filtered samples
    filtered_samples_list <- unique(npx_data_new[[firstrow]][grepl(filtered_samples, npx_data_new[[firstrow]], ignore.case = F, perl = F, fixed = F, useBytes = F)])

    ### 处理对照样本
    if (length(filtered_samples_list) == 1) { # Sample Control
      npx_clean_longdata <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list),]
    } else if (length(filtered_samples_list) == 2) { # [1] "SC"   "SC-2"
      npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != filtered_samples_list[1]),]
      npx_clean_longdata <- npx_data_new2[which(npx_data_new2[[firstrow]] != filtered_samples_list[2]),]
    } else if (length(filtered_samples_list) >= 3) { # 1, 2, 3, 4, ...
      df_temp <- data.frame()
      for (f in filtered_samples_list) {
        npx_data_new2 <- npx_data_new[which(npx_data_new[[firstrow]] != f),]
        df_temp <- rbind(df_temp, npx_data_new2)
      }
      npx_clean_longdata <- df_temp
    } else {
      print("没有检测到对照样本，不做处理！")
      npx_clean_longdata <- npx_data_new
    }

  }
  ### --------------------------------------------------------------------------
  ### summary
  print(summary(npx_clean_longdata))
  # SampleID           UniProt               NPX
  # Length:129448      Length:129448      Min.   :-9.2970
  # Class :character   Class :character   1st Qu.:-0.2419
  # Mode  :character   Mode  :character   Median : 0.2425
  # Mean   : 0.2938
  # 3rd Qu.: 0.7697
  # Max.   :10.4329
  ### --------------------------------------------------------------------------
  samplelist <- unique(npx_clean_longdata$SampleID)
  print(paste0("Olink-NPX数据样本个数为：", length(samplelist)))

  proteinlist <- unique(npx_clean_longdata$UniProt)
  print(paste0("Olink-NPX数据蛋白个数为：", length(proteinlist)))

  firstcol_list <- samplelist

  data_part_df <- npx_clean_longdata[which(npx_clean_longdata[[firstrow]] == npx_clean_longdata[[firstrow]][1]),]
  colnames(data_part_df)[3] <- npx_clean_longdata[[firstrow]][1]
  df <- data_part_df[!duplicated(data_part_df[[firstcol]]),]
  rownames(df) <- df[[firstcol]]

  df_new <- as.data.frame(df$UniProt, stringsAsFactors = F)
  colnames(df_new) <- "UniProt"

  # df <- data.frame()
  for (fc in firstcol_list) {
    data_part <- npx_clean_longdata[which(npx_clean_longdata[[firstrow]] == fc),]
    colnames(data_part)[3] <- fc
    data_part_new <- data_part[!duplicated(data_part[[firstcol]]),]
    rownames(data_part_new) <- data_part_new[[firstcol]]
    data_part_new2 <- as.data.frame(data_part_new[,c(fc)], stringsAsFactors = FALSE)
    colnames(data_part_new2) <- fc
    df_new <- cbind(df_new, data_part_new2)
    print(dim(df_new))
  }

  # print(dim(df_new))
  if (!is.null(npx_data_path) && is.null(raw_npx_data_path) && is.null(target_npx_data)) { # npx_data_path不为空，即读取的是NPX的长数据
  output_wide_npx_file <- strsplit(basename(npx_data_path), "[.]")[[1]][1]
  #### 获取NPX宽数据的行数（蛋白、基因）和列数（样本）
  npx_wide_nrow <- NROW(df_new)
  npx_wide_ncol <- NCOL(df_new)
  export(df_new, file = paste0(output_file_path,"/", output_wide_npx_file, "_WideData_", npx_wide_nrow, "x", npx_wide_ncol, ".xlsx"))
  } else if (is.null(npx_data_path) && !is.null(raw_npx_data_path) && is.null(target_npx_data)) { # raw_npx_data_path不为空，即读取的是NPX的raw数据（非长数据或宽数据）
    output_wide_npx_file <- strsplit(basename(raw_npx_data_path), "[.]")[[1]][1]
    #### 获取NPX宽数据的行数（蛋白、基因）和列数（样本）
    npx_wide_nrow <- NROW(df_new)
    npx_wide_ncol <- NCOL(df_new)
    export(df_new, file = paste0(output_file_path,"/", output_wide_npx_file, "_WideData_", npx_wide_nrow, "x", npx_wide_ncol, ".xlsx"))

  } else if (is.null(npx_data_path) && is.null(raw_npx_data_path) && !is.null(target_npx_data)) { # 使用R环境中的npx数据
    #### 获取NPX宽数据的行数（蛋白、基因）和列数（样本）
    npx_wide_nrow <- NROW(df_new)
    npx_wide_ncol <- NCOL(df_new)
    export(df_new, file = paste0(output_file_path,"/Olink_NPX_WideData_", npx_wide_nrow, "x", npx_wide_ncol,".xlsx"))
  } else {
    break
  }
  #
  cat(blue("***************************************************************************************\n"))
  cat(bgYellow("[Olink项目NPX数据长变宽转换完成！]\n"))
  cat(bgRed(paste0("Olink宽数据路径为：\n", bold(output_file_path,"/*_WideData.xlsx\n"))))
  cat(red("****************************************************************************************\n"))
}

################################################################################
################################################################################
## 执行函数：olink_wide2long
olink_wide2long(npx_widedata_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/Olink_NPX_WideData.xlsx",
                            target_npx_widedata = NULL,
                            firstcol = "UniProt",
                            firstrow = "SampleID",
                            cellvalue = "NPX",
                            output_file_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/")

## olink_wide2long
olink_wide2long <- function(npx_widedata_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/Olink_NPX_WideData.xlsx",
                            target_npx_widedata = NULL,
                            firstcol = "UniProt",
                            firstrow = "SampleID",
                            cellvalue = "NPX",
                            output_file_path = "/Users/mingma_10000455/P/OLINK项目/东方肝胆医院Olink项目数据质控/") {
  ## 安装和导入R包
  ## 安装和导入R包
  analibpackages <- c("rio", "openxlsx", "crayon",
                      "magrittr", "reshape2",
                      "dplyr", "data.table", "tibble",
                      "parallel", "tidyr", "purrr",
                      "tidyverse", "stringr")
  for (a in analibpackages) {
    if (!requireNamespace(a, quietly = T)) {
      BiocManager::install(a)
    } else {
      #print(paste0("[软件包: '", a, "'已经安装完毕！]"))
      #print(packageVersion(a))
    }
  }

  sapply(analibpackages, library, character.only = T, quietly = T)

  ### devtools ###
  if (!requireNamespace("devtools", quietly = T)) {
    install.packages("devtools")
  } else {
    library("devtools", quietly = T)
  }


  ### OlinkAnalyze ###
  if (!requireNamespace("OlinkAnalyze", quietly = T)) {
    devtools::install_github(repo = 'Olink-Proteomics/OlinkRPackage/OlinkAnalyze',
                             build_vignettes = TRUE)
  } else {
    library("OlinkAnalyze", quietly = T)
  }
  cat(blue("*****************************************************************\n"))
  cat(bgYellow("[开始执行Olink项目宽-长数据转换任务：]\n"))
  cat(red("*****************************************************************\n"))
  ## 读取NPX宽数据
  if (!is.null(npx_widedata_path) && is.null(target_npx_widedata)) { #  读取NPX宽数据文件
    npx_widedata <- import(npx_widedata_path)

  } else if (is.null(npx_widedata_path) && !is.null(target_npx_widedata)) { #
    npx_widedata <- target_npx_widedata
  } else {
    next
  }
  ## 执行宽-长数据转换
  npx_widedata_long <- melt(npx_widedata,
                            id.vars = "UniProt",
                            value.name = "NPX")
  colnames(npx_widedata_long)[2] <- "SampleID"

  ## 保存数据
  export(npx_widedata_long, file = paste0(output_file_path,"/Olink_NPX_LongData.xlsx"))
  cat(blue("*****************************************************************\n"))
  cat(bgYellow("[Olink项目NPX数据宽变长转换完成！]\n"))
  cat(bgRed(paste0("宽数据路径为：\n", bold(output_file_path,"/Olink_NPX_LongData.xlsx\n"))))
  cat(red("*****************************************************************\n"))
}
