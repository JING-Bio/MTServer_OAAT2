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
  rio::export(npx_widedata_long, file = paste0(output_file_path,"/Olink_NPX_LongData.xlsx"))
  cat(blue("*****************************************************************\n"))
  cat(bgYellow("[Olink项目NPX数据宽变长转换完成！]\n"))
  cat(bgRed(paste0("宽数据路径为：\n", bold(output_file_path,"/Olink_NPX_LongData.xlsx\n"))))
  cat(red("*****************************************************************\n"))
}

olink_wide2long(npx_widedata_path = "/Users/mingma_10000455/P/OLINK项目/医创云康/医创云康88血浆_Raw_NPX_Data.xlsx",
                            target_npx_widedata = NULL,
                            firstcol = "UniProt",
                            firstrow = "SampleID",
                            cellvalue = "NPX",
                            output_file_path = "/Users/mingma_10000455/P/OLINK项目/医创云康/")
