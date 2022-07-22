#!//glusterfs/home/local_jing_xx/anaconda3/bin/R
source("/glusterfs/home/local_jing_xx/software/OAAT2/OlinkQualityControl.R", encoding = "utf-8")
## 获取外部参数：sample_grouping_information.xlsx  伯豪_Expansion_Long_NPX_Data.xlsx
args = commandArgs(TRUE)
if( length(args) != 6) {
        cat("Useage: /glusterfs/home/local_jing_xx/anaconda3/bin/Rscript  runOlinkQualityControl.R [*NPX_Data.xlsx] [sample_grouping_information.xlsx] [longtype] [SC] [pdf] [SC]")
        q()
}
### 参数传递
npx_data_files = args[1]
sample_grouping_info = args[2]
npx_data_type = args[3]
sample_control_name = args[4]
output_images_file_format = args[5]
filtered_samples = args[6]

## 执行olink_quality_control()函数
path <- getwd()
olink_quality_control(path = path, # Olink NPX Data QC Path
                      npx_data_type = npx_data_type, # has three types: rawtype, longtype, widetype
                      npx_data_files = npx_data_files, # NPX Data Object or in Path of NPX Data(raw/long) Excel Files(or other files format, like csv/tsv/xls/txt, and so on, all support)
                      sample_grouping_info = sample_grouping_info, # Samples grouping info or other additional info, like QC_warning
                      sample_control_name = sample_control_name, # Sample Control/SAMPLE CONTROL/SC/sc/Sample_Control/SAMPLE_CONTROL
                      output_images_file_format = output_images_file_format, # the output images file format.
                      # output_data_file_format = "xlsx", # the output data file format.
                      olink_bridgeselect = F, # Olink different batch sample plant data bridgeselect function, the arguments is logical.
                      bridge_samples_number = 8, # If olink_bridgeselect = T, this arguments is work, otherwise is invalid.
                      olink_normalization = F,
                      olink_project_type = "T96", # Target 48/Target 96(T48/T96)/Expore 384(E384)/Expore 1536(E1536)/Expore 3072(E3072)
                      folder_system = T, # Whether Need a Folder System
                      heatmap_color_style = "light", # dark/light
                      panel_pca = T, # Whether Need PCA for every Panel
                      panel_heatmap = T, # Whether Need Heatmap for every Panel
                      firstcol = "UniProt", # NPX Wide Data 1st Column
                      firstrow = "SampleID", # NPX Wide Data 1st Row
                      cellvalue = "NPX", # NPX Wide Data Cell Value
                      filtered_samples = filtered_samples, # Control Samples Regular Name
                      qc_warning = T, # Whether Plotting QC_warning Samples in PCA as diff Shape
                      panel_corrplot = T, #
                      corrplot_method = "circle", #
                      corrplot_type = "upper", #
                      corrplot_style = T # Corrplot Style
                      )

