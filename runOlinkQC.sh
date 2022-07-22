#!/bin/bash
echo `date`
echo "开始执行Olink项目数据质控流程：[ OAAT2-QC ]"
/glusterfs/home/local_jing_xx/anaconda3/bin/Rscript /glusterfs/home/local_jing_xx/software/OAAT2/runOlinkQualityControl.R
