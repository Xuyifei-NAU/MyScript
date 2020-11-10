#!/bin/bash
# 切换目录
  wd=`pwd`
  cd ${wd}/result/salmon
# 生成列表
  ls gene*.count > count.list
  ls gene*.TPM > TPM.list
# 生成name的文件
  for i in `cut -f 1 TPM.list`;do
    cut -f 1 ${i} > name_${i} 
    cut -f 2 ${i} > value_${i} 
  done

  for i in `cut -f 1 count.list`;do
    cut -f 1 ${i} > name_${i} 
    cut -f 2 ${i} > value_${i} 
  done
# 比较name顺序是否一样
  touch log
  for i in  `cut -f 1 TPM.list`;do
    diff name_`head -1 TPM.list` name_${i} >>log
  done

  for i in  `cut -f 1 count.list`;do
    diff name_`head -1 count.list` name_${i} >>log
  done
# 生成merge文件
  paste name_`head -1 TPM.list` `ls value*.TPM`>gene.TPM
  paste name_`head -1 count.list` `ls value*.count`>gene.count
# 删除temp文件,返回目录
  rm name*
  rm value*
  rm `cat *list`
  rm *list
  cd ${wd}
