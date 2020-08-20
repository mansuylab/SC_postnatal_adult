############################
## Author: Deepak Tanwar  ##
## Date: 20200414	  ##
############################

# This makefile will sort the bam files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1_val_1_bismark_bt2_pe.bam)
targets := $(patsubst $(source_dir)/%_1_val_1_bismark_bt2_pe.bam, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_1_val_1_bismark_bt2_pe.bam
	samtools sort -@ 16 $< $(basename $@) && samtools index $@
