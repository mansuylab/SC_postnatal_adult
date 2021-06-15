############################
## Author: Deepak Tanwar  ##
## Date: 20200414	  ##
############################

# This makefile will sort the bam files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_trimmed_bismark_bt2.bam)
targets := $(patsubst $(source_dir)/%_trimmed_bismark_bt2.bam, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_trimmed_bismark_bt2.bam
	samtools sort -@ 16 $< $(basename $@) && samtools index $@
