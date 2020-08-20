############################
## Author: Deepak Tanwar  ##
## Date: 20180720	  ##
############################

# This makefile will extract the methylation information from bismark alignment files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1_val_1_bismark_bt2_pe.bam)
targets := $(patsubst $(source_dir)/%_1_val_1_bismark_bt2_pe.bam, $(target_dir)/%, $(files))

all: $(targets)

$(target_dir)/%: $(source_dir)/%_1_val_1_bismark_bt2_pe.bam
	mkdir -p $@ && bismark_methylation_extractor --version > ./log/$(notdir $@).log && bismark_methylation_extractor $< -s -o $@ --multicore 8 --bedGraph --cytosine_report --genome_folder ./input/mm10/ 2>> ./log/$(notdir $@).log >> ./log/$(notdir $@).log && pigz -3 -p 8 $@/*.txt
