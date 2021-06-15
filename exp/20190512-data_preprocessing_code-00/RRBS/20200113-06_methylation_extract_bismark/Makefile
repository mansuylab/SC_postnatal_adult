############################
## Author: Deepak Tanwar  ##
## Date: 20200113	  ##
############################

# This makefile will extract the methylation information from bismark alignment files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_trimmed_bismark_bt2.bam)
targets := $(patsubst $(source_dir)/%_trimmed_bismark_bt2.bam, $(target_dir)/%, $(files))

all: $(targets)

$(target_dir)/%: $(source_dir)/%_trimmed_bismark_bt2.bam
	mkdir -p $@ && bismark_methylation_extractor --version > ./log/$(notdir $@).log && bismark_methylation_extractor $< -s -o $@ --multicore 8 --bedGraph --cytosine_report --genome_folder ./input/mm10/ 2>> ./log/$(notdir $@).log >> ./log/$(notdir $@).log && pigz -3 -p 8 $@/*.txt
