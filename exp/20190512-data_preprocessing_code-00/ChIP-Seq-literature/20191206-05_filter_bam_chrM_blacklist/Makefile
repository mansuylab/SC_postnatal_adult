############################
## Author: Deepak Tanwar  ##
## Date: 20191206	  ##
############################

# This makefile will run the deepTools for filtering bamfile for blacklist regions and remove chrM.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%.bam
	alignmentSieve --version 2> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log && alignmentSieve --blackListFileName ./input/atac_blacklist_mm10.bed --minMappingQuality 30 --filterMetrics $(basename $@).filterMetric.txt -p 16 -b $< -o $(basename $@)_tmp.bam 2>> ./log/$(basename $(notdir $@)).log && samtools sort -@ 16 $(basename $@)_tmp.bam $(basename $@) && rm -rf $(basename $@)_tmp.bam && samtools index $@ && rm -rf /tmp/_deep*
