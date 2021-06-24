##############################
## Author: Deepak Tanwar    ##
## Date: 20210617	    ##
##############################

# This makefile will run alignment (subjunc) with Rsubread

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_trimmed.fq.gz)
targets := $(patsubst $(source_dir)/%_trimmed.fq.gz, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_trimmed.fq.gz
	./alignment_Rsubread.R $< 2>./log/$(basename $(notdir $@)).stderr > ./log/$(basename $(notdir $@)).stdout && rm $@.bai && samtools index $@
