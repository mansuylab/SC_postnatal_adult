##############################
## Author: Deepak Tanwar    ##
## Date: 20191210	    ##
##############################

# This makefile will run alignment (subjunc) with Rsubread

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1_trimmed.fq.gz)
targets := $(patsubst $(source_dir)/%_1_trimmed.fq.gz, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_1_trimmed.fq.gz
	 ./20191210-alignment_Rsubread.R $< 2>./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log
