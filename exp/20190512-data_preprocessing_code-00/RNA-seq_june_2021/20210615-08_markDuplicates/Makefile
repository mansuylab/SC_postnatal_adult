############################
## Author: Deepak Tanwar  ##
## Date: 20210623	  ##
############################

# This makefile will mark duplicated reads using picard tools.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.markDup.bam, $(files))

all: $(targets)

$(target_dir)/%.markDup.bam: $(source_dir)/%.bam
	sambamba markdup -t 4 -p $< $@ 2>./log/$(notdir $(basename $<)).stderr > ./log/$(notdir $(basename $<)).stdout
