############################
## Author: Deepak Tanwar  ##
## Date: 20190816	  ##
############################

# This makefile will mark duplicated reads using picard tools.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.dupPicard.bam, $(files))

all: $(targets)

$(target_dir)/%.dupPicard.bam: $(source_dir)/%.bam
	picard MarkDuplicates I=$< O=$@ M=$(basename $@).txt 2>./log/$(notdir $(basename $<)).log
