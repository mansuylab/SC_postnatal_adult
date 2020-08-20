###########################
## Author: Deepak Tanwar ##
## Date: 20191008        ##
###########################

## This makefile will run MACS2 for peak-calling.

source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%, $(files))

all: $(targets)

$(target_dir)/%: $(source_dir)/%.bam
	mkdir $@ && macs2 callpeak -g 2744254612 -f BAMPE --outdir $@ --name $(notdir $@) --tempdir $@ -t $< 2> ./log/$(notdir $@).log
