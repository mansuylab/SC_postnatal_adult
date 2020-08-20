############################
## Author: Deepak Tanwar  ##
## Date: 20190816         ##
############################

# This makefile will run featureCounts with Rsubread

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.RData, $(files))

all: $(targets)

$(target_dir)/%.RData: $(source_dir)/%.bam
	./20190816-featureCounts_Rsubread.R $< 8 2>./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log
