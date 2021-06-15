############################
## Author: Deepak Tanwar  ##
## Date: 20191008	  ##
############################

# This makefile will run the TrimGalore software to improve the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*R1.fq.gz)
targets := $(patsubst $(source_dir)/%_R1.fq.gz, $(target_dir)/%, $(files))

all: $(targets)

$(target_dir)/%: $(source_dir)/%_R1.fq.gz
	mkdir $@ && trim_galore -j 5 --paired -q 30 --length 30 -o $@ --trim-n --retain_unpaired $< ./input/$(basename $(basename $(notdir $@)))_R2.fq.gz > ./log/$(notdir $@).log 2>>./log/$(notdir $@).log
