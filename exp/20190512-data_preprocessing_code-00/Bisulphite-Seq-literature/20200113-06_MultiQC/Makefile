##############################
## @Author : Deepak Tanwar  ##
## @Date :   20200113       ##
##############################

## MultiQC report of Bisulphite-Seq data

all: task1
task1:
	multiqc --version > log && multiqc -d -f -i "MultiQC report for SC Bisulphite-Seq data (literature)</b>" -b "<i>--by</i> <br> <b>Deepak Tanwar</b>" -c multiqc_config.yaml -o ./output -z -v ./input  2>> log
