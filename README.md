# META ANALYSIS tools
Tools for doing x way meta-analysis

# COVID-19 HGI meta-analysis workflow

This repository is used for COVID-19 Host Genetics Initiative GWAS meta-analysis. [WDL](https://github.com/openwdl/wdl) workflows and Google Compute Engine are used for computing. The workflows consist of cleaning/munging input files to the same format and running a meta-analysis.

1. The WDL workflow in [wdl/clean_sumstats.wdl](wdl/clean_sumstats.wdl) and [wdl/clean_sumstats.json](wdl/clean_sumstats.json) is used to convert submitted SAIGE summary stat files to tab-delimited files sorted by chromosome position. Chromosomes are renamed so that e.g. "chr1" and "01" become "1" and "X" becomes "23". Scientific notation is converted to decimal notation for base pair positions. The output of this step is for each study a [bgzipped](http://www.htslib.org/doc/bgzip.html) tab-delimited summary stat file and its [tabix](http://www.htslib.org/doc/tabix.html) index, as well as manhattan and qq plots.
1.1. Before running the above workflow, scripts in scripts/format can be used to format non-SAIGE summary stats to SAIGE format.

2. The WDL workflow in [wdl/meta.wdl](wdl/meta.wdl), [wdl/meta.sub.wdl](wdl/meta.sub.wdl) and [wdl/meta.json](wdl/meta.json) is used to run meta-analysis with inverse-variance weighted betas. The analysis program is in [scripts/meta_analysis.py](scripts/meta_analysis.py). The output is a bgzipped, tab-delimited summary stat file with summary stats of each individual study and meta-analysis stats across all studies that have the variant. Leave-one-out meta-analysis can also be performed so that stats are available using all but one study in the meta-analysis. Manhattan and qq plots are also created.

## Docker image

A Docker image is available in [Docker Hub](https://hub.docker.com/repository/docker/covid19hg/meta)

To build a Docker image:

```
git clone https://github.com/covid19-hg/META_ANALYSIS
docker build -t covid19hg/meta:TAG_NAME -f META_ANALYSIS/docker/Dockerfile .
```
Replace TAG_NAME with the desired tag. The latest commit hash is used as tag name for versioning in Docker Hub.

# General info

## Variant matching across studies
Variants are matched using chr pos ref and alt. For this reason results need to first be lifted over to the same build.
scripts/lift.py can be used to liftover results first if needed. Strands and flips are taken into account when matching variants.

IMPORTANT: Studies need to be ordered by chr (1-22, x,y,mt) and position. Chromosome can be indicated with numbers 1-25 or chr1-22, chrX,chrY,chrMT and they will be internally coded to numerical values.

## Running single trait meta-analysis
scripts/meta_analysis.py is the main script for running meta-analysis for a single trait. Meta-analyses to be performed are specified with json
configuration file. Example configuration in data/conf.json. Script will try to align using both strands as well as by flipping ref vs. alt.

First parameter should be a path to a json configuration file with these elements:
```
            "name":"STUDY_NAME",
            "file":"/path/to/sumstats.gz",
            "n_cases": 6570 , # number of cases. Used only if sample size weighted meta-analysis is used
            "n_controls": 48378, # number of controls. Used only if sample size weighted meta-analysis is used
            "chr":"CHR", #chromosome column name in the file
            "pos":"POS", #position column name in the file
            "ref":"Allele1", #reference allele column name in the file
            "alt":"Allele2", #alternate allele column name in the file (effect is for this allele)
            "effect":"BETA", #effect size column name in the file
            "effect_type":"beta", #is the effect column beta or or. In case of OR the value will be log transformed to beta.
            "pval":"p.value" # effect size column name in the file
            "se":"SE" <- this parameter is optional. If given for compared studies additional p-value will be added using this as a weight for z-score.
```
Second parameter is output prefix where results are written. 

scripts/meta_analysis.py supports 3 different meta-analysis methods N: purely weight by sample size and use z-score from p-value,
variance: weight z-score from p-value by variance, inv_var: regulare inverse variance weighted betas meta-analysis.
inv_var is recommeneded if betas and variances are comparable. In case of combining data from different models (e.g.) linear vs. logistic you should use sample size weighted meta.
