# META ANALYSIS tools
Tools for doing x way meta-analysis

# COVID-19 HGI meta-analysis workflow

This repository is used for COVID-19 Host Genetics Initiative GWAS meta-analysis. [WDL](https://github.com/openwdl/wdl) workflows and Google Compute Engine are used for computing. The workflows consist of cleaning/munging input files to the same format and running a meta-analysis.

1. The WDL workflow in [wdl/munge_sumstats.wdl](wdl/munge_sumstats.wdl) and [wdl/munge_sumstats.json](wdl/munge_sumstats.json) is used to filter and convert submitted SAIGE summary stat files to a unified format. INFO and AF filtering are done to the each summary stat file. Stats in build 37 are automatically lifted to build 38. Alleles are harmonized (matching ref/alt alleles, effect direction) using gnomAD 3.0 genomes as reference and fold change of AF to gnomAD AF for the population is added to the stats. Chromosomes are renamed so that e.g. "chr1" and "01" become "1" and "X" becomes "23". Scientific notation is converted to decimal notation for base pair positions. The output of this step is for each study a [bgzipped](http://www.htslib.org/doc/bgzip.html) tab-delimited summary stat file and its [tabix](http://www.htslib.org/doc/tabix.html) index, as well as manhattan and qq plots and AF-gnomAD_AF plots.
1.1. Before running the above workflow, scripts in scripts/format can be used to format non-SAIGE summary stats to SAIGE format.

2. The WDL workflow in [wdl/meta.wdl](wdl/meta.wdl), [wdl/meta.sub.wdl](wdl/meta.sub.wdl) and [wdl/meta.json](wdl/meta.json) is used to run meta-analysis with inverse-variance weighted betas. The analysis program is in [scripts/meta_analysis.py](scripts/meta_analysis.py). The output is a bgzipped, tab-delimited summary stat file with summary stats of each individual study and meta-analysis stats across all studies that have each variant. Leave-one-out meta-analysis can also be performed so that stats are available using all but one study in the meta-analysis. Manhattan and qq plots are also created.

## Docker image

A Docker image is available in [Docker Hub](https://hub.docker.com/repository/docker/covid19hg/meta)

To build a Docker image:

```
git clone https://github.com/covid19-hg/META_ANALYSIS
docker build -t covid19hg/meta:TAG_NAME -f META_ANALYSIS/docker/Dockerfile .
```
Replace TAG_NAME with the desired tag. The latest commit hash is used as tag name for versioning in Docker Hub.
