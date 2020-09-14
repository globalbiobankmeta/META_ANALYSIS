FROM gcr.io/finngen-refinery-dev/bioinformatics:0.6

RUN echo "hello world2"
ADD install_packages.R /usr/local/bin/
RUN Rscript /usr/local/bin/install_packages.R

RUN echo "update8"
ADD META_ANALYSIS /META_ANALYSIS

RUN python3 /META_ANALYSIS/scripts/harmonize_merge_gnomadAF.py --help
RUN python3 /META_ANALYSIS/scripts/harmonize_postGWASQC.py --help


ADD annovar /annovar
#RUN perl /annovar/annotate_variation.pl --help 


RUN Rscript /annovar/mergeAnno.R --help
