task run_range {

    String docker
    String pheno
    String method
    File conf
    # not used but needed to localize files - locally referenced to in the conf file
    Array[File] summary_stats
    String opts
    String chrom
    Int n_studies = length(summary_stats)

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - run meta"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"
        echo "method: ${method}"
        echo "options: ${opts}"
        echo "chrom: ${chrom}"
        echo "conf: ${conf}"
        echo "n studies: ${n_studies}"
        printf "${sep='\n' summary_stats}\n"

        /META_ANALYSIS/scripts/meta_analysis.py ${opts} --chrom ${chrom} ${conf} "${pheno}_${method}_chr${chrom}_meta" ${method} && \
        echo "`date` done"

    >>>

    output {
        File out = pheno + "_" + method + "_chr" + chrom + "_meta.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task gather {

    String docker
    Array[File] meta_stats
    String pheno
    String method
    File conf
    Array[String] summary_stats
    String opts
    Int n_studies = length(summary_stats)
    Int n_pieces = length(meta_stats)
    Int min_n_studies
    Int loglog_ylim

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - gather results"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"
        echo "method: ${method}"
        echo "options: ${opts}"
        echo "conf: ${conf}"
        echo "n result pieces: ${n_pieces}"
        echo "n studies: ${n_studies}"
        printf "${sep='\n' summary_stats}\n"

        echo "`date` gathering result pieces into one"
        cat <(gunzip -c ${meta_stats[0]} | head -1) \
            <(for file in ${sep=" " meta_stats}; do
                  gunzip -c $file | awk '
                  NR==1 {for (i=1;i<=NF;i++) a[$i]=i}
                  NR>1 && $a["all_meta_N"] >= ${min_n_studies}
                  ';
              done) \
        | bgzip > ${pheno}_${method}_meta.gz

        echo "`date` plotting qq and manhattan"
        gunzip -c ${pheno}_${method}_meta.gz | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) a[$i]=i;}
        {print $a["#CHR"],$a["POS"],$a["all_${method}_meta_p"]}
        ' > ${pheno}_${method}_meta_p
        qqplot.R --file ${pheno}_${method}_meta_p --bp_col "POS" --chrcol "#CHR" --pval_col "all_${method}_meta_p" --loglog_ylim ${loglog_ylim}

        echo "`date` done"

    >>>

    output {
        File out = pheno + "_" + method + "_meta.gz"
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task add_rsids_af {

    File file
    String base = basename(file, ".gz")
    File file_ref
    String method
    Float p_thresh
    String docker

    command <<<

        python3 <<EOF | bgzip > ${base}.gz

        import gzip, numpy

        fp_ref = gzip.open('${file_ref}', 'rt')
        ref_has_lines = True
        ref_chr = 1
        ref_pos = 0
        ref_h_idx = {h:i for i,h in enumerate(fp_ref.readline().strip().split('\t'))}

        with gzip.open('${file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            n_idx = []
            af_idx = []
            for i,h in enumerate(header.split('\t')):
                if h.endswith('_AF_Allele2'):
                    af_idx.append(i)
                if h.endswith('_N') and not h.endswith('_meta_N'):
                    n_idx.append(i)
            if len(n_idx) != len(af_idx):
                raise Exception('Unexpected header in ${file} (' + str(len(n_idx)) + ') _N fields, (' + str(len(af_idx)) + ') AF fields')
            print(header + '\tall_meta_sample_N\tall_meta_AF\trsid')
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = int(s[h_idx['#CHR']])
                pos = int(s[h_idx['POS']])
                ref = s[h_idx['REF']]
                alt = s[h_idx['ALT']]
                ref_vars = []
                while ref_has_lines and int(ref_chr) < chr or (int(ref_chr) == chr and ref_pos < pos):
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#chr']]
                        ref_pos = int(ref_line[ref_h_idx['pos']])
                    except IndexError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chr and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#chr']]
                        ref_pos = int(ref_line[ref_h_idx['pos']])
                    except IndexError:
                        ref_has_lines = False

                rsid = 'NA'
                for r in ref_vars:
                    if r[ref_h_idx['ref']] == ref and r[ref_h_idx['alt']] == alt:
                        rsid = r[ref_h_idx['rsid']]
                        break

                af_total = 0
                n_sum = 0
                for i,idx in enumerate(af_idx):
                    if s[idx] != 'NA':
                        af_total = af_total + float(s[idx]) * int(s[n_idx[i]])
                        n_sum = n_sum + int(s[n_idx[i]])
                af_total = af_total / n_sum

                print(line + '\t' + str(n_sum) + '\t' + numpy.format_float_scientific(af_total, precision=3) + '\t' + rsid)
        EOF

        echo "`date` filtering p-value ${p_thresh}"
        gunzip -c ${base}.gz | awk '
        NR==1 {for (i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["all_${method}_meta_p"] < ${p_thresh}
        ' > ${base}_${p_thresh}.txt

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ${base}.gz
        echo "`date` done"

    >>>

    output {
        File out = base + ".gz"
        File out_tbi = base + ".gz.tbi"
        File sign = base + "_" + p_thresh + ".txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task filter_cols {

    String docker
    File file
    String pheno
    String method
    Float p_thresh
    String out_template
    String outfile = sub(out_template, "\\{PHENO\\}", pheno)

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - filter columns"
        echo "docker: ${docker}"
        echo "file: ${file}"
        echo "pheno: ${pheno}"
        echo "out_template: ${out_template}"

        echo "`date` filtering columns"
        gunzip -c ${file} | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) { a[$i]=i; if (i<6||$i~"_AF_Allele2$"||$i~"_AF_fc$"||$i~"_N$"||$i~"^all_"||$i~"^rsid") use[i]=1 }}
        NR>=1 {printf $1; for(i=2;i<=NF;i++) if(use[i]==1) printf "\t"$i; printf "\n"}' | \
        bgzip > ${outfile}

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ${outfile}

        echo "`date` filtering p-value ${p_thresh}"
        gunzip -c ${outfile} | awk '
        NR==1 {for (i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["all_${method}_meta_p"] < ${p_thresh}
        ' > ${outfile}_${p_thresh}.txt

        echo "`date` done"

    >>>

    output {
        File out = outfile
        File out_tbi = outfile + ".tbi"
        File sign = outfile + "_" + p_thresh + ".txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task lift {

    String docker
    File file
    String base = basename(file, ".txt.gz")
    String method
    Float p_thresh

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - liftover to 37"
        echo "docker: ${docker}"
        echo "file: ${file}"
        echo "method: ${method}"
        echo "p_thresh: ${p_thresh}"

        TEMP=/cromwell_root
        TMP=/cromwell_root
        TMPDIR=/cromwell_root

        mv ${file} ${base}.gz

        echo "`date` liftover"
        time /META_ANALYSIS/scripts/lift.py -chr "#CHR" -pos POS -ref REF -alt ALT \
        -chain_file /liftover/hg38ToHg19.over.chain.gz -tmp_path /cromwell_root/ \
        ${base}.gz > ${base}.lift.out 2> ${base}.lift.err

        gunzip -c ${base}.gz.lifted.gz | \
        cut -f2- | awk '
        BEGIN { FS=OFS="\t" }
        NR==1 { for (i=1;i<=NF;i++) a[$i]=i; print $0 }
        NR>1 {
            temp=$a["#CHR"]; $a["#CHR"]=$a["anew_chr"]; $a["anew_chr"]=temp; temp=$a["POS"]; $a["POS"]=$a["anew_pos"]; $a["anew_pos"]=temp;
            sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]);
            if ($a["#CHR"] ~ /^[0-9]+$/) {
                $a["SNP"] = $a["#CHR"]":"$a["POS"]":"$a["REF"]":"$a["ALT"]
                print $0
            }
        }' | bgzip > ${base}.b37.txt.gz

        echo "`date` tabixing"
        tabix -s1 -b2 -e2 ${base}.b37.txt.gz

        echo "`date` filtering p-value ${p_thresh}"
        gunzip -c ${base}.b37.txt.gz | awk '
        NR==1 {for (i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["all_${method}_meta_p"] < ${p_thresh}
        ' > ${base}.b37_${p_thresh}.txt

        echo "`date` done"

    >>>

    output {
        File lift_out = base + ".lift.out"
        File lift_err = base + ".lift.err"
        File out = base + ".b37.txt.gz"
        File out_tbi = base + ".b37.txt.gz.tbi"
        File sign = base + ".b37_" + p_thresh + ".txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow run_meta {

    String pheno
    Int min_n_studies
    String conf
    String method
    String opts
    Array[String] summary_stats

    scatter (chr in range(23)) {
        call run_range {
            input: pheno=pheno, method=method, opts=opts, conf=conf, summary_stats=summary_stats, chrom=chr+1
        }
    }

    call gather {
        input: pheno=pheno, method=method, opts=opts, conf=conf, summary_stats=summary_stats, meta_stats=run_range.out, min_n_studies=min_n_studies
    }

    call add_rsids_af {
        input: file=gather.out, method=method
    }

    call filter_cols {
        input: pheno=pheno, file=add_rsids_af.out, method=method
    }

    call lift {
        input: file=filter_cols.out, method=method
    }
}
