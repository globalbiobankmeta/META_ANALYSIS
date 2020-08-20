task plot {

    File stats
    String base = basename(stats, ".gz")
    String col
    Int loglog_ylim
    String docker

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - plot"
        echo "${stats}"
        echo "${col}"
        gunzip -c ${stats} | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) a[$i]=i;}
        {print $a["#CHR"],$a["POS"],$a["${col}"]}' \
        > ${col}
        qqplot.R --file ${col} --bp_col "POS" --chrcol "#CHR" --pval_col ${col} --loglog_ylim ${loglog_ylim}
        echo "`date` done"

    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow qq_manha {

    File colfile
    Array[String] cols = read_lines(colfile)
    String stats

    scatter (col in cols) {
        call plot {
            input: stats=stats, col=col
        }
    }
}
