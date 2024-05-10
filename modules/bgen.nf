process create_snp_channel {
    container 'roskamsh/bgen_env:0.3.0'
    input:
        path BQTLS
        path EQTLS

    output:
        path "snps.csv"

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        def pull_chr(string):
            return string.split(":")[0]
        
        def remove_chr_prefix(s):
            if s.startswith('chr'):
                return s[3:]
            else:
                return s

        bqtls = pd.read_csv("${BQTLS}")
        eqtls = pd.read_csv("${EQTLS}")

        # Label type
        bqtls["TYPE"] = "bQTL"
        eqtls["TYPE"] = "eQTL"
        all_snps = pd.concat([bqtls,eqtls],axis=0)

        # Label CHRs in dataframe
        all_snps["CHR"] = all_snps["ID"].apply(pull_chr)
        all_snps["CHR"] = all_snps["CHR"].apply(remove_chr_prefix)

        all_snps.to_csv('snps.csv', index = False)
        """
}

process generate_info_score {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'
    publishDir("${params.OUTDIR}/info_scores", pattern: "*.snpstats") 

    input:
        tuple val(chr), val(prefix), path(files)

    output:
        tuple val(chr), path("chr${chr}.snpstats")

    script:
        """
        qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -snp-stats -osnp chr${chr}.snpstats
        """
}

process filter_low_quality_snps {
    container 'roskamsh/bgen_env:0.3.0'

    input:
        tuple val(chr), path(snpstats)

    output:
        path("chr${chr}_info_score0.6.txt")

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=9)
        df = df[df["info"] >= 0.6]

        snps = df[["rsid","chromosome","position","alleleA","alleleB","HW_exact_p_value",
                "minor_allele_frequency","minor_allele","major_allele","info"]]

        # Write table
        snps.to_csv("chr${chr}_info_score0.6.txt", index = False)
        """
}

process compile_info_score {
    container 'roskamsh/bgen_env:0.3.0'

    input:
        path files

    output:
        path 'info_score0.6_all_chrs.csv'

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    files = "${files}".strip('[]').split(' ')

    dfs = []
    for file in files:
        df = pd.read_csv(file)
        dfs.append(df)
    
    compiled = pd.concat(dfs, ignore_index=True)

    compiled.to_csv("info_score0.6_all_chrs.csv", index = False)
    """
}

workflow info_score {
    main:
        Channel
            .fromPath(params.BQTLS, checkIfExists: true)
            .set { bqtls }
        Channel
            .fromPath(params.EQTLS, checkIfExists: true)
            .set { eqtls }
        Channel
            .fromFilePairs(params.BGEN_FILES, size: -1, checkIfExists: true)
            .multiMap {
                prefix, files ->
                prefix: prefix
                files: [files]
            }
            .set { bgen_files_ch }

        snps_csv = create_snp_channel(bqtls, eqtls)
        chrs_ch = Channel.of(1..22)

        chrs_ch
            .combine(bgen_files_ch.prefix)
            .combine(bgen_files_ch.files)
            .set { chr_bgen_ch }
            
        generate_info_score(chr_bgen_ch)

        filter_low_quality_snps(generate_info_score.out)
    
        filter_low_quality_snps.out
            .collect()
            .set { all_info_scores }

        compile_info_score(all_info_scores) 

        // Configure input channel for SNPs
        snps_csv
            .splitCsv(header:true)
            .map { row -> tuple(row.CHR, row.ID, row.TYPE, row.TF) }
            .combine(compile_info_score.out)
            .set { snps_info_score_ch }
        
    emit:
        snps_info_score_ch
        
}