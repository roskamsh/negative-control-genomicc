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

        all_snps = all_snps.drop_duplicates()
        all_snps.to_csv('snps.csv', index = False)
        """
}

process define_exclusion_regions {
    container 'roskamsh/bgen_env:0.3.0'

    input:
        path exclusion_regions

    output:
        path "exclusion_regions.txt"

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        exclusion_file = pd.read_csv("${exclusion_regions}")
        exclusion_file["lower_bound"] = exclusion_file["lower_bound"].apply(int)
        exclusion_file["upper_bound"] = exclusion_file["lower_bound"].apply(int)
        exclusion_ranges =  exclusion_file['CHR'].astype(str) + ':' + exclusion_file['lower_bound'].astype(str) + '-' + exclusion_file['upper_bound'].astype(str)

        # Format like QCtool requires, space-separated
        exclusion_ranges_formatted = " ".join(exclusion_ranges)

        with open('exclusion_regions.txt', 'w') as file:
            # Write the string to the file
            file.write(exclusion_ranges_formatted) 
        """
}

process generate_info_score {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'
    publishDir("${params.OUTDIR}/info_scores", pattern: "*.snpstats", mode: "copy") 

    input:
        tuple val(chr), val(prefix), path(files), path(exclusion_regions)

    output:
        tuple val(chr), path("chr${chr}.snpstats")

    script:
        """
        ranges=\$( cat ${exclusion_regions} )
        qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -excl-range \$ranges -snp-stats -osnp chr${chr}.snpstats
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

        df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=10)
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
        Channel
            .fromPath(params.EXCLUSION_REGIONS, checkIfExists: true)
            .set { exclusion_regions }

        snps_csv = create_snp_channel(bqtls, eqtls)
        chrs_ch = Channel.of(1..22)

        define_exclusion_regions(exclusion_regions)

        chrs_ch
            .combine(bgen_files_ch.prefix)
            .combine(bgen_files_ch.files)
            .combine(define_exclusion_regions.out)
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