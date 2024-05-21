process get_maf_match {
    container 'roskamsh/bgen_env:0.3.0'

    input:
    tuple val(CHR), val(SNP), val(TYPE), val(TF), path(INFO_SCORES), val(REPLICATE)

    output:
    path "${SNP}_${TYPE}_${TF}_rep${REPLICATE}_negative_control_variant.csv"

    script:
    """
    #!/usr/bin/env python

    import numpy as np
    import pandas as pd
    import requests
    import time

    def parse_snp_format(snp_id):
        parts = snp_id.split(':')
        if len(parts) == 4:
            chromosome, position, ref, alt = parts
            if chromosome.startswith("chr"):
                chromosome = chromosome[3:]
            return chromosome, int(position), ref, alt
        else:
            raise ValueError("Invalid SNP format. Please use 'chr:pos:ref:alt'.")

    def check_snp_in_functional_region(snp_id, assembly_version='GRCh38'):
        chromosome, position, _, _ = parse_snp_format(snp_id)
        print(f"Checking if SNP {snp_id} lies in functional region using ENSEMBL...")
        # ENSEMBL API endpoint for fetching information about a SNP with specified assembly version
        ensembl_api_url = f"https://rest.ensembl.org/overlap/region/human/{chromosome}:{position}-{position}?feature=regulatory;feature=exon;feature=other_regulatory"
        # Wrapper for retry strategy in case connection is refused because of too many requests
        num_retries = 5
        for _ in range(num_retries):
            try:
                response = requests.get(ensembl_api_url, headers={"Content-Type": "application/json"})
                if response.status_code in [200, 404]:
                    ## Escape for loop if returns a successful response
                    break
            except requests.exceptions.ConnectionError:
                print("Connection refused by the server..")
                print("Let me sleep for 5 seconds")
                print("ZZzzzz...")
                time.sleep(5) 
                pass
        if response.status_code == 200:
            snp_info = response.json()
            if len(snp_info) > 0:
                return True
            else:
                # If SNP is not in a functional region
                return False
        else:
            # Error handling if request fails
            print(f"Failed to fetch SNP information. Status code: {response.status_code}")
            return False

    snp = "${SNP}"
    qtltype = "${TYPE}"
    tf = "${TF}"
    replicate = int("${REPLICATE}")
    all_snps = pd.read_csv("${INFO_SCORES}")
    relative_threshold = 0.01
    rng = 123*replicate 

    all_snps_maf = all_snps.minor_allele_frequency.values

    print(f"Finding candidate MAF-matched SNPs for {snp}.")
    variant_maf = all_snps[all_snps.rsid == snp].minor_allele_frequency.values[0]
    match_maf = abs((variant_maf - all_snps_maf)/variant_maf) <= relative_threshold
    variants_that_match = all_snps[match_maf].copy()
    ncandidate = variants_that_match.shape[0]
    print(f"{ncandidate} MAF-matched variants found from initial search. Will sample from this list and then check for functional annotation.")
    is_functional = True
    while is_functional:
        sampled_variant = variants_that_match['rsid'].sample(n=1, random_state=rng).values[0]
        is_functional = check_snp_in_functional_region(sampled_variant)
        rng = rng + 1
    print(f"Found MAF-matched SNP {sampled_variant}.")
    out = pd.DataFrame([[snp,replicate,sampled_variant,variant_maf,qtltype,tf]], columns = ['input_snp','replicate','negative_control_snp','maf_to_match','type','tf'])
    out.to_csv(f"{snp}_{qtltype}_{tf}_rep{replicate}_negative_control_variant.csv", index = False)
    """
}

process compile_results {
    container 'roskamsh/bgen_env:0.3.0'
    publishDir("${params.OUTDIR}/snps", pattern: "*.csv", mode: "copy")

    input:
    path files

    output:
    path "negative_control_bqtls.csv", emit: bqtls
    path "negative_control_transactors.csv", emit: transactors

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import glob

    files = glob.glob("./*_negative_control_variant.csv")
    dfs = []

    for file in files:
        df = pd.read_csv(file)
        dfs.append(df)

    final = pd.concat(dfs, ignore_index = True)
    eqtls = final[final.type == 'eQTL'].copy()
    bqtls = final[final.type == 'bQTL'].copy()
    eqtls.to_csv("negative_control_transactors.csv", index = False)
    bqtls.to_csv("negative_control_bqtls.csv", index = False)
    """

}

process create_ld_file {
    container 'roskamsh/bgen_env:0.3.0'
    publishDir("${params.OUTDIR}/snps", mode: "copy")

    input:
    path bqtls
    path transactors

    output:
    path "snps_for_ld_block_removal.csv"

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    bqtls = pd.read_csv("${bqtls}")
    eqtls = pd.read_csv("${transactors}")

    def get_chr_pos(id):
        chr, pos, _, _ = id.split(":")
        if chr.startswith("chr"):
            chr = chr[3:]
        return chr, pos

    out = pd.concat([bqtls,eqtls])
    out[["CHR","POS"]] = out["negative_control_snp"].apply(lambda x: pd.Series(get_chr_pos(x)))

    out = out[["negative_control_snp","CHR","POS"]]
    out.columns = ["RSID","CHR","POS"]
    final = out.drop_duplicates()
    final.to_csv("snps_for_ld_block_removal.csv", index = False)
    """

}

workflow generate_negative_control {
    take:
    snps_info_score

    main:
    nreplicates = Channel.from(1..params.NREPLICATES)

    // Create combined channel which is the cartesian product of snps_info_score and nreplicates
    snps_info_score.combine(nreplicates)
        .set { snps_info_replicate }
    
    get_maf_match(snps_info_replicate)

    get_maf_match.out
        .collect()
        .set { negative_control_files }

    compile_results(negative_control_files)

    create_ld_file(compile_results.out.bqtls, compile_results.out.transactors)
}
