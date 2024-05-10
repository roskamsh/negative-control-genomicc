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

    # Wrapper for retry startegy in case connection is refused
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

snps = ['chr9:133438003:C:G', 'chr9:133438003:C:G', 'chr9:21415688:G:A', 'chr9:21415688:G:A', 
        'chr9:21415688:G:A', 'chr9:133257570:C:T', 'chr9:133438003:C:G', 'chr9:133438003:C:G']
all_snps = pd.read_csv("chr9_info_score0.6.txt")
relative_threshold = 0.01
rng = 123

all_snps_maf = all_snps.minor_allele_frequency.values

for snp in snps:
    print(f"Finding candidate MAF-matched SNPs for {snp}.")
    variant_maf = all_snps[all_snps.rsid == snp].minor_allele_frequency.values[0]
    match_maf = abs((variant_maf - all_snps_maf)/variant_maf) <= relative_threshold
    variants_that_match = all_snps[match_maf].copy()
    ncandidate = variants_that_match.shape[0]
    print(f"{ncandidate} MAF-matched variants found from initial search. Now checking for functional annotation.")
    # Now - query ensembl
    variants_that_match['is_in_functional_region'] = variants_that_match["rsid"].apply(check_snp_in_functional_region)
    to_sample_from = variants_that_match[np.logical_not(variants_that_match.is_in_functional_region)]
    sampled_variant = to_sample_from['rsid'].sample(n=1, random_state=rng).values[0]
    out = pd.DataFrame([[snp,sampled_variant,variant_maf]], columns = ['input_snp','negative_control_snp','maf_to_match'])
    out.to_csv(f"{snp}_negative_control_variant.csv", index = False)
