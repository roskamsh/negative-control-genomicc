import pandas as pd
import requests

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
    
    # ENSEMBL API endpoint for fetching information about a SNP with specified assembly version
    ensembl_api_url = f"https://rest.ensembl.org/overlap/region/human/{chromosome}:{position}-{position}?feature=regulatory;feature=exon;feature=other_regulatory"

    # Make a request to the ENSEMBL API
    response = requests.get(ensembl_api_url, headers={"Content-Type": "application/json"})

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

df = pd.read_csv("chr21.snpstats", delimiter="\t", skiprows=9)
df = df[df["info"] >= 0.6]

snps = df[["rsid","chromosome","position","alleleA","alleleB","HW_exact_p_value",
           "minor_allele_frequency","minor_allele","major_allele","info"]]

## Query ENSEMBL and check if in regulatory region
snps['is_in_functional_region'] = snps["rsid"].apply(check_snp_in_functional_region)

# Write rsids
snps.to_csv("chr21_info_score0.6.txt", index = False)