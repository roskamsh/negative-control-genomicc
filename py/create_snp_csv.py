import pandas as pd

def pull_chr(string):
    return string.split(":")[0]

def remove_chr_prefix(s):
    if s.startswith('chr'):
        return s[3:]
    else:
        return s

if __name__=="__main__":
    bqtls = pd.read_csv('data/final_bqtls_for_targene_one_snp_removed.csv')
    eqtls = pd.read_csv('data/final_eqtls_for_targene.csv')

    # label QTL
    bqtls["TYPE"] = "bQTL"
    eqtls["TYPE"] = "eQTL"
    all_snps = pd.concat([bqtls,eqtls],axis=0)

    # Extract all CHRs needed and make CSV file with each unique chromosome per row
    all_snps["CHR"] = all_snps["ID"].apply(pull_chr)
    all_snps["CHR"] = all_snps["CHR"].apply(remove_chr_prefix)

    all_snps.to_csv('snps.csv', index = False)