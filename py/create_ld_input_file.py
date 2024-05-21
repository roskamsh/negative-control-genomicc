import pandas as pd
import numpy as np
import requests

bqtls = pd.read_csv("output/snps/negative_control_bqtls.csv")
eqtls = pd.read_csv("output/snps/negative_control_transactors.csv")

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
final.to_csv("output/snps/snps_for_ld_block_removal.csv", index = False)