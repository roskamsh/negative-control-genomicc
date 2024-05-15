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
eqtls.to_csv("negative_control_transactors.csv")
bqtls.to_csv("negative_control_bqtls.csv")