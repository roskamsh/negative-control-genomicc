import pandas as pd

assembly = "hg38"
exclusion_file = pd.read_csv("data/genomicc_associated_SNPs_BGEN_ID_with_LD_blocks.csv")
exclusion_file["lower_bound"] = exclusion_file["lower_bound"].apply(int)
exclusion_file["upper_bound"] = exclusion_file["lower_bound"].apply(int)
exclusion_ranges =  exclusion_file['CHR'].astype(str) + ':' + exclusion_file['lower_bound'].astype(str) + '-' + exclusion_file['upper_bound'].astype(str)

exclusion_ranges_formatted = " ".join(exclusion_ranges)

# Add FlashPCA2 exclusion regions to this
if (assembly in ["hg19","grch37"]):
    exclusion_ranges_formatted = exclusion_ranges_formatted + " 5:44000000-51500000 6:25000000-33500000 8:8000000-12000000 11:45000000-57000000"
elif (assembly in ["hg38","grch38"]):
    exclusion_ranges_formatted = exclusion_ranges_formatted + " 5:43999898-52204166 6:24999772-33532223 8:8142478-12142491 11:44978449-57232526"
else:
    raise ValueError("No valid assembly given. Please choose either hg19/grch37 or hg38/grch38.")

with open('exclusion_regions.txt', 'w') as file:
    # Write the string to the file
    file.write(exclusion_ranges_formatted)