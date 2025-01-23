
#!/usr/bin/env python3
# Harmonizing GWAS summary statistics for meta-analyses 

import argparse
import pandas as pd
from multiprocessing import Pool
import os
import logging

# Argument parsing
parser = argparse.ArgumentParser(description="Harmonize population data.")
parser.add_argument("--pop_a_file", required=True, help="File path for Population A data (TSV format).")
parser.add_argument("--pop_b_file", required=True, help="File path for Population B data (TSV format).")
parser.add_argument("--output_dir", required=True, help="Directory to save harmonized output files.")
parser.add_argument("--log_dir", required=False, default=os.getcwd(), help="Directory to save the log file (default: current directory).")
parser.add_argument("--beta_threshold", required=False, default=100, type=float, help="Threshold for extreme BETA values.")
parser.add_argument("--snp_col", required=False, default="SNP", help="Name of the SNP column (default: 'SNP').")
parser.add_argument("--chr_col", required=False, default="CHR", help="Name of the chromosome column (default: 'CHR').")
parser.add_argument("--pos_col", required=False, default="POS", help="Name of the position column (default: 'POS').")
parser.add_argument("--a1_col", required=False, default="A1", help="Name of the allele 1 column (default: 'A1').")
parser.add_argument("--a2_col", required=False, default="A2", help="Name of the allele 2 column (default: 'A2').")
parser.add_argument("--beta_col", required=False, default="BETA", help="Name of the BETA column (default: 'BETA').")
parser.add_argument("--se_col", required=False, default="SE", help="Name of the SE column (default: 'SE').")
parser.add_argument("--p_col", required=False, default="P", help="Name of the P column (default: 'P').")
parser.add_argument("--eaf_col", required=False, default="EAF", help="Name of the EAF column (default: 'EAF').")
parser.add_argument("--chr_pos_col", required=False, default="CHR_POS", help="Name of the CHR_POS column (default: 'CHR_POS').")
args = parser.parse_args()

# Create the log file path in the specified directory
log_file_path = os.path.join(args.log_dir, "harmonization.log")

# Set up logging
logging.basicConfig(
    filename=log_file_path,
    filemode='w',  # Overwrite the log file each run
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logging.getLogger().addHandler(console_handler)
logging.info("Starting harmonization process")

# Load population data from .tsv files with error handling
try:
    pop_A = pd.read_csv(args.pop_a_file, sep='\t', low_memory=False)
    pop_B = pd.read_csv(args.pop_b_file, sep='\t', low_memory=False)
    logging.info("Successfully loaded population data files.")
except Exception as e:
    logging.error(f"Error reading input files: {e}")
    exit(1)

# Check for required columns
required_columns = [args.snp_col, args.chr_col, args.pos_col, args.a1_col, args.a2_col, args.beta_col, args.se_col, args.p_col, args.eaf_col, args.chr_pos_col]
for df, pop_name in zip([pop_A, pop_B], ['A', 'B']):
    for col in required_columns:
        if col not in df.columns:
            logging.error(f"Missing column '{col}' in {pop_name} population data.")
            raise ValueError(f"Missing column '{col}' in {pop_name} population data.")

# Check for duplicate SNPs and remove them
for pop, name in zip([pop_A, pop_B], ['A', 'B']):
    if pop.duplicated(subset=[args.snp_col]).any():
        logging.warning(f"Duplicates detected in population {name}. Removing duplicates.")
        pop.drop_duplicates(subset=[args.snp_col], inplace=True)

# Check for extreme BETA values (outliers)
for pop, name in zip([pop_A, pop_B], ['A', 'B']):
    if (pop[args.beta_col].abs() > args.beta_threshold).any():
        logging.warning(f"Extreme BETA values detected in population {name}. Proceed with caution.")

# Check for missing chromosome and position values
for pop, name in zip([pop_A, pop_B], ['A', 'B']):
    if pop[[args.chr_col, args.pos_col]].isnull().any().any():
        logging.warning(f"Missing chromosome or position values in population {name}.")

# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
    logging.info(f"Created output directory: {args.output_dir}")


# Count the number of SNPs present in both pop_A and pop_B
common_snps = pop_A[pop_A[args.snp_col].isin(pop_B[args.snp_col])]
num_common_snps = common_snps.shape[0]
num_snps_popA= pop_A.shape[0]
num_snps_popB= pop_B.shape[0]

logging.info(f"Number of SNPs present in both Population A and Population B: {num_common_snps}")
logging.info(f"Number of SNPs present in Population A: {num_snps_popA}")
logging.info(f"Number of SNPs present in Population B: {num_snps_popB}")


def are_complementary(allele1, allele2):
    """
    Check if two alleles are complementary based on DNA base-pairing rules.
    
    Args:
    allele1 (str): The first allele.
    allele2 (str): The second allele.
    
    Returns:
    bool: True if the alleles are complementary (A <-> T, C <-> G), False otherwise.
    """
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complements.get(allele1) == allele2


def harmonize_for_chromosome(pop1, pop2, snp_col, chr_col, pos_col, a1_col, a2_col, beta_col, chrom):
    """
    Harmonize allele data between two populations for a specific chromosome.
    
    Args:
    pop1 (pd.DataFrame): Population 1 data.
    pop2 (pd.DataFrame): Population 2 data to be harmonized.
    snp_col (str): Column name for SNP IDs.
    chr_col (str): Column name for chromosome.
    pos_col (str): Column name for position.
    a1_col (str): Column name for allele 1.
    a2_col (str): Column name for allele 2.
    beta_col (str): Column name for BETA values.
    chrom (str): Chromosome to harmonize.
    
    Returns:
    pd.DataFrame: Harmonized population 2 data.
    pd.DataFrame: SNPs that failed to harmonize.
    pd.DataFrame: SNPs with mismatched chromosome or position.
    int: Count of successfully harmonized SNPs.
    pd.DataFrame: Details of successful harmonizations.
    """
    try:
        pop1_chrom = pop1[pop1[chr_col] == chrom]
        pop2_chrom = pop2[pop2[chr_col] == chrom]
        
        pop2_harmonized = pop2_chrom.copy()
        failed_to_harmonize = []
        failed_chr_pos_match = []
        successfully_harmonized = []
        harmonized_count = 0

        for i, row in pop2_harmonized.iterrows():
            if row[snp_col] in pop1_chrom[snp_col].values:
                row1 = pop1_chrom[pop1_chrom[snp_col] == row[snp_col]].iloc[0]
                
                if row[chr_col] != row1[chr_col] or row[pos_col] != row1[pos_col]:
                    failed_chr_pos_match.append(row)
                    continue
                
                if row[a1_col] == row1[a1_col] and row[a2_col] == row1[a2_col]:
                    harmonized_count += 1
                    successfully_harmonized.append((row[snp_col], "Direct Match"))
                elif row[a1_col] == row1[a2_col] and row[a2_col] == row1[a1_col]:
                    pop2_harmonized.at[i, a1_col], pop2_harmonized.at[i, a2_col] = row[a2_col], row[a1_col]
                    pop2_harmonized.at[i, beta_col] = -row[beta_col]
                    harmonized_count += 1
                    successfully_harmonized.append((row[snp_col], "Flipped Alleles"))
                elif (are_complementary(row[a1_col], row1[a1_col]) and are_complementary(row[a2_col], row1[a2_col])) or \
                     (are_complementary(row[a1_col], row1[a2_col]) and are_complementary(row[a2_col], row1[a1_col])):
                    pop2_harmonized.at[i, a1_col] = row[a2_col]
                    pop2_harmonized.at[i, a2_col] = row[a1_col]
                    pop2_harmonized.at[i, beta_col] = -row[beta_col]
                    harmonized_count += 1
                    successfully_harmonized.append((row[snp_col], "Complementary Alleles"))
                else:
                    failed_to_harmonize.append(row)
        
        return pop2_harmonized, pd.DataFrame(failed_to_harmonize), pd.DataFrame(failed_chr_pos_match), harmonized_count, successfully_harmonized
    except Exception as e:
        logging.error(f"Error in harmonizing chromosome {chrom}: {e}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), 0, []


def harmonize_data(pop_A, pop_B, snp_col, chr_col, pos_col, a1_col, a2_col, beta_col, num_processes=4):
    """
    Harmonize allele data between two populations across all chromosomes using parallel processing.
    
    Args:
    pop_A (pd.DataFrame): Population A data.
    pop_B (pd.DataFrame): Population B data to be harmonized.
    snp_col (str): Column name for SNP IDs.
    chr_col (str): Column name for chromosome.
    pos_col (str): Column name for position.
    a1_col (str): Column name for allele 1.
    a2_col (str): Column name for allele 2.
    beta_col (str): Column name for BETA values.
    num_processes (int, optional): Number of processes for parallel execution. Default is 4.
    
    Returns:
    pd.DataFrame: Harmonized population B data.
    pd.DataFrame: SNPs that failed to harmonize.
    pd.DataFrame: SNPs with mismatched chromosome or position.
    int: Count of successfully harmonized SNPs.
    """
    chromosomes = pop_A[chr_col].unique()
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(
            harmonize_for_chromosome,
            [(pop_A, pop_B, snp_col, chr_col, pos_col, a1_col, a2_col, beta_col, chrom) for chrom in chromosomes]
        )
    
    # Combine results
    harmonized_data = pd.concat([result[0] for result in results], ignore_index=True)
    failed_to_harmonize = pd.concat([result[1] for result in results], ignore_index=True)
    failed_chr_pos_match = pd.concat([result[2] for result in results], ignore_index=True)
    harmonized_count = sum(result[3] for result in results)
    successfully_harmonized = [item for sublist in [result[4] for result in results] for item in sublist]
    
    logging.info(f"Successfully harmonized {harmonized_count} SNPs across all chromosomes.")
    return harmonized_data, failed_to_harmonize, failed_chr_pos_match, harmonized_count

# Harmonize data and save results
harmonized_data, failed_to_harmonize, failed_chr_pos_match, harmonized_count = harmonize_data(
    pop_A, pop_B, args.snp_col, args.chr_col, args.pos_col, args.a1_col, args.a2_col, args.beta_col
)

# Save the harmonized data and failed SNPs
harmonized_data.to_csv(os.path.join(args.output_dir, "harmonized_data.tsv"), sep='\t', index=False)
failed_to_harmonize.to_csv(os.path.join(args.output_dir, "failed_to_harmonize.tsv"), sep='\t', index=False)
failed_chr_pos_match.to_csv(os.path.join(args.output_dir, "failed_chr_pos_match.tsv"), sep='\t', index=False)

logging.info("Harmonization process complete.")
