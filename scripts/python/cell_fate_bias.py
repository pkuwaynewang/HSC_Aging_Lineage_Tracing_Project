 import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp
import sys
import warnings
warnings.filterwarnings('ignore')

def clonal_fate_bias (text_file):
    df = pd.read_csv(text_file, sep='\t', header=0)
    df['all_cell_number'] = df['clone_size'].sum()
    df['fate_number'] = df['fate_size'].sum()
    df["P1"] = stats.hypergeom(M=df["all_cell_number"], n=df["fate_number"], N=df["clone_size"]).cdf(df["fate_size"])
    df["P2"] = stats.hypergeom(M=df["all_cell_number"], n=df["fate_number"], N=df["clone_size"]).sf(df["fate_size"])
    df["P_min"] = df[['P1', 'P2']].min(axis=1)
    def hypothesis_testing(pv):
        qv = statsmodels.sandbox.stats.multicomp.multipletests(
            list(pv), alpha=0.05, method="fdr_bh")[1]
        return qv
    
    df_list = []
    for x in df["clone_size"].unique():
        df_tmp = df[(df["clone_size"] == x)]
        if len(df_tmp) > 0:
            df_tmp["Q_value"] = hypothesis_testing(df_tmp["P_min"])
            df_list.append(df_tmp)
    df_final = pd.concat(df_list)
    
    df_final["bias"] = -np.log10(df_final["Q_value"])
    df_final['fate_bias'] = np.where(df_final['P2'] < df_final['P1'], 
                                 1 * df_final['bias'], 
                                -1 * df_final['bias'])
    output_file = text_file.split('.')[0] + '_fate_bias_results.txt'
    df_final.to_csv(output_file, sep='\t', index=False, header=True)
    #logg.info("Data saved in new files named [bias.files.txt]")  
    return output_file  
# Main block to run the function with a command line argument
if __name__ == "__main__":
    # Check if the script has been given an argument
    if len(sys.argv) != 2:
        print("Usage: python fate.bias.py file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    result_file = clonal_fate_bias(input_file)
    print(f"Your analysis result is saved in {result_file}") 
