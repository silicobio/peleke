import pandas as pd
from abnumber import Chain
import os

## Load antigens and generated sequences from Excel file
test_antigens_df = pd.read_excel('test_cases.xlsx', sheet_name='test_cases')
generated_seqs_df = pd.read_excel('test_cases.xlsx', sheet_name='generated_seqs')

## Join dataframes on test_antigens_df
test_cases_df = pd.merge(generated_seqs_df, test_antigens_df, on='antigen_id', suffixes=('_gen', '_test'))
test_cases_df.drop(columns=['model', 'pdb_id', 'source', 'antigen_name', 'antigen_source', 'antigen_ids', 'highlighted_epitope_seqs'], inplace=True)

## Define numbering test
def test_numbering(seq:str="", chain_type:str="") -> str:
    try:
        chain = Chain(seq, scheme="chothia")
        if chain_type == "heavy":
            if not chain.is_heavy_chain():
                return 'FAIL'
            else:
                return 'PASS'
        elif chain_type == "light":
            if not chain.is_light_chain():
                return 'FAIL'
            else:
                return 'PASS'
        else:
            return 'PASS'
    except Exception as e:
        print(f"Error with sequence {seq}: {e}")
        return 'FAIL'

## Run tests and store results
i = 1
for idx, row in test_cases_df.iterrows():
    print(f"Running test case {i} of {len(test_cases_df)}")
    h_chain_seq = row['h_chain']
    l_chain_seq = row['l_chain']
    ## Test numbering
    test_cases_df.at[idx, 'test_h_chain_numbering'] = test_numbering(seq=h_chain_seq, chain_type="heavy")
    test_cases_df.at[idx, 'test_l_chain_numbering'] = test_numbering(seq=l_chain_seq, chain_type="light")
    i += 1

## Save results to new Excel file
output_file = 'numbering_test_results.csv'
test_cases_df.to_csv(output_file, index=False)