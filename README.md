# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h4 align="right">Trey Pridgen, Nicholas Santolla, Prbhuv Nigam, and Colby T. Ford</h4><h3 align="right">Silico Biosciences</h3>


## Generate Antibody Sequences

```bash
python scripts/generate.py --antigen "MKT[LLI]LAV[AA]A..." --model "peleke-phi-4"
```

## Training Dataset

The training dataset consists of paired antigen and antibody sequences, where the antigen is the target for which the antibody is generated. This was curated from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab). Using PandaProt, epitope residues were highlighted in the antibody sequences using `[ ]`, which were used to tune the model to generate antibodies sequences that fold and bind to specific epitopes on the desired antigen. Note that mutli-chain antigen sequences are delimited by `|` in the `antigen_ids` column, and the heavy and light chain antibody sequences are delimited by `|` in the `antibody_sequences` column.


<!-- | **pdb_id** | **h_chain_id** | **l_chain_id** | **antigen_ids** | **h_chain_seq** | **l_chain_seq** | **antigen_seqs** |    **antibody_sequences**    | **highlighted_epitope_sequences** |          **epitope_residues**         |
|:----------:|:--------------:|:--------------:|:---------------:|:---------------:|:---------------:|:----------------:|:----------------------------:|:---------------------------------:|:-------------------------------------:|
| 8xa4       | C              | D              | A\|B            | QLQLQESGPG...   | EIVLTQSPGT...   | SCNGLYYQGSC...   | QLQLQESGPG...\|EIVLTQSPGT... | ...LTTWLI[D][Y]V[E][D][T]WGS...   | A:ARG 176\|A:ASP 146\|A:ASP 150...    |
| 9cph       | H              | L              | A               | EVQLVESGGG...   | AQMTQSPSSL...   | KIEEGKLVIWI...   | EVQLVESGGG...\|AQMTQSPSSL... | ...FQDKL[Y][P][F]TW[D][A]VRY...   | A:ALA 1116\|A:ALA 1122\|A:ALA 1128... |
| 9d7i       | H              | G              | E               | VQLQESGPGV...   | YELTQPPSVS...   | LWVTVYYGVPV...   | VQLQESGPGV...\|YELTQPPSVS... | ...DVVQI[N]NKEYRL[I]NC[N][T]...   | E:ARG 429\|E:ARG 469\|E:ASN 177...    |
| 9d7i       | J              | I              | C               | VQLQESGPGV...   | YELTQPPSVS...   | LWVTVYYGVPV...   | VQLQESGPGV...\|YELTQPPSVS... | ...YRL[I]NC[N][T]SACTQACPKVS...   | C:ARG 469\|C:ASN 197\|C:ASN 280...    |
| 9d7o       | H              | G              | E               | QVQLQESGPG...   | YELTQPPSVS...   | LWVTVYYGVPV...   | QVQLQESGPG...\|YELTQPPSVS... | ...VQIN[E]SNKEYRL[I]NC[N]TSA...   | E:ARG 429\|E:ARG 469\|E:ASN 197...    | -->


|      **Column Name**     |                             **Description**                            |            **Example**           |
|:------------------------:|:----------------------------------------------------------------------:|:--------------------------------:|
| pdb_id                   | The PDB ID on Protein Data Bank                                        | 8xa4                             |
| h_chain_id               | The chain ID of the antibody's heavy chain                             | C                                |
| l_chain_id               | The chain ID of the antibody's light chain                             | D                                |
| antigen_ids              | A \|-delimited list of chain IDs of the antigen chain(s)               | A\|B                             |
| h_chain_seq              | The heavy chain amino acid sequence                                    | QLQLQESGPG…                      |
| l_chain_seq              | The light chain amino acid sequence                                    | EIVLTQSPGT…                      |
| antigen_seqs             | The antigen sequence(s), \|-delimited                                  | SCNGL...\|SCNGL…                 |
| antibody_seqs            | The heavy and light chain sequences, \|-delimited                      | QLQLQ...\|EIVLT…                 |
| h_chain_fv_seq           | The heavy chain sequence, trimmed to the Fv portion                    | QLQLQESGPG…                      |
| l_chain_fv_seq           | The light chain sequence, trimmed to the Fv portion                    | EIVLTQSPGT…                      |
| antibody_fv_seqs         | The heavy and light chain Fv sequences, \|-delimited                   | QLQLQ...\|EIVLT…                 |
| highlighted_epitope_seqs | The antigen sequence(s) with epitope residues encased in [ ]           | ...WLI[D][Y]V[E][D][T]WGS…       |
| epitope_residues         | The list of epitope residues, \|-delimited in a "chain: AA #"   format | A:ARG 176\|A:ASP 146\|A:ASP 150… |

- See the prepared training dataset: [data/sabdab/sabdab_training_dataset.csv](data/sabdab/sabdab_training_dataset.csv)
- Our data preparation scripts:
    1. Get sequences from PDB structures: [data/sabdab/01_get_structure_seqs.ipynb](data/sabdab/01_get_structure_seqs.ipynb)
    2. Detect contacts and highlight epitopes: [data/sabdab/02b_pandaprot_parallel.ipynb](data/sabdab/02b_pandaprot_parallel.ipynb)
    3. Generate the training dataset: [data/sabdab/03_generate_dataset.ipynb](data/sabdab/03_generate_dataset.ipynb)






## Fine-Tuning

```bash
python scripts/finetune.py
```