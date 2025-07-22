# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h4 align="right">Trey Pridgen, Nicholas Santolla, Prbhuv Nigam, and Colby T. Ford</h4><h3 align="right">Silico Biosciences</h3>


## Generate Antibody Sequences

```bash
python scripts/generate.py --antigen "MKT[LLI]LAV[AA]A..." --model "peleke-phi-4"
```

## Training Dataset

The training dataset consists of paired antigen and antibody sequences, where the antigen is the target for which the antibody is generated. This was curated from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab). Using PandaProt, epitope residues were highlighted in the antibody sequences using `[ ]`, which were used to tune the model to generate antibodies sequences that fold and bind to specific epitopes on the desired antigen.


| **pdb_id** | **h_chain_id** | **l_chain_id** | **antigen_ids** | **h_chain_seq** | **l_chain_seq** | **antigen_seqs** |    **antibody_sequences**    | **highlighted_epitope_sequences** |          **epitope_residues**         |
|:----------:|:--------------:|:--------------:|:---------------:|:---------------:|:---------------:|:----------------:|:----------------------------:|:---------------------------------:|:-------------------------------------:|
| 8xa4       | C              | D              | A\|B            | QLQLQESGPG...   | EIVLTQSPGT...   | SCNGLYYQGSC...   | QLQLQESGPG...\|EIVLTQSPGT... | ...LTTWLI[D][Y]V[E][D][T]WGS...   | A:ARG 176\|A:ASP 146\|A:ASP 150...    |
| 9cph       | H              | L              | A               | EVQLVESGGG...   | AQMTQSPSSL...   | KIEEGKLVIWI...   | EVQLVESGGG...\|AQMTQSPSSL... | ...FQDKL[Y][P][F]TW[D][A]VRY...   | A:ALA 1116\|A:ALA 1122\|A:ALA 1128... |
| 9d7i       | H              | G              | E               | VQLQESGPGV...   | YELTQPPSVS...   | LWVTVYYGVPV...   | VQLQESGPGV...\|YELTQPPSVS... | ...DVVQI[N]NKEYRL[I]NC[N][T]...   | E:ARG 429\|E:ARG 469\|E:ASN 177...    |
| 9d7i       | J              | I              | C               | VQLQESGPGV...   | YELTQPPSVS...   | LWVTVYYGVPV...   | VQLQESGPGV...\|YELTQPPSVS... | ...YRL[I]NC[N][T]SACTQACPKVS...   | C:ARG 469\|C:ASN 197\|C:ASN 280...    |
| 9d7o       | H              | G              | E               | QVQLQESGPG...   | YELTQPPSVS...   | LWVTVYYGVPV...   | QVQLQESGPG...\|YELTQPPSVS... | ...VQIN[E]SNKEYRL[I]NC[N]TSA...   | E:ARG 429\|E:ARG 469\|E:ASN 197...    |

See the dataset: [data/sabdab/sabdab_training_dataset.csv](data/sabdab/sabdab_training_dataset.csv)






## Fine-Tuning

```bash
python scripts/finetune.py
```