# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h4 align="right">Nicholas Santolla, Trey Pridgen, Prbhuv Nigam, and Colby T. Ford</h4><h3 align="right">Silico Biosciences</h3>


## About
`Peleke-1` is a suite of antibody language models that were fine-tuned to generate antibody sequences that specifically target given antigen sequences. By leveraging advanced protein language models and general large language models, each `peleke-1` model aims to streamline the process of *in silico* antibody design.

## Generate Antibody Sequences

```python
model_name = 'silicobio/peleke-phi-4'
config = PeftConfig.from_pretrained(model_name)

tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)

model = AutoModelForCausalLM.from_pretrained(config.base_model_name_or_path, torch_dtype=torch.bfloat16, trust_remote_code=True).cuda()
model.resize_token_embeddings(len(tokenizer))
model = PeftModel.from_pretrained(model, model_name).cuda()

```

Currently, the supported models are:
- [`peleke-phi-4`](https://huggingface.co/silicobio/peleke-phi-4), based on [Microsoft's Phi-4](https://huggingface.co/microsoft/phi-4) model.
- [`peleke-llama-3.1-8b-instruct`](https://huggingface.co/silicobio/peleke-llama-3.1-8b-instruct), based on [Meta's Llama 3.1 8B Instruct](https://huggingface.co/meta-llama/Llama-3.1-8B) model.
- [`peleke-mistral-7b-instruct-v0.2`](https://huggingface.co/silicobio/peleke-mistral-7b-instruct-v0.2), based on [Mistral's 7B Instruct v0.2](https://huggingface.co/mistralai/Mistral-7B-Instruct-v0.2) model.

You can also fine-tune your own `peleke-1`-like model by following the fine-tuning logic that can be found under [scripts/](scripts/).

### Tokenization

The `peleke-1` suite of models expects an amino acid sequence of an antigen protein as an input. Epitope residues should be enclosed in `<epi>` and `</epi>>` tokens.
However, if you prefer to use square brackets `[ ]`, which are easier, use the following function:

```python
def format_prompt(antigen_sequence):
    epitope_seq = re.sub(r'\[([A-Z])\]', r'<epi>\1</epi>', antigen_sequence)
    formatted_str = f"Antigen: {epitope_seq}<|im_end|>\nAntibody:"
    return formatted_str
```

For example, `AAM[K][R]HGL[D][N][Y]RG` will get formatted as `AAM<epi>K</epi><epi>R</epi>HGL<epi>D</epi><epi>N</epi><epi>Y</epi>RG`, using `<epi>` and `</epi>` as special tags.


## Training Dataset

The training dataset consists of paired antigen and antibody sequences, where the antigen is the target for which the antibody is generated. This was curated from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab). Using [PandaProt](https://github.com/pritampanda15/PandaProt), epitope residues were highlighted in the antibody sequences using `[ ]`, which helped to tune the model to generate antibodies sequences that fold and bind to specific epitopes on the desired antigen. Note that multi-chain antigen sequences are delimited by `|` in the `antigen_ids` column, and the heavy and light chain antibody sequences are delimited by `|` in the `antibody_sequences` column. We also provide the Fv portions of the antibodies chains, which (for length consistency) were used to tune the models.


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

- See the prepared training dataset:
    - In this repo: [data/sabdab/sabdab_training_dataset.csv](data/sabdab/sabdab_training_dataset.csv)
    - On Hugging Face: [silicobio/peleke_antibody-antigen_sabdab](https://huggingface.co/datasets/silicobio/peleke_antibody-antigen_sabdab)
- Our data preparation scripts:
    1. Get sequences from PDB structures: [data/sabdab/01_get_structure_seqs.ipynb](data/sabdab/01_get_structure_seqs.ipynb)
    2. Detect contacts and highlight epitopes: [data/sabdab/02b_pandaprot_parallel.ipynb](data/sabdab/02b_pandaprot_parallel.ipynb)
    3. Generate the training dataset: [data/sabdab/03_generate_dataset.ipynb](data/sabdab/03_generate_dataset.ipynb)

