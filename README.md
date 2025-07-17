# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h4 align="right">Trey Pridgen, Nicholas Santolla, Prbhuv Nigam, and Colby T. Ford</h4><h3 align="right">Silico Biosciences</h3>


## Training Dataset

The training dataset consists of paired antigen and antibody sequences, where the antigen is the target for which the antibody is generated. This was curated from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab).


| antigen               | antibody                                           |
| --------------------- | -------------------------------------------------- |
| `MKT[LLI]LAV[AA]A...` | `QVQLVQSGAEVKKPGAS...\|DIQMTQSPSSLSASVGDRVTITC...` |

See the dataset: [data/sabdab/sabdab_sequences.csv](data/sabdab/sabdab_sequences.csv)



## Generate Antibody Sequences

```bash
python scripts/generate.py --antigen "MKT[LLI]LAV[AA]A..." --model "peleke-phi-4"
```


## Fine-Tuning

```bash
python scripts/finetune.py
```