# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h3 align="right">Silico Biosciences</h3>


## Training Dataset

The training dataset consists of paired antigen and antibody sequences, where the antigen is the target for which the antibody is generated. Thiss was curated from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab).


| antigen               | antibody                                          |
| --------------------- | ------------------------------------------------- |
| `MKT[LLI]LAV[AA]A...` | `QVQLVQSGAEVKKPGAS...\|DIQMTQSPSSLSASVGDRVTITC...` |


## Fine-Tuning

```bash
python scripts/finetune.py
```

## Generate Antibody Sequences

```bash
python scripts/generate.py --antigen "MKT[LLI]LAV[AA]A..."
```