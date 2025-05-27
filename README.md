# peleke 🦋

Fine-Tuned Protein Language Models for Targeted Antibody Sequence Generation.

<h3 align="right">Silico Biosciences</h3>


## Training Dataset

| antigen               | antibody                                          |
| --------------------- | ------------------------------------------------- |
| `MKT[LLI]LAV[AA]A...` | `QVQLVQSGAEVKKPGAS...\|DIQMTQSPSSLSASVGDRVTITC...` |


## Fine Tuning

```bash
python scripts/finetune.py
```

## Generate Antibody Sequences

```bash
python scripts/generate.py --antigen "MKT[LLI]LAV[AA]A..."
```