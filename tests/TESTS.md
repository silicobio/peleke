
1. Numbering Test:

```bash
/opt/anaconda3/bin/python3 run_numbering_test.py
```

This script tests each generated sequence to see if it can be numbered using the Chothia antibody numbering system.


2. Structural Prediction:

- Run `post_analyses.ipynb` in a Jupyter Notebook environment with access to the ESM multimer model and an ESM token.


<!-- ```bash
docker pull cford38/openmm:cuda12.5.0
docker run -v ./tests:/mnt/tests --name openmm --rm -it cford38/openmm:cuda12.5.0 /bin/bash
cd /mnt/tests
python run_amber_relax.py /mnt/tests/structures/predicted_complexes/ /mnt/tests/structures/relaxed_complexes/
``` -->


3. Binding Affinity Prediction:

```bash
# docker pull cford38/haddock:3-2024.10.0b6
docker run -v ./tests:/mnt/tests --name haddock3 --rm -it cford38/haddock:3-2024.10.0b6 /bin/bash
cd /mnt/tests
python run_binding_affinity.py
```