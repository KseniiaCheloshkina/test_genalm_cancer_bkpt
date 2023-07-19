## Data Sources
Genome sequence - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
Chromosome lengths - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

Methods to get DNA sequence by coordinates in genome:
* pyfaidx:
```python
from pyfaidx import Fasta
seq = Fasta("data/ncbi_dataset/ncbi_dataset/data/GCA_000001405.15/GCA_000001405.15_GRCh38_genomic.fna")
```
* pybedtools
* ucsc
```bash
curl http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr1:4321,4330
```
* api ucsc
```bash
curl -L 'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=4321;end=4330'
```

## Using code

Run tests:
```bash
python -m pytest
```