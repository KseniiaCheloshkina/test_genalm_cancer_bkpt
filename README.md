## Data Sources
Breakpoints data - https://github.com/KseniiaCheloshkina/cancer_breakpoints_hotspots_prediction/tree/master/data/raw%20breakpoints
Genome sequence (hg38) - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz (Dec 2013)
Chromosome lengths - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

Get DNA sequence by coordinates - other methods (curl):
* ucsc
```bash
curl http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr1:4321,4330
```
* api ucsc
```bash
curl -L 'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=4321;end=4330'
```

## Data preprocessing pipeline

1) Creates bed files to convert hg17 coordinates of breakppoints to hg38 coordinates
``` bash
python src/convert_to_bed_to_get_hg38_coordinates.py
```
2) Conversion is done in https://genome.ucsc.edu/cgi-bin/hgLiftOver

3) Removes blacklisted regions from the list of breakpoints
``` bash
python src/filter_bad_breakpoints.py
```

4) Find genome windows of specified length around breakpoint (positive examples) or randomly - uniform by genome and chromosome (negative examples)
``` bash
python src/generate_windows.py --win_len 512
python src/generate_windows.py --win_len 4000
```

5) Get DNA sequences for coordinates using BEDTOOLS - for 512 and 4000 window length 
```bash
apt-get update
apt-get install bedtools
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
bedtools getfasta -fi hg38.fa -bed positive_all_cancers_4000.bed -tab -fo pos_4000.bed
```

6) Collect final dataset with class balance 1:`n_times_neg_more` (positive: negative). Removes excluded regions from negatives
``` bash
python src/create_datasets.py --n_times_neg_more 1 --win_len 512
python src/create_datasets.py --n_times_neg_more 1 --win_len 4000
```
## Using code

Run tests:
```bash
python -m pytest
```