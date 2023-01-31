# README

Link data directory and ignore from git tracking
```bash
ln -s /data/CARD/projects/proteomicsXprogression_PD/data-allreleases data
echo 'data' >> .gitignore
```

Generate summary table [`mds_scores.tsv`](reports/mds_scores.tsv) for MDS scores from [`combine_mds_scores.R`](src/combine_mds_scores.R)
```bash
module load R/3.6.3
Rscript src/combine_mds_scores.R
```

Call [`format_protein_data.R`](src/format_protein_data.R) to read and reformat various protein quantification matrices, concatenating into a single long-form file [`protein_quantification.tsv`](reports/protein_quantification.tsv)
```bash
module load R/3.6.3
Rscript src/combine_mds_scores.R
```