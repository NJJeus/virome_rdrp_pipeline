# USAGE
```bash
bash ./prepare_datasets.sh
snakemake -s neordrp_pipeline.py --use-conda -j 10 -c 10 --config input_directory="TEST_INPUT/" hmm_folder="HMMs/" hmm_info="Annotation_of_Seed_RdRp_datasets.tsv" suffix='.fa' hmm_threads=10 output_directory="TEST_OUTPUT"
```
