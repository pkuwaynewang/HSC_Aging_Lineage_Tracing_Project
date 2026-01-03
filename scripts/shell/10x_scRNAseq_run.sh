snakemake \
-s /path/to/snakemake_DARLIN/snakefiles/snakefile_integrate_CARLIN.py  \
--configfile /config/10x_scRNAseq_config.yaml \
--core 10 \
--config sbatch=0 \
-R CARLIN 
