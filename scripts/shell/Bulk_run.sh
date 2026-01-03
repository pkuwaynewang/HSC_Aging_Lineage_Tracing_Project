# run.sh
#
# This script runs a single execution of the pipeline
# using parameters specified in config/config.yaml.
# It assumes that the configuration file has already
# been prepared for the desired sample.


snakemake \
-s /path/to/snakefile_integrate_CARLIN.py \
--configfile /config/Bluk_analysis_config.yaml \
--core 10 \
--config \
sbatch=0 \
-R CARLIN 
