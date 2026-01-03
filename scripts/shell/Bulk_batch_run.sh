for i in $(cat /data/file_name.txt)
do
    sample_name="${i}"
    echo ${sample_name}
    sed -i "s/SampleList: \['[^']*'\]/SampleList: ['${sample_name}']/" /config/Bluk_analysis_config.yaml
    bash /scripts/shell/Bulk_run.sh
done
