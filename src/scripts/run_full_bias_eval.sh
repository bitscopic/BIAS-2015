# This script runs the full bias evaluation pipeline on the eRepo Nirvana output.
# You can uncomment out the Nirvana command and run it if you don't already have Nirvana output.
# You will need to edit the file paths to match your inputs and outputs and script locations.

#usage: run_full_bias_eval.sh output_prefix

output_prefix=$1
output_file="${output_prefix}_out.tsv"

# dotnet /home/bdn/software/Nirvana-v3.18.1/Nirvana.dll \
# --cache /home/bdn/data/Cache/GRCh37/Both --sd /home/bdn/data/SupplementaryAnnotation/GRCh37 \
# --ref /home/bdn/data/References/Homo_sapiens.GRCh37.Nirvana.dat \
# --in ./bias_sorted_vcf_file.vcf --o erepo_nirvana


python3 /home/bdn/repos/alternates/bias-2015/bias-2015.py \
    /home/bdn/bias_runs/bias_val/erepo_nirvana.json.gz \
    /home/bdn/repos/bias-2015/config/requiredPaths.json \
    ${output_file}

python3 /home/bdn/repos/bias-2015/src/scripts/evaluate_bias_erepo_validation.py \
    /home/bdn/bias_runs/bias_val/erepo_tabbed.txt \
    $output_file \
    /home/bdn/bias_runs/bias_val/bias_sorted_vcf_file.vcf \
    bias
