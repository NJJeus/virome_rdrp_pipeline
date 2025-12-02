import os
from glob import glob

# --------------------------
# Configuration
# --------------------------

# Load configuration with clear variable names
OUTPUT_DIR = config["output_directory"]
INPUT_DIR = config["input_directory"]
HMM_DIR = config["hmm_folder"]
HMM_INFO = config["hmm_info"]
READ_SUFFIX = config["suffix"]  # Single suffix for single-end reads
THREADS = config["hmm_threads"]


def get_samples():
    """Get sample names from input files, excluding special files - SINGLE-END VERSION"""
    files = glob(os.path.join(INPUT_DIR, "*" + READ_SUFFIX))
    samples = [os.path.basename(f).replace(READ_SUFFIX, "") for f in files]
    # Filter out any non-sample files (like reports, etc.)
    return [s for s in samples if not s.startswith(('report', 'general', 'summary'))]

# --------------------------
# Sample Detection
# --------------------------

# Get sample names from input files (single-end version)
SAMPLES, = glob_wildcards(os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX))
# --------------------------
# Rule Definitions
# --------------------------

rule all:
    input:
        # Individual sample reports first
        expand(os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv"), sample=SAMPLES)
        # Then general summary files

# --------------------------
# HMM Analysis
# --------------------------

rule hmm_scan:
    input:
        os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX)
    output:
        hits = os.path.join(OUTPUT_DIR, "tmp_hmm_results/{sample}.csv")
    params:
        hmm_dir = HMM_DIR
    threads: 10
    priority: 38
    conda:
        'envs/hmm_scan.yaml'
    shell:
        """
        python scripts/analyze_seq_hmm.py \
            -f {input} \
            -o {output.hits} \
            -m {params.hmm_dir} \
            -t {threads} \
            --batch 10000
        """

rule generate_hmm_report:
    input:
        hits = rules.hmm_scan.output.hits,
        contigs = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX)
    output:
        report = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv")
    params:
        hmm_info = HMM_INFO
    conda:
        'envs/hmm_scan.yaml'
    priority: 40
    shell:
        """
        python scripts/generate_hmm_report.py \
            -i {input.hits} \
            -o {output.report} \
            -m {params.hmm_info} \
        """


