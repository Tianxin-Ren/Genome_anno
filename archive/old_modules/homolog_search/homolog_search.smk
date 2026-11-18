# Miscellaneous Analysis Module: Homolog Search
# Includes: Batch Homolog Search (BLASTP + HMMER)

import os

# 1. Setup & Path Logic
# ------------------------------------------------------------------------------
c = config.get("homolog_search", {}).get("homolog_search", {})
SEARCH_PROJ = os.path.abspath(config.get("homolog_search", {}).get("project_dir", "."))

QUERY_DIR = os.path.join(SEARCH_PROJ, c.get("query_dir", "data/queries"))
QUERY_EXT = c.get("query_ext", ".faa")
MAG_FAA_DIR = os.path.join(SEARCH_PROJ, c.get("mag_faa_dir", "results/mag_annotation"))
SEARCH_OUT = os.path.join(SEARCH_PROJ, c.get("output_dir", "results/homolog_search"))
SEARCH_THREADS = c.get("threads", 8)
SEARCH_EVALUE = c.get("evalue", 1e-2)
SEARCH_ENV = config.get("homolog_search", {}).get("conda_envs", {}).get("bio", "bio")

# Identify all query files
if os.path.exists(QUERY_DIR):
    QUERIES, = glob_wildcards(os.path.join(QUERY_DIR, "{query}" + QUERY_EXT))
else:
    QUERIES = []

# Helper to find all .faa files for database
def get_mag_faa_inputs(wildcards):
    if os.path.isdir(MAG_FAA_DIR):
        import glob
        return glob.glob(os.path.join(MAG_FAA_DIR, "**/*.faa"), recursive=True)
    elif os.path.isfile(MAG_FAA_DIR):
        return [MAG_FAA_DIR]
    return []

# 2. Target Registration
# ------------------------------------------------------------------------------
TARGETS.extend(expand(os.path.join(SEARCH_OUT, "{query}_blast_results.tsv"), query=QUERIES))
TARGETS.extend(expand(os.path.join(SEARCH_OUT, "{query}_jackhmmer_results.tsv"), query=QUERIES))

# 3. Rules
# ------------------------------------------------------------------------------

rule cat_mags_protein:
    input:
        mags = get_mag_faa_inputs
    output:
        combined = temp(os.path.join(SEARCH_OUT, "tmp/all_mags_combined.faa"))
    shell:
        """
        rm -f {output.combined}
        mkdir -p $(dirname {output.combined})
        for f in {input.mags}; do
            mag_name=$(basename "$f" .faa)
            sed "s/^>/>${mag_name}|/" "$f" >> {output.combined}
        done
        """

rule blast_homologs:
    input:
        query = os.path.join(QUERY_DIR, "{query}" + QUERY_EXT),
        db_faa = os.path.join(SEARCH_OUT, "tmp/all_mags_combined.faa")
    output:
        result = os.path.join(SEARCH_OUT, "{query}_blast_results.tsv")
    log:
        os.path.join(SEARCH_OUT, "logs/blast_{query}.log")
    threads: SEARCH_THREADS
    params:
        env = SEARCH_ENV
    shell:
        """
        if [ ! -f {input.db_faa}.pdb ]; then
            conda run -p {params.env} makeblastdb -in {input.db_faa} -dbtype prot -out {input.db_faa}.db > /dev/null
        fi
        
        conda run -p {params.env} blastp -query {input.query} \
               -db {input.db_faa}.db \
               -out {output.result} \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
               -evalue {SEARCH_EVALUE} \
               -num_threads {threads} 2>&1 | tee {log}
        """

rule hmm_homologs:
    input:
        query = os.path.join(QUERY_DIR, "{query}" + QUERY_EXT),
        db_faa = os.path.join(SEARCH_OUT, "tmp/all_mags_combined.faa")
    output:
        result = os.path.join(SEARCH_OUT, "{query}_jackhmmer_results.tsv")
    log:
        os.path.join(SEARCH_OUT, "logs/jackhmmer_{query}.log")
    threads: SEARCH_THREADS
    params:
        env = SEARCH_ENV
    shell:
        "conda run -p {params.env} jackhmmer --cpu {threads} -E {SEARCH_EVALUE} --tblout {output.result} {input.query} {input.db_faa} > /dev/null 2> {log}"
