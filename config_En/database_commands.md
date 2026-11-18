# Database Commands

Run only the blocks you need. These commands are intentionally separated from
environment installation so existing databases are not touched accidentally.

Set the shared database root first:

```bash
DBROOT="/dssg/home/acct-clswfp/clswfp-icdli1/Database"
THREADS="${THREADS:-16}"
```

## 1. Confirm GTDB-Tk database

Your check showed GTDB marker files already exist. Use this command to confirm
the installed release with GTDB-Tk itself:

```bash
conda run -n gtdbtk_2.7.1 env GTDBTK_DATA_PATH="$DBROOT/GTDB/archive/gtdb-tk/release232" gtdbtk check_install
```

Only if that check fails, download GTDB-Tk data:

```bash
mkdir -p "$DBROOT/GTDB/archive/gtdb-tk/release232"
cd "$DBROOT/GTDB/archive/gtdb-tk/release232"
curl -L -O https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz
```

## 2. Download Rfam for Infernal

This is the main missing RNA database from your check.

```bash
mkdir -p "$DBROOT/Rfam/archive/15.1"
cd "$DBROOT/Rfam/archive/15.1"
curl --fail --location --retry 5 --retry-delay 10 --remote-name \
  https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
curl --fail --location --retry 5 --retry-delay 10 --remote-name \
  https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
gunzip -f Rfam.cm.gz
conda run -n infernal_1.1.5 cmpress Rfam.cm
```

## 3. Download Swiss-Prot and build MMseqs2 DB

Run this only if `Database/swiss-prot` has no usable FASTA/DB files.

```bash
DBROOT="/dssg/home/acct-clswfp/clswfp-icdli1/Database"
THREADS="${THREADS:-16}"
mkdir -p "$DBROOT/swiss-prot" "$DBROOT/mmseqs2_db"
cd "$DBROOT/swiss-prot"
curl -L -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
conda run -n mmseqs2_17.b804f mmseqs createdb \
  "$DBROOT/swiss-prot/uniprot_sprot.fasta.gz" \
  "$DBROOT/mmseqs2_db/swiss-prot"
conda run -n mmseqs2_17.b804f mmseqs createindex \
  "$DBROOT/mmseqs2_db/swiss-prot" \
  "$DBROOT/mmseqs2_db/tmp_swiss-prot" \
  --threads "$THREADS" --split-memory-limit 8G
```

## 4. Optional: CheckM2 database

Your check showed this already exists. Only run if missing:

```bash
conda run -n checkm2_1.0.2 checkm2 database --download --path "$DBROOT/CheckM2_database"
```

## 5. Optional: eggNOG database

Your check showed this already exists. Only run if missing:

```bash
mkdir -p "$DBROOT/eggnog"
conda run -n eggnog-mapper download_eggnog_data.py -y -M --data_dir "$DBROOT/eggnog"
```

## 6. Optional: KOfam database

Your check showed this already exists. Only run if missing:

```bash
mkdir -p "$DBROOT/kofam"
cd "$DBROOT/kofam"
curl -L -O https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
curl -L -O https://www.genome.jp/ftp/db/kofam/ko_list.gz
tar -xzf profiles.tar.gz
gunzip -f ko_list.gz
```

## 7. Optional large database: UniRef90

Do not run unless you explicitly want this large database.

```bash
mkdir -p "$DBROOT/UniRef90/current_release" "$DBROOT/mmseqs2_db"
cd "$DBROOT/UniRef90/current_release"
curl -L -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz
conda run -n mmseqs2_17.b804f mmseqs createdb \
  "$DBROOT/UniRef90/current_release/uniref90.fasta.gz" \
  "$DBROOT/mmseqs2_db/uniref90"
conda run -n mmseqs2_17.b804f mmseqs createindex \
  "$DBROOT/mmseqs2_db/uniref90" \
  "$DBROOT/mmseqs2_db/tmp_uniref90" \
  --threads "$THREADS"
```

## 8. Optional large database: RefSeq bacteria/archaea WP proteins

Do not run unless you explicitly want this large database.

```bash
mkdir -p "$DBROOT/RefSeq/protein"
cd "$DBROOT/RefSeq/protein"
wget -c -nd -r -l1 -A "bacteria.wp_protein.*.protein.faa.gz" ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/
wget -c -nd -r -l1 -A "archaea.wp_protein.*.protein.faa.gz" ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/
zcat bacteria.wp_protein.*.protein.faa.gz archaea.wp_protein.*.protein.faa.gz > refseq_bacteria_archaea.faa
conda run -n mmseqs2_17.b804f mmseqs createdb \
  refseq_bacteria_archaea.faa \
  refseq_bacteria_archaea.mmseqs
conda run -n mmseqs2_17.b804f mmseqs createindex \
  refseq_bacteria_archaea.mmseqs \
  tmp_refseq \
  --threads "$THREADS"
```

## 9. Create BioB local resource directories

Create directories for your curated BioB HMMs and references.

```bash
mkdir -p "$DBROOT/BioB/hmm" \
  "$DBROOT/BioB/reference_proteins" \
  "$DBROOT/BioB/references/bacteria" \
  "$DBROOT/BioB/references/archaea"
```

To run the formal BioB HMM evidence channel, place curated HMMER3 profile HMM
files under `$DBROOT/BioB/hmm/source_hmms/`, then concatenate and press them:

```bash
mkdir -p "$DBROOT/BioB/hmm/source_hmms"

# After copying curated *.hmm files into source_hmms:
mapfile -d '' BIOB_HMMS < <(find "$DBROOT/BioB/hmm/source_hmms" -type f -name "*.hmm" -print0 | sort -z)
if [ "${#BIOB_HMMS[@]}" -eq 0 ]; then
  echo "No curated BioB *.hmm files found in $DBROOT/BioB/hmm/source_hmms" >&2
  exit 1
fi
cat "${BIOB_HMMS[@]}" > "$DBROOT/BioB/hmm/biob_core_functions.hmm"
conda run -n hmmer_3.4 hmmpress -f "$DBROOT/BioB/hmm/biob_core_functions.hmm"
```

If the curated resource is already a single combined HMM file, copy or symlink
it directly:

```bash
ln -sf /path/to/curated/biob_core_functions.hmm \
  "$DBROOT/BioB/hmm/biob_core_functions.hmm"
conda run -n hmmer_3.4 hmmpress -f "$DBROOT/BioB/hmm/biob_core_functions.hmm"
```
