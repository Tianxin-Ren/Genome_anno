# External Genomes Annotation Dispatcher
# ------------------------------------------------------------------------------

include: "external_common.smk"

if AN_C.get("batch_hmm", False):
    # Define both standard and batch rules. They will be partitioned using wildcard constraints.
    include: "external_standard.smk"
    include: "external_batch.smk"
else:
    include: "external_standard.smk"
