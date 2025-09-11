import pandas as pd
from Bio import SeqIO
from tqdm.auto import tqdm



## SNAKEMAKE I/O
reference_file = snakemake.input.reference
fastq_file = snakemake.input.fastq
mapped_output = snakemake.output.mapped
unmapped_output = snakemake.output.unmapped
max_effective_mismatches = snakemake.params.max_effective_mismatches
max_mismatches = snakemake.params.max_mismatches
quality_threshold = snakemake.params.quality_threshold
part = snakemake.wildcards.part  # 'barcodes' or 'guides'
sample = snakemake.wildcards.sample

## UTILS
def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def is_true_mismatch(base1, base2, qual, qual_threshold):
    if base1 != base2 and qual >= qual_threshold:
        return True
    return False

## LOAD REFERENCE
assert part in ('barcodes', 'guides'), f"part must be 'barcodes' or 'guides', got {part}"
reference_name_key = "barcode" if part == 'barcodes' else "gRNA"
reference = pd.read_csv(reference_file, index_col=reference_name_key).guide_name.to_dict()

## RUN MAPPING
res = {}
unmatched = []
for record in tqdm(SeqIO.parse(fastq_file, "fastq")):
    qual = record.letter_annotations["phred_quality"]
    seq = str(record.seq)
    best_match = None
    for key in reference.keys():
        mismatches = hamming_distance(seq, key)
        if mismatches > max_mismatches:
            continue  # Skip, too many mismatches mismatches)
        effective_mismatches = sum(
                is_true_mismatch(a, b, q, quality_threshold)
                for a, b, q in zip(seq, key, qual)
            )
        if effective_mismatches > max_effective_mismatches:
            continue  # Skip, too many effective mismatches
        if best_match is None or (mismatches, effective_mismatches) < (best_match[1], best_match[2]):
            best_match = (key, mismatches, effective_mismatches)
    if best_match:
        res[record.id] = best_match
    else:
        unmatched.append(record.id)

# SAVE OUTPUT
res_df = pd.DataFrame.from_dict(res, orient='index', columns=['match', 'mismatches', 'effective_mismatches'])
res_df.to_csv(mapped_output)
unmatched_df = pd.DataFrame(unmatched, columns=['read_id'])
unmatched_df.to_csv(unmapped_output, index=False)