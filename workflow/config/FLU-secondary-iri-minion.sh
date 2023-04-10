# FLU-secondary-iri.sh
# Source: https://github.com/IRI-UW-Bioinformatics/flu-ngs
SKIP_E=1
SORT_PROG="LABEL BLAT"
RESIDUAL_ASSEMBLY_FACTOR=400
DO_SECONDARY=1


###################################
# Parameters set in FLU-minion.sh #
###################################

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=0			# average or median threshold for QUALITY reads
MIN_LEN=150				# minimum read length for QUALITY reads
INS_T=0.75				# threshold for insertion refinement
DEL_T=0.75				# threshold for deletion refinement
MIN_RP=3				# minimum read pattern count to continue
MIN_RC=3				# minimum read count to continue

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8			# minimum average variant quality, does not apply to deletions

# SORT_PROG="BLAT"  # DJP: Differ from FLU-minion.sh. Default to LABEL as in FLU-secondary-iri.sh
ALIGN_PROG="SAM BLAT"

SSW_M=2			# smith-waterman match score
SSW_X=3			# smith-waterman mismatch penalty
SSW_O=6			# smith-waterman gap open penalty
SSW_E=1			# smith-waterman gap extension penalty
