#!/bin/bash
# Run featureCounts on Ulixertinib/Ctrl BAMs (same params as PDC/PDO pipeline)
set -uo pipefail
ROOT=/home/susu/Rprojects/PCP202607
TURN=$ROOT/outputs/turn_20260802065502_636832
GTF=$ROOT/bulk2/refs/gencode.v44.primary_assembly.annotation.gtf
BAM_DIR=$TURN/03_Ulixertinib_vs_Ctrl/BAM
OUT=$TURN/03_Ulixertinib_vs_Ctrl/featureCounts
LOG=$TURN/03_Ulixertinib_vs_Ctrl/logs
mkdir -p "$OUT" "$LOG"

SAMPLES=(p3 p3u p4 p4u p5 p5u p6 p6u)
PIDS=()
MAX_PARALLEL=2

for s in "${SAMPLES[@]}"; do
  while (( ${#PIDS[@]} >= MAX_PARALLEL )); do
    NEW=()
    for pid in "${PIDS[@]}"; do
      if kill -0 "$pid" 2>/dev/null; then NEW+=("$pid"); fi
    done
    PIDS=("${NEW[@]}")
    [ ${#PIDS[@]} -ge $MAX_PARALLEL ] && sleep 5
  done
  echo "[$(date -Is)] start $s"
  nohup bash -c "
    /home/susu/Miniforge3/conda/envs/omics/bin/featureCounts \
      -T 16 -p -B -C --countReadPairs \
      -a '$GTF' -t exon -g gene_id \
      --extraAttributes gene_name,gene_type \
      -o '$OUT/${s}.featureCounts.tsv' \
      '$BAM_DIR/${s}.bam' \
      > '$LOG/${s}.featureCounts.log' 2>&1
    echo \"[\$(date -Is)] $s rc=\$?\" >> '$LOG/featureCounts_run.log'
  " >/dev/null 2>&1 &
  PIDS+=($!)
done
wait
echo "[$(date -Is)] ALL DONE"