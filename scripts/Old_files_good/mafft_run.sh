mkdir -p data/processed/03_ml_orthogroups/06_alignments
LOG=data/processed/03_ml_orthogroups/logs/mafft_run.log
: > "$LOG"

for f in data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta; do
  og=$(basename "$f" .fasta)
  n=$(grep -c '^>' "$f")

  if [ "$n" -lt 3 ]; then
    echo "[SKIP] $og : only $n sequences" | tee -a "$LOG"
    continue
  fi

  # Choose a safe algorithm by size
  if [ "$n" -le 200 ]; then
    algo="--maxiterate 1000 --localpair"   # L-INS-i (high accuracy)
  else
    algo="--auto"                           # scalable
  fi

  echo "[RUN] $og : $n sequences : mafft $algo" | tee -a "$LOG"
  mafft $algo "$f" > data/processed/03_ml_orthogroups/06_alignments/${og}.aln.fasta
done
