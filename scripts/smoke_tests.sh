#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_BASE="${TMPDIR:-/tmp}"
OUT_DIR="$(mktemp -d "${TMP_BASE%/}/gwasTools-smoke.XXXXXX")"

cleanup() {
  rm -rf "${OUT_DIR}"
}
trap cleanup EXIT

mkdir -p "${OUT_DIR}"
cd "${ROOT_DIR}"

run() {
  printf '+'
  for arg in "$@"; do
    printf ' %q' "${arg}"
  done
  printf '\n'
  "$@"
}

check_file() {
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    echo "Expected non-empty output missing: ${path}" >&2
    return 1
  fi
}

expect_fail_contains() {
  local label="$1"
  local pattern="$2"
  shift 2
  local log_path="${OUT_DIR}/${label//[^A-Za-z0-9]/_}.log"
  set +e
  "$@" >"${log_path}" 2>&1
  local status=$?
  set -e
  if [[ ${status} -eq 0 ]]; then
    echo "Expected command to fail (${label}), but it succeeded" >&2
    return 1
  fi
  if ! grep -q "${pattern}" "${log_path}"; then
    echo "Expected failure (${label}) to contain: ${pattern}" >&2
    sed -n '1,80p' "${log_path}" >&2
    return 1
  fi
}

run Rscript -e 'required <- c("optparse", "data.table", "RColorBrewer", "plotrix", "pwr"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) { message("Missing R packages: ", paste(missing, collapse = ", ")); quit(status = 1) }'

run Rscript QQplot.r --help >/dev/null
run Rscript ManhattanPlot.r --help >/dev/null
run Rscript PowerPlot.r --help >/dev/null

run Rscript QQplot.r \
  --input tests/fixtures/minimal_gwas.tsv \
  --prefix "${OUT_DIR}/smoke" \
  --maf MAF \
  --pvalue PVALUE \
  --maintitle "Smoke QQ" \
  --width 600 \
  --height 600 \
  --pointsize 12
check_file "${OUT_DIR}/smoke_QQ.png"

run Rscript QQplot.r \
  --input tests/fixtures/edge_gwas.tsv \
  --prefix "${OUT_DIR}/edge" \
  --maf MAF \
  --pvalue PVALUE \
  --allele_frequency_type auto \
  --maintitle "Edge QQ" \
  --width 600 \
  --height 600 \
  --pointsize 12
check_file "${OUT_DIR}/edge_QQ.png"

expect_fail_contains "qq declared maf rejects eaf" "MAF column contains values > 0.5" \
  Rscript QQplot.r \
    --input tests/fixtures/edge_gwas.tsv \
    --prefix "${OUT_DIR}/bad_maf" \
    --maf MAF \
    --pvalue PVALUE \
    --allele_frequency_type MAF

run Rscript ManhattanPlot.r \
  --input tests/fixtures/minimal_gwas.tsv \
  --prefix "${OUT_DIR}/smoke" \
  --pos POS \
  --pvalue PVALUE \
  --maintitle "Smoke Manhattan" \
  --width 900 \
  --height 600 \
  --pointsize 12 \
  --disable_thinning T
check_file "${OUT_DIR}/smoke_Manhattan.png"

run Rscript ManhattanPlot.r \
  --input tests/fixtures/edge_gwas.tsv \
  --prefix "${OUT_DIR}/edge" \
  --pos POS \
  --pvalue PVALUE \
  --maintitle "Edge Manhattan" \
  --width 900 \
  --height 600 \
  --pointsize 12 \
  --thinning_mode none
check_file "${OUT_DIR}/edge_Manhattan.png"
check_file "${OUT_DIR}/edge_ManhattanPlot_500kb_CandidateRegions.txt"

cat >"${OUT_DIR}/legacy_chr.tsv" <<'EOF'
CHR	POS	MAF	PVALUE
1	1000	0.10	0.50
2	2000	0.20	0.01
EOF

run Rscript ManhattanPlot.r \
  --input "${OUT_DIR}/legacy_chr.tsv" \
  --prefix "${OUT_DIR}/legacy_chr" \
  --pos POS \
  --pvalue PVALUE \
  --width 900 \
  --height 600 \
  --pointsize 12 \
  --disable_thinning T
check_file "${OUT_DIR}/legacy_chr_Manhattan.png"

cat >"${OUT_DIR}/mixed_chrom.tsv" <<'EOF'
CHROM	POS	MAF	PVALUE
1	1000	0.10	0.50
chr1	2000	0.20	0.01
X	3000	0.30	0.20
23	4000	0.40	0.30
EOF

run Rscript ManhattanPlot.r \
  --input "${OUT_DIR}/mixed_chrom.tsv" \
  --prefix "${OUT_DIR}/mixed_chrom" \
  --pos POS \
  --pvalue PVALUE \
  --width 900 \
  --height 600 \
  --pointsize 12 \
  --disable_thinning T
check_file "${OUT_DIR}/mixed_chrom_Manhattan.png"

cat >"${OUT_DIR}/single_chrom.tsv" <<'EOF'
CHROM	POS	MAF	PVALUE
1	1000	0.10	0.50
1	2000	0.20	0.01
EOF

run Rscript ManhattanPlot.r \
  --input "${OUT_DIR}/single_chrom.tsv" \
  --prefix "${OUT_DIR}/single_chrom" \
  --pos POS \
  --pvalue PVALUE \
  --width 900 \
  --height 600 \
  --pointsize 12 \
  --disable_thinning T
check_file "${OUT_DIR}/single_chrom_Manhattan.png"

cat >"${OUT_DIR}/invalid_pvalue.tsv" <<'EOF'
CHROM	POS	MAF	PVALUE
1	1000	0.10	0
1	2000	0.20	1.2
EOF

expect_fail_contains "qq invalid pvalue" "Invalid p-values" \
  Rscript QQplot.r \
    --input "${OUT_DIR}/invalid_pvalue.tsv" \
    --prefix "${OUT_DIR}/invalid_qq" \
    --maf MAF \
    --pvalue PVALUE

expect_fail_contains "manhattan invalid pvalue" "Invalid p-values" \
  Rscript ManhattanPlot.r \
    --input "${OUT_DIR}/invalid_pvalue.tsv" \
    --prefix "${OUT_DIR}/invalid_manhattan" \
    --pos POS \
    --pvalue PVALUE

cat >"${OUT_DIR}/unsupported_contig.tsv" <<'EOF'
CHROM	POS	MAF	PVALUE
GL000220.1	1000	0.10	0.5
EOF

expect_fail_contains "manhattan unsupported contig" "Unsupported chromosome" \
  Rscript ManhattanPlot.r \
    --input "${OUT_DIR}/unsupported_contig.tsv" \
    --prefix "${OUT_DIR}/bad_contig" \
    --pos POS \
    --pvalue PVALUE

expect_fail_contains "manhattan labelpeaks offline" "Nearest-gene labels require" \
  Rscript ManhattanPlot.r \
    --input tests/fixtures/edge_gwas.tsv \
    --prefix "${OUT_DIR}/labelpeaks" \
    --pos POS \
    --pvalue PVALUE \
    --gene_data_dir "${OUT_DIR}/gene_positions" \
    --labelPeaks T \
    --disable_thinning T
if [[ -e "${OUT_DIR}/gene_positions" ]]; then
  echo "Nearest-gene cache directory should not be created without --allow_gene_download T" >&2
  exit 1
fi

run Rscript PowerPlot.r \
  --prefix "${OUT_DIR}/smoke" \
  --cases 100,200 \
  --controls 150,300 \
  --minMAF 0.01 \
  --alpha 5E-8 \
  --stepsize 0.01 \
  --width 700 \
  --height 600 \
  --pointsize 12
check_file "${OUT_DIR}/smoke_PowerAnalysis.png"

run Rscript PowerPlot.r \
  --trait_type quantitative \
  --prefix "${OUT_DIR}/quantitative" \
  --samples 100,200 \
  --minMAF 0.01 \
  --alpha 5E-8 \
  --stepsize 0.01 \
  --width 700 \
  --height 600 \
  --pointsize 12
check_file "${OUT_DIR}/quantitative_PowerAnalysis.png"

expect_fail_contains "power mismatched sample sizes" "same number" \
  Rscript PowerPlot.r \
    --prefix "${OUT_DIR}/bad_power" \
    --cases 100,200 \
    --controls 100 \
    --stepsize 0.01

expect_fail_contains "quantitative power invalid samples" "greater than --covariates + 2" \
  Rscript PowerPlot.r \
    --trait_type quantitative \
    --prefix "${OUT_DIR}/bad_quantitative_power" \
    --samples 2 \
    --covariates 1 \
    --stepsize 0.01

echo "Smoke tests passed. Temporary outputs were written under ${OUT_DIR} and removed."
