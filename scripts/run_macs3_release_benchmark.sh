#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-${HOME}/benchmarks/macs3-release-benchmark}"
REPEATS="${REPEATS:-3}"
TIME_MODE="${TIME_MODE:-auto}"
BASELINE_REF="${BASELINE_REF:-v3.0.4}"

BASELINE_CMD="${BASELINE_CMD:-macs3}"
CURRENT_CMD="${CURRENT_CMD:-macs3}"

SRC_DIR="${ROOT}/src"
BASELINE_TREE="${BASELINE_TREE:-${SRC_DIR}/MACS-baseline}"
CURRENT_TREE="${CURRENT_TREE:-${SRC_DIR}/MACS-current}"
CURRENT_REF="${CURRENT_REF:-$(git -C "${CURRENT_TREE}" rev-parse --short HEAD)}"
DATA_DIR="${ROOT}/data"
RESULTS_DIR="${ROOT}/results"
LOGS_DIR="${ROOT}/logs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TREATMENT="${DATA_DIR}/CTCF_12878_5M.bed.gz"
CONTROL="${DATA_DIR}/Input_12878_5M.bed.gz"
RUNS_TSV="${RESULTS_DIR}/benchmark_runs.tsv"
FAILED=0

quote_arg() {
  printf "%q" "$1"
}

time_mode() {
  if [[ "${TIME_MODE}" != "auto" ]]; then
    printf "%s\n" "${TIME_MODE}"
    return
  fi

  case "$(uname -s)" in
    Darwin) printf "macos\n" ;;
    Linux) printf "linux\n" ;;
    *)
      echo "Unsupported OS for timing: $(uname -s)" >&2
      return 1
      ;;
  esac
}

time_command() {
  case "$(time_mode)" in
    macos) printf "/usr/bin/time -l" ;;
    linux) printf "/usr/bin/time -v" ;;
    *)
      echo "Unsupported TIME_MODE=${TIME_MODE}" >&2
      return 1
      ;;
  esac
}

parse_wall_seconds() {
  local file="$1"

  if [[ "$(time_mode)" == "linux" ]]; then
    awk -F': ' '/Elapsed \(wall clock\) time/ { print elapsed_to_seconds($2) }
      function elapsed_to_seconds(value, parts, n, days, hours, minutes, seconds, dayparts) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
        if (value ~ /-/) {
          split(value, dayparts, "-")
          days = dayparts[1] + 0
          value = dayparts[2]
        }
        n = split(value, parts, ":")
        if (n == 3) {
          hours = parts[1] + 0
          minutes = parts[2] + 0
          seconds = parts[3] + 0
        } else if (n == 2) {
          minutes = parts[1] + 0
          seconds = parts[2] + 0
        } else {
          seconds = value + 0
        }
        return (days * 86400) + (hours * 3600) + (minutes * 60) + seconds
      }' "${file}" | tail -1
  else
    awk '/[[:space:]]real$/ || $2 == "real" { value = $1 } END { print value }' "${file}"
  fi
}

parse_peak_rss() {
  local file="$1"

  if [[ "$(time_mode)" == "linux" ]]; then
    awk -F': ' '/Maximum resident set size/ { value = $2 * 1024 } END { printf "%.0f\n", value }' "${file}"
  else
    awk '/maximum resident set size/ { value = $1 } END { print value }' "${file}"
  fi
}

prepare_inputs() {
  mkdir -p "${DATA_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}"

  cp "${CURRENT_TREE}/test/CTCF_12878_5M.bed.gz" "${TREATMENT}"
  cp "${CURRENT_TREE}/test/Input_12878_5M.bed.gz" "${CONTROL}"

  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${TREATMENT}" "${CONTROL}" > "${LOGS_DIR}/input_sha256.txt"
  else
    shasum -a 256 "${TREATMENT}" "${CONTROL}" > "${LOGS_DIR}/input_sha256.txt"
  fi
}

record_metadata() {
  {
    echo "## baseline (${BASELINE_REF})"
    git -C "${BASELINE_TREE}" rev-parse HEAD
    git -C "${BASELINE_TREE}" describe --tags --always --dirty
    echo
    echo "## current (${CURRENT_REF})"
    git -C "${CURRENT_TREE}" rev-parse HEAD
    git -C "${CURRENT_TREE}" describe --tags --always --dirty
  } > "${LOGS_DIR}/git_revisions.txt"

  {
    echo "## baseline (${BASELINE_REF})"
    bash -lc "${BASELINE_CMD} --version"
    echo
    echo "## current (${CURRENT_REF})"
    bash -lc "${CURRENT_CMD} --version"
  } > "${LOGS_DIR}/version_checks.txt"
}

append_run_row() {
  local revision="$1"
  local ref="$2"
  local rep="$3"
  local order="$4"
  local cmd="$5"
  local status="$6"
  local wall="$7"
  local rss="$8"
  local log="$9"

  cmd="${cmd//$'\t'/ }"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${revision}" "${ref}" "${rep}" "${order}" "${cmd}" "${status}" \
    "${wall}" "${rss}" "${log}" >> "${RUNS_TSV}"
}

run_one() {
  local revision="$1"
  local ref="$2"
  local rep="$3"
  local order="$4"
  local executable="$5"
  local outdir="${RESULTS_DIR}/${revision}_run_${rep}"
  local logbase="${LOGS_DIR}/${revision}_rep${rep}"
  local treatment_q
  local control_q
  local outdir_q
  local cmd
  local status
  local wall
  local rss

  rm -rf "${outdir}"
  mkdir -p "${outdir}"
  treatment_q="$(quote_arg "${TREATMENT}")"
  control_q="$(quote_arg "${CONTROL}")"
  outdir_q="$(quote_arg "${outdir}")"
  cmd="${executable} callpeak -t ${treatment_q} -c ${control_q} -f BED -g hs -n macs3_ctcf_5m --outdir ${outdir_q} -q 0.01"

  echo "Running ${revision} (${ref}) repeat ${rep}, order ${order}"
  set +e
  $(time_command) bash -lc "${cmd}" \
    > "${logbase}.stdout.txt" \
    2> "${logbase}.time_stderr.txt"
  status=$?
  set -e

  wall="$(parse_wall_seconds "${logbase}.time_stderr.txt")"
  rss="$(parse_peak_rss "${logbase}.time_stderr.txt")"
  append_run_row "${revision}" "${ref}" "${rep}" "${order}" "${cmd}" \
    "${status}" "${wall}" "${rss}" "${logbase}.time_stderr.txt"

  if [[ "${status}" -ne 0 ]]; then
    FAILED=1
  fi
}

run_benchmarks() {
  local rep

  printf "revision\tref\trep\torder\tcommand\tstatus\tseconds_wall\tmax_rss_bytes\tlog\n" > "${RUNS_TSV}"

  # Alternate the order to reduce bias from runner warm-up and drift.
  for rep in $(seq 1 "${REPEATS}"); do
    if (( rep % 2 == 1 )); then
      run_one baseline "${BASELINE_REF}" "${rep}" 1 "${BASELINE_CMD}"
      run_one current "${CURRENT_REF}" "${rep}" 2 "${CURRENT_CMD}"
    else
      run_one current "${CURRENT_REF}" "${rep}" 1 "${CURRENT_CMD}"
      run_one baseline "${BASELINE_REF}" "${rep}" 2 "${BASELINE_CMD}"
    fi
  done

  python3 "${SCRIPT_DIR}/summarize_macs3_release_benchmark.py" "${ROOT}"
}

main() {
  if ! [[ "${REPEATS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "REPEATS must be a positive integer, got: ${REPEATS}" >&2
    return 2
  fi

  prepare_inputs
  record_metadata
  run_benchmarks

  if [[ "${FAILED}" -ne 0 ]]; then
    echo "At least one benchmark command failed; inspect benchmark_runs.tsv and logs." >&2
    return 1
  fi
}

main "$@"
