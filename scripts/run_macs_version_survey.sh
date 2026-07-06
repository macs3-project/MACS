#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-${HOME}/benchmarks/macs-version-survey}"
REPEATS="${REPEATS:-3}"
MACS_REPO_URL="${MACS_REPO_URL:-https://github.com/macs3-project/MACS.git}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIME_MODE="${TIME_MODE:-auto}"

MACS1_REF="${MACS1_REF:-origin/macs_v1}"
MACS2_REF="${MACS2_REF:-origin/macs_v2}"
MACS3_REF="${MACS3_REF:-v3.0.4}"

MACS1_CMD="${MACS1_CMD:-macs14}"
MACS2_CMD="${MACS2_CMD:-macs2}"
MACS3_CMD="${MACS3_CMD:-macs3}"

SRC_DIR="${ROOT}/src"
MAIN_REPO="${SRC_DIR}/MACS"
MACS1_TREE="${SRC_DIR}/MACS-v1"
MACS2_TREE="${SRC_DIR}/MACS-v2"
MACS3_TREE="${MACS3_TREE:-${SRC_DIR}/MACS-v3.0.4}"
DATA_DIR="${ROOT}/data"
RESULTS_DIR="${ROOT}/results"
LOGS_DIR="${ROOT}/logs"

TREATMENT="${DATA_DIR}/CTCF_12878_5M.bed.gz"
CONTROL="${DATA_DIR}/Input_12878_5M.bed.gz"
RUNS_TSV="${RESULTS_DIR}/benchmark_runs.tsv"

quote_arg() {
  printf "%q" "$1"
}

prepare_dirs() {
  mkdir -p "${SRC_DIR}" "${ROOT}/envs" "${DATA_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}"
}

prepare_repo() {
  if [[ -d "${MAIN_REPO}/.git" ]]; then
    git -C "${MAIN_REPO}" fetch --all --tags --prune
  else
    git clone --recurse-submodules "${MACS_REPO_URL}" "${MAIN_REPO}"
    git -C "${MAIN_REPO}" fetch --all --tags --prune
  fi
}

ensure_worktree() {
  local path="$1"
  local ref="$2"

  if [[ -e "${path}/.git" ]]; then
    return 0
  fi

  if [[ -e "${path}" ]]; then
    echo "Refusing to create worktree because ${path} exists but is not a Git worktree." >&2
    return 1
  fi

  git -C "${MAIN_REPO}" worktree add "${path}" "${ref}"
  if [[ -f "${path}/.gitmodules" ]]; then
    git -C "${path}" submodule update --init --recursive
  fi
}

record_revisions() {
  {
    for path in "${MACS1_TREE}" "${MACS2_TREE}" "${MACS3_TREE}"; do
      echo "## ${path}"
      git -C "${path}" rev-parse HEAD
      git -C "${path}" describe --tags --always --dirty
    done
  } > "${LOGS_DIR}/git_revisions.txt"
}

copy_data() {
  local source_tree="${MACS3_TREE}"

  if [[ ! -f "${source_tree}/test/CTCF_12878_5M.bed.gz" || ! -f "${source_tree}/test/Input_12878_5M.bed.gz" ]]; then
    source_tree="${MAIN_REPO}"
  fi

  cp "${source_tree}/test/CTCF_12878_5M.bed.gz" "${DATA_DIR}/"
  cp "${source_tree}/test/Input_12878_5M.bed.gz" "${DATA_DIR}/"

  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${TREATMENT}" "${CONTROL}" > "${LOGS_DIR}/input_sha256.txt"
  else
    shasum -a 256 "${TREATMENT}" "${CONTROL}" > "${LOGS_DIR}/input_sha256.txt"
  fi
}

record_versions() {
  {
    echo "## MACS1"
    bash -lc "${MACS1_CMD} --help 2>&1 | head -5" || true
    echo

    echo "## MACS2"
    bash -lc "${MACS2_CMD} --version 2>&1 || ${MACS2_CMD} --help 2>&1 | head -5" || true
    echo

    echo "## MACS3"
    bash -lc "${MACS3_CMD} --version 2>&1 || ${MACS3_CMD} --help 2>&1 | head -5" || true
  } > "${LOGS_DIR}/version_checks.txt"
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

time_mode() {
  if [[ "${TIME_MODE}" != "auto" ]]; then
    printf "%s\n" "${TIME_MODE}"
    return
  fi

  case "$(uname -s)" in
    Darwin) printf "macos\n" ;;
    Linux) printf "linux\n" ;;
    *)
      echo "Unsupported OS for timing: $(uname -s). Set TIME_MODE=macos or TIME_MODE=linux." >&2
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

append_run_row() {
  local version="$1"
  local rep="$2"
  local cmd="$3"
  local status="$4"
  local wall="$5"
  local rss="$6"
  local log="$7"

  cmd="${cmd//$'\t'/ }"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${version}" "${rep}" "${cmd}" "${status}" "${wall}" "${rss}" "${log}" >> "${RUNS_TSV}"
}

run_one() {
  local version="$1"
  local rep="$2"
  local cmd="$3"
  local outdir="${RESULTS_DIR}/${version}_run_${rep}"
  local logbase="${LOGS_DIR}/${version}_rep${rep}"
  local status
  local wall
  local rss

  rm -rf "${outdir}"
  mkdir -p "${outdir}"

  echo "Running ${version} rep ${rep}"
  set +e
  $(time_command) bash -lc "${cmd}" \
    > "${logbase}.stdout.txt" \
    2> "${logbase}.time_stderr.txt"
  status=$?
  set -e

  wall="$(parse_wall_seconds "${logbase}.time_stderr.txt")"
  rss="$(parse_peak_rss "${logbase}.time_stderr.txt")"

  append_run_row "${version}" "${rep}" "${cmd}" "${status}" "${wall}" "${rss}" "${logbase}.time_stderr.txt"
}

run_benchmarks() {
  local treatment_q
  local control_q
  local outdir_q
  local cmd
  local rep

  treatment_q="$(quote_arg "${TREATMENT}")"
  control_q="$(quote_arg "${CONTROL}")"

  printf "version\trep\tcommand\tstatus\tseconds_wall\tmax_rss_bytes\tlog\n" > "${RUNS_TSV}"

  for rep in $(seq 1 "${REPEATS}"); do
    outdir_q="$(quote_arg "${RESULTS_DIR}/macs1_run_${rep}")"
    cmd="cd ${outdir_q} && ${MACS1_CMD} -t ${treatment_q} -c ${control_q} -f BED -g hs -n macs1_ctcf_5m"
    run_one "macs1" "${rep}" "${cmd}" || true
  done

  for rep in $(seq 1 "${REPEATS}"); do
    outdir_q="$(quote_arg "${RESULTS_DIR}/macs2_run_${rep}")"
    cmd="${MACS2_CMD} callpeak -t ${treatment_q} -c ${control_q} -f BED -g hs -n macs2_ctcf_5m --outdir ${outdir_q} -q 0.01"
    run_one "macs2" "${rep}" "${cmd}" || true
  done

  for rep in $(seq 1 "${REPEATS}"); do
    outdir_q="$(quote_arg "${RESULTS_DIR}/macs3_run_${rep}")"
    cmd="${MACS3_CMD} callpeak -t ${treatment_q} -c ${control_q} -f BED -g hs -n macs3_ctcf_5m --outdir ${outdir_q} -q 0.01"
    run_one "macs3" "${rep}" "${cmd}" || true
  done

  python3 "${SCRIPT_DIR}/summarize_macs_version_survey.py" "${ROOT}"
  echo "Wrote ${RUNS_TSV}, ${RESULTS_DIR}/summary.tsv, and ${RESULTS_DIR}/output_peak_counts.tsv"
}

main() {
  prepare_dirs
  prepare_repo
  ensure_worktree "${MACS1_TREE}" "${MACS1_REF}"
  ensure_worktree "${MACS2_TREE}" "${MACS2_REF}"
  ensure_worktree "${MACS3_TREE}" "${MACS3_REF}"
  record_revisions
  copy_data
  record_versions
  run_benchmarks
}

main "$@"
