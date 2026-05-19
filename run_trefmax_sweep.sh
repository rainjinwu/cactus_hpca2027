#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/wurunjin/cactus-ae"
SIM="${ROOT}/memsim/sim_dramsim3"
TRACE="${ROOT}/traces/omnetpp_17.mtf.gz"
CONFIG_DIR="${ROOT}/DRAMsim3/configs/cactus"
RESULT_DIR="${ROOT}/RESULTS_TREFMAX_SWEEP_RATE8_250M"
INST_LIMIT="${INST_LIMIT:-250000000}"
RATEMODE="${RATEMODE:-8}"
JOBS="${JOBS:-4}"

mkdir -p "${RESULT_DIR}" "${CONFIG_DIR}/sweep"

make_config() {
    local policy="$1"
    local interval="$2"
    local src="${CONFIG_DIR}/DDR5_32Gb_mop4_cactus_fig15_${policy}_trefmax_trhd250.ini"
    local dst="${CONFIG_DIR}/sweep/DDR5_32Gb_mop4_cactus_fig15_${policy}_trefmax${interval}_trhd250.ini"

    cp "${src}" "${dst}"
    perl -0pi -e "s/tref_interval = [0-9]+/tref_interval = ${interval}/" "${dst}"
    echo "${dst}"
}

run_one() {
    local policy="$1"
    local interval="$2"
    local cfg
    cfg="$(make_config "${policy}" "${interval}")"
    local out="${RESULT_DIR}/mop4_cactus_fig15_${policy}_trefmax${interval}_trhd250_omnetpp_17.mtf.gz.out"
    local tmp="${out}.tmp"

    echo "[run] policy=${policy} tref_interval=${interval} -> ${out}"
    "${SIM}" -inst_limit "${INST_LIMIT}" -ratemode "${RATEMODE}" \
        -dramsim3cfg "${cfg}" "${TRACE}" > "${tmp}" 2>&1
    LC_ALL=C tr -d '\000' < "${tmp}" > "${tmp}.clean"
    mv -f "${tmp}.clean" "${tmp}"
    mv -f "${tmp}" "${out}"
}

wait_for_slot() {
    while [ "$(jobs -rp | wc -l)" -ge "${JOBS}" ]; do
        sleep 5
    done
}

for interval in 1 2 4 8; do
    for policy in direct random; do
        wait_for_slot
        run_one "${policy}" "${interval}" &
    done
done
wait

csv="${RESULT_DIR}/trefmax_sweep_omnetpp_trhd250.csv"
{
    echo "policy,tref_interval,avg_ipc,sys_cycles,alerts_ch0,alerts_ch1,rfmab_ch0,rfmab_ch1,tref_max_ch0,tref_max_ch1,tref_skip_ch0,tref_skip_ch1,refab_ch0,refab_ch1"
    for interval in 1 2 4 8; do
        for policy in direct random; do
            out="${RESULT_DIR}/mop4_cactus_fig15_${policy}_trefmax${interval}_trhd250_omnetpp_17.mtf.gz.out"
            avg_ipc="$(awk '/CORE_[0-9][0-9]_IPC/ {sum += $3; n++} END {if (n) printf "%.6f", sum / n; else printf "NA"}' "${out}")"
            sys_cycles="$(awk '/SYS_CYCLES/ {print $3}' "${out}")"
            alerts="$(awk '/num_alerts/ {printf "%s%s", sep, $3; sep="/"}' "${out}")"
            rfmab="$(awk '/num_rfmab_cmds/ {printf "%s%s", sep, $3; sep="/"}' "${out}")"
            tref_max="$(awk '/num_tref_max_mitigs/ {printf "%s%s", sep, $3; sep="/"}' "${out}")"
            tref_skip="$(awk '/num_tref_rfmab_skips/ {printf "%s%s", sep, $3; sep="/"}' "${out}")"
            refab="$(awk '/num_refab_cmds/ {printf "%s%s", sep, $3; sep="/"}' "${out}")"
            echo "${policy},${interval},${avg_ipc},${sys_cycles},${alerts/\//,},${rfmab/\//,},${tref_max/\//,},${tref_skip/\//,},${refab/\//,}"
        done
    done
} > "${csv}"

echo "[done] ${csv}"
