#!/usr/bin/env bash
# regression_tests/lib/regression_checks.sh
#
# Shared check functions for MadVoro THUNDER regression tests.

REGRESSION_CHECK_MSG=""

check_no_fatal_markers() {
    local stdout_log="$1"
    local stderr_log="$2"
    local fatal_patterns="Segmentation fault|SIGSEGV|std::bad_alloc|terminate called|Aborted|MPI_ABORT|MadVoroException"

    for log in "$stdout_log" "$stderr_log"; do
        [[ -f "$log" ]] || continue
        if grep -qE "${fatal_patterns}" "$log" 2>/dev/null; then
            local match
            match="$(grep -m1 -E "${fatal_patterns}" "$log")"
            REGRESSION_CHECK_MSG="Fatal marker in $(basename "$log"): ${match}"
            return 1
        fi
    done
    return 0
}

is_nonempty_and_newer() {
    local filepath="$1"
    local start_epoch="$2"
    if [[ ! -f "$filepath" ]]; then
        REGRESSION_CHECK_MSG="Missing output file: $filepath"
        return 1
    fi
    if [[ ! -s "$filepath" ]]; then
        REGRESSION_CHECK_MSG="Empty output file: $filepath"
        return 1
    fi
    if [[ "$start_epoch" != "0" ]]; then
        local mtime
        mtime="$(stat -c %Y "$filepath" 2>/dev/null || stat -f %m "$filepath" 2>/dev/null)" || true
        if [[ -n "$mtime" && "$mtime" -lt "$start_epoch" ]]; then
            REGRESSION_CHECK_MSG="Stale output file (mtime before suite start): $filepath"
            return 1
        fi
    fi
    return 0
}

check_voronoi_volume_case() {
    local run_dir="$1"
    local start_epoch="$2"
    local stdout_log="$3"
    local stderr_log="$4"
    local metrics_file="${run_dir}/voronoi_volume_metrics.txt"
    local rel_error
    local pass_flag
    local max_rel_error="${VORONOI_VOLUME_MAX_REL_ERROR:-1e-10}"

    check_no_fatal_markers "$stdout_log" "$stderr_log" || return 1
    is_nonempty_and_newer "$metrics_file" "$start_epoch" || return 1

    rel_error=$(awk '$1 == "rel_error" { print $2 }' "$metrics_file")
    pass_flag=$(awk '$1 == "pass" { print $2 }' "$metrics_file")

    if [[ -z "$rel_error" || -z "$pass_flag" ]]; then
        REGRESSION_CHECK_MSG="failed to parse voronoi volume metrics"
        return 1
    fi
    if [[ "$pass_flag" != "1" ]]; then
        REGRESSION_CHECK_MSG="voronoi_volume test reported pass=${pass_flag}"
        return 1
    fi
    if ! awk -v r="$rel_error" -v t="$max_rel_error" 'BEGIN { exit !(r < t) }'; then
        REGRESSION_CHECK_MSG="voronoi_volume rel_error exceeds threshold (${rel_error} >= ${max_rel_error})"
        return 1
    fi

    REGRESSION_CHECK_MSG="PASS (rel_error=${rel_error})"
    return 0
}

check_voronoi_volume_periodic_case() {
    local run_dir="$1"
    local start_epoch="$2"
    local stdout_log="$3"
    local stderr_log="$4"
    local metrics_file="${run_dir}/voronoi_volume_periodic_metrics.txt"
    local rel_error
    local pass_flag
    local periodic_flag
    local max_rel_error="${VORONOI_VOLUME_MAX_REL_ERROR:-1e-10}"

    check_no_fatal_markers "$stdout_log" "$stderr_log" || return 1
    is_nonempty_and_newer "$metrics_file" "$start_epoch" || return 1

    rel_error=$(awk '$1 == "rel_error" { print $2 }' "$metrics_file")
    pass_flag=$(awk '$1 == "pass" { print $2 }' "$metrics_file")
    periodic_flag=$(awk '$1 == "periodic" { print $2 }' "$metrics_file")

    if [[ -z "$rel_error" || -z "$pass_flag" || -z "$periodic_flag" ]]; then
        REGRESSION_CHECK_MSG="failed to parse periodic voronoi volume metrics"
        return 1
    fi
    if [[ "$periodic_flag" != "1" ]]; then
        REGRESSION_CHECK_MSG="voronoi_volume_periodic metrics must have periodic=1"
        return 1
    fi
    if [[ "$pass_flag" != "1" ]]; then
        REGRESSION_CHECK_MSG="voronoi_volume_periodic test reported pass=${pass_flag}"
        return 1
    fi
    if ! awk -v r="$rel_error" -v t="$max_rel_error" 'BEGIN { exit !(r < t) }'; then
        REGRESSION_CHECK_MSG="voronoi_volume_periodic rel_error exceeds threshold (${rel_error} >= ${max_rel_error})"
        return 1
    fi

    REGRESSION_CHECK_MSG="PASS (periodic rel_error=${rel_error})"
    return 0
}

check_voronoi_parallel_check_case() {
    local run_dir="$1"
    local start_epoch="$2"
    local stdout_log="$3"
    local stderr_log="$4"

    check_no_fatal_markers "$stdout_log" "$stderr_log" || return 1

    if grep -q "voronoi_parallel_check PASS" "$stdout_log" 2>/dev/null; then
        local seed
        seed="$(grep -m1 "voronoi_parallel_check PASS" "$stdout_log" | grep -oP 'seed=\K[0-9]+')"
        REGRESSION_CHECK_MSG="PASS (seed=${seed:-?})"
        return 0
    fi

    if grep -q "voronoi_parallel_check FAIL" "$stdout_log" 2>/dev/null; then
        REGRESSION_CHECK_MSG="$(grep -m1 "voronoi_parallel_check FAIL" "$stdout_log")"
        return 1
    fi

    REGRESSION_CHECK_MSG="No PASS/FAIL marker found for voronoi_parallel_check"
    return 1
}

check_voronoi_mock_mesh_periodic_case() {
    local run_dir="$1"
    local start_epoch="$2"
    local stdout_log="$3"
    local stderr_log="$4"
    local metrics_file="${run_dir}/voronoi_mock_mesh_periodic_metrics.txt"
    local rel_error_before
    local rel_error_after
    local pass_flag
    local did_rebalance
    local max_rel_error="${VORONOI_VOLUME_MAX_REL_ERROR:-1e-10}"

    check_no_fatal_markers "$stdout_log" "$stderr_log" || return 1
    is_nonempty_and_newer "$metrics_file" "$start_epoch" || return 1

    rel_error_before=$(awk '$1 == "rel_error_before" { print $2 }' "$metrics_file")
    rel_error_after=$(awk '$1 == "rel_error_after" { print $2 }' "$metrics_file")
    pass_flag=$(awk '$1 == "pass" { print $2 }' "$metrics_file")
    did_rebalance=$(awk '$1 == "did_rebalance" { print $2 }' "$metrics_file")

    if [[ -z "$rel_error_before" || -z "$rel_error_after" || -z "$pass_flag" || -z "$did_rebalance" ]]; then
        REGRESSION_CHECK_MSG="failed to parse periodic MockMesh metrics"
        return 1
    fi
    if [[ "$did_rebalance" != "1" ]]; then
        REGRESSION_CHECK_MSG="periodic MockMesh test did not rebalance"
        return 1
    fi
    if [[ "$pass_flag" != "1" ]]; then
        REGRESSION_CHECK_MSG="voronoi_mock_mesh_periodic test reported pass=0"
        return 1
    fi
    if ! awk -v r="$rel_error_before" -v t="$max_rel_error" 'BEGIN { exit !(r < t) }'; then
        REGRESSION_CHECK_MSG="volume before MockMesh exceeds threshold (${rel_error_before} >= ${max_rel_error})"
        return 1
    fi
    if ! awk -v r="$rel_error_after" -v t="$max_rel_error" 'BEGIN { exit !(r < t) }'; then
        REGRESSION_CHECK_MSG="volume after MockMesh exceeds threshold (${rel_error_after} >= ${max_rel_error})"
        return 1
    fi
    if ! grep -q "voronoi_mock_mesh_periodic PASS=1" "$stdout_log" 2>/dev/null; then
        REGRESSION_CHECK_MSG="stdout missing voronoi_mock_mesh_periodic PASS=1 marker"
        return 1
    fi

    REGRESSION_CHECK_MSG="Periodic MockMesh rebalance passed (rel_error_after=${rel_error_after})"
    return 0
}
