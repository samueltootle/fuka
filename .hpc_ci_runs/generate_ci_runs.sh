#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_SCRIPT="$script_dir/env_template.sh"
SLURM_TEMPLATE="$script_dir/slurm_submission_template.sh"

FUKA_SOLVER_DIR="${HOME_KADATH}/codes/FUKA"

# 1. Ensure all solvers are built and up to date

# . $HOME_KADATH/.local_dev_scripts/build_fuka_solvers.sh

# We only have SLURM for now
source "$SLURM_TEMPLATE"

update_root() {
    local root="EMAIL"
    local param="$2"
    local val="$3"
    local filename="$4"
    sed -E -i "1,$ {0,/$root/! { /\{/,/\}/ s/($param) (.*)/\1 $val/}}"  $filename
}

update_bin_param() {
    local param="EMAIL"
    local val="$2"
    local co=$3
    local filename="$4"
    sed -E -i "1,$ {0,/binary/! {/$co/,/\}/ s/($param) (.*)/\1 $val/}}" $filename
}

generate_bh_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BH/bh_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BH/bh_first_run/

    $FUKA_SOLVER_DIR/BH/bin/Release/solve

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BH" \
        "initial_bh.info" \
        "bh_first_run" \
        "200" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BH/bh_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BH/bh_second_run/

    $FUKA_SOLVER_DIR/BH/bin/Release/solve

    update_root "bh" "chi" "0.85" "initial_bh.info"
    update_root "bh" "mch" "1.0" "initial_bh.info"
    update_root "bh" "res" "11" "initial_bh.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BH" \
        "initial_bh.info" \
        "bh_second_run" \
        "200" \
        "24:00:00" \
        "EMAIL"
}

generate_ns_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS/ns_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS/ns_first_run/

    $FUKA_SOLVER_DIR/NS/bin/Release/solve

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS" \
        "initial_ns.info" \
        "ns_first_run" \
        "200" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS/ns_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS/ns_second_run/

    $FUKA_SOLVER_DIR/NS/bin/Release/solve

    update_root "ns" "chi" "0.6" "initial_ns.info"
    update_root "ns" "madm" "2.3" "initial_ns.info"
    update_root "ns" "res" "11" "initial_ns.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS" \
        "initial_ns.info" \
        "ns_second_run" \
        "200" \
        "24:00:00" \
        "EMAIL"
}

generate_bbh_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BBH/bbh_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BBH/bbh_first_run/

    $FUKA_SOLVER_DIR/BBH/bin/Release/solve

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BBH" \
        "initial_bbh.info" \
        "bbh_first_run" \
        "400" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BBH/bbh_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BBH/bbh_second_run/

    $FUKA_SOLVER_DIR/BBH/bin/Release/solve

    update_root "binary" "res" "11" "initial_bbh.info"
    update_bin_param "res" "9" "bh1" "initial_bbh.info"
    update_bin_param "mch" "0.1" "bh1" "initial_bbh.info"
    update_bin_param "chi" "-0.5" "bh1" "initial_bbh.info"

    update_bin_param "res" "9" "bh2" "initial_bbh.info"
    update_bin_param "mch" "0.9" "bh2" "initial_bbh.info"
    update_bin_param "chi" "0.85" "bh2" "initial_bbh.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BBH" \
        "initial_bbh.info" \
        "bbh_second_run" \
        "400" \
        "24:00:00" \
        "EMAIL"
}

generate_bns_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BNS/bns_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BNS/bns_first_run/

    $FUKA_SOLVER_DIR/BNS/bin/Release/solve

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BNS" \
        "initial_bns.info" \
        "bns_first_run" \
        "400" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BNS/bns_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BNS/bns_second_run/

    $FUKA_SOLVER_DIR/BNS/bin/Release/solve

    update_root "binary" "res" "11" "initial_bns.info"

    update_bin_param "chi" "0.0" "ns1" "initial_bns.info"
    update_bin_param "madm" "1.18" "ns1" "initial_bns.info"
    update_bin_param "res" "11" "ns1" "initial_bns.info"

    update_bin_param "chi" "0.52" "ns2" "initial_bns.info"
    update_bin_param "madm" "2.42" "ns2" "initial_bns.info"
    update_bin_param "res" "11" "ns2" "initial_bns.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BNS" \
        "initial_bns.info" \
        "bns_second_run" \
        "400" \
        "24:00:00" \
        "EMAIL"
}

generate_bhns_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BHNS/bhns_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BHNS/bhns_first_run/

    $FUKA_SOLVER_DIR/BHNS/bin/Release/solve

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BHNS" \
        "initial_bhns.info" \
        "bhns_first_run" \
        "400" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/BHNS/bhns_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/BHNS/bhns_second_run/

    $FUKA_SOLVER_DIR/BHNS/bin/Release/solve

    update_root "binary" "res" "11" "initial_bhns.info"

    update_bin_param "chi" "0.0" "ns1" "initial_bhns.info"
    update_bin_param "madm" "1.18" "ns1" "initial_bhns.info"
    update_bin_param "res" "11" "ns1" "initial_bhns.info"

    update_bin_param "chi" "0.52" "ns2" "initial_bhns.info"
    update_bin_param "mch" "2.42" "ns2" "initial_bhns.info"
    update_bin_param "res" "9" "ns2" "initial_bhns.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/BHNS" \
        "initial_bhns.info" \
        "bhns_second_run" \
        "400" \
        "24:00:00" \
        "EMAIL"
}

generate_ns_diffrot_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS_DIFFROT/ns_diffrot_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS_DIFFROT/ns_diffrot_first_run/

    $FUKA_SOLVER_DIR/NS_DIFFROT/bin/Release/solve

    update_root "ns" "chi" "0.1" "initial_ns.info"
    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS_DIFFROT" \
        "initial_ns.info" \
        "ns_diffrot_first_run" \
        "200" \
        "24:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS_DIFFROT/ns_diffrot_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS_DIFFROT/ns_diffrot_second_run/

    $FUKA_SOLVER_DIR/NS_DIFFROT/bin/Release/solve

    update_root "ns" "chi" "0.1" "initial_ns.info"
    update_root "ns" "madm" "2.3" "initial_ns.info"
    update_root "ns" "res" "13" "initial_ns.info"

    update_root "ns" "A_ratio" "1" "initial_ns.info"
    update_root "ns" "R_ratio" "0.6" "initial_ns.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS_DIFFROT" \
        "initial_ns.info" \
        "ns_diffrot_second_run" \
        "200" \
        "24:00:00" \
        "EMAIL"
}

generate_ns_isotropic_ci_runs() {
    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS_isotropic/ns_isotropic_first_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS_isotropic/ns_isotropic_first_run/

    $FUKA_SOLVER_DIR/NS_isotropic/bin/Release/solve

    update_root "ns" "chi" "0.1" "initial_2dns.info"
    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS_isotropic" \
        "initial_2dns.info" \
        "ns_isotropic_first_run" \
        "20" \
        "01:00:00" \
        "EMAIL"

    cd ..

    mkdir -p $HOME_KADATH/.hpc_ci_runs/CIs/NS_isotropic/ns_isotropic_second_run/
    cd $HOME_KADATH/.hpc_ci_runs/CIs/NS_isotropic/ns_isotropic_second_run/

    $FUKA_SOLVER_DIR/NS_isotropic/bin/Release/solve

    update_root "ns" "chi" "0.1" "initial_2dns.info"
    update_root "ns" "madm" "2.3" "initial_2dns.info"
    update_root "ns" "res" "13" "initial_2dns.info"

    update_root "ns" "A_ratio" "1" "initial_2dns.info"
    update_root "ns" "R_ratio" "0.6" "initial_2dns.info"

    update_root "stages" "differential_rotation" "on" "initial_2dns.info"

    generate_submission \
        `pwd` \
        "$FUKA_SOLVER_DIR/NS_isotropic" \
        "initial_2dns.info" \
        "ns_isotropic_second_run" \
        "20" \
        "01:00:00" \
        "EMAIL"
}

EMAIL=$1
generate_bh_ci_runs EMAIL
generate_ns_ci_runs EMAIL
generate_ns_diffrot_ci_runs EMAIL
generate_ns_isotropic_ci_runs EMAIL
generate_bbh_ci_runs EMAIL
generate_bns_ci_runs EMAIL
generate_bhns_ci_runs EMAIL
