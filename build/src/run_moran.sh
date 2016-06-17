#!/usr/bin/env bash

function run_moran_model {
    pop=$1
    tm=$2
    tm_no_scientific=$(printf '%.0f' $tm)
    mut=$3
    fn="moran_model_n${pop}_t${tm}_m${mut}.txt"
    ./moran_model_boundary_mutation -n $pop -t $tm_no_scientific \
                                    -m $mut -f $fn
}
export -f run_moran_model

parallel "run_moran_model {1} {2} {3}" \
         ::: 0600 1000 \
         ::: 1e7 \
         ::: 0.10
