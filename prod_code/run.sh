#!/bin/bash

set -e

LIBFQ=`cat ./run.config | sed -n ${SGE_TASK_ID}p`

/global/homes/a/andreopo/memtime /global/projectb/sandbox/rqc/sliang/sag/prod_code/mapone.sh $LIBFQ
