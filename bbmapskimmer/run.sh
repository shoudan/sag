#!/bin/bash

set -e

LIBFQ=`cat ./libs_all | sed -n ${SGE_TASK_ID}p`

/global/homes/a/andreopo/memtime /global/projectb/sandbox/rqc/sliang/sag/bbmapskimmer/mapone.sh $LIBFQ
