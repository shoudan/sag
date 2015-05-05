#!/bin/bash

set -e

LIBFQ=`cat ./run.config | sed -n ${SGE_TASK_ID}p`

/global/homes/a/andreopo/memtime ./mapone.sh $LIBFQ
