#!/bin/sh

if [ ! -e logs ]
then
mkdir logs
fi

wdir=`pwd`
pdir=`dirname $wdir`

while read lib; do
while read ref; do
echo $lib $ref
done < ${wdir}/libs_all
done < ${wdir}/libs_all
