#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

function exitOnError {

    ERR=$?

    if [ "$ERR" != 0 ]; then

        echo
        echo FATAL ERROR
        echo $ERR
        echo
        exit $ERR

    fi

}

if [ ! -e "README.txt" ]; then
    echo
    echo "    Please create README.txt file for #REQUEST"
    echo
    exit
fi

PROJNO=$(basename $PWD | sed 's/Proj_//')

PTYPE=$(basename $(dirname $PWD))

if [ "$PTYPE" == "variant" ]; then
    echo "Creating" $PROJNO "as a" $PTYPE "project"

    python3 ~/Code/LIMS/LimsETL/getProjectFiles.py $PROJNO
    exitOnError

    Rscript --no-save ~/Code/LIMS/LimsETL/makeVariantProject.R *_metadata.yaml
    exitOnError

    SPECIES=$(cat _request  | fgrep Species | sed 's/.*: //')

    echo "SPECIES =" $SPECIES

    if [ "$SPECIES" == "Mouse" ]; then

        $SDIR/findMousePools.R *_sample_mapping.txt
        exitOnError
        cat mapping.pool >> *_sample_mapping.txt
    fi

    Rscript --no-save $SDIR/generateMousePairGroupFiles.R
    exitOnError

    $SDIR/makeProjectFiles.R
    exitOnError

    echo
    echo "    Done - Make sure to check Pairing file"
    echo "    Fix if necessary and then run"
    echo "    - validateProjectFiles.sh"
    echo

    exit

fi

echo "Unknown PTYPE =" $PTYPE

