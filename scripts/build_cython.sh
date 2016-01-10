#!/bin/bash

#Get location of dir of this script

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    [[ $SOURCE != /* ]] && SOURCE="$SCRIPT_DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

echo $DIR/..

PDDIR=$DIR/../src/pd
#if pd exists
if [[ -d $PDDIR ]] ;
then
	cd $PDDIR
	python setup.py build_ext --inplace
	cd ..
fi

