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

#The source dir is one step up
DIR=$DIR/..


#Uncomment private source repo dependencies
sed -i 's/^\(from loprop.*\)/#\1/g' $DIR/*.py
sed -i 's/^\(import h5py*\)/#\1/g' $DIR/*.py
sed -i 's/^\(from gaussian*\)/#\1/g' $DIR/*.py
sed -i 's/^\(import gaussian*\)/#\1/g' $DIR/*.py

