#!/usr/bin/env bash

sed 's/\(from loprop.*\)/#\1/g' *py
sed 's/\(import h5py*\)/#\1/g' *py
sed 's/\(from gaussian*\)/#\1/g' *py
sed 's/\(import gaussian*\)/#\1/g' *py
