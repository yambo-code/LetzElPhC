#!/bin/bash

## Usage sh ph_collect.sh folder_prefix
PREFIX=$1
# Assigns the first command-line argument to the PREFIX variable

if [ -z "$PREFIX" ]; then
    echo "Usage: $0 PREFIX"
    exit 1
fi
# Checks if PREFIX is empty and exits with an error message if no argument is provided

N=$(find . -maxdepth 1 -type d -name "${PREFIX}[0-9]*" | wc -l)
# Finds and counts all directories in the current path that match the prefix followed by a number

echo "Found $N folders with prefix ${PREFIX}"
# Prints the total number of matching folders found

cp -r "${PREFIX}1/_ph0" .
# Copies the base _ph0 directory from the first folder into the current working directory

for (( X=2; X<=N; X++ )); do
    cp -r "${PREFIX}${X}/_ph0/"*.q_${X} _ph0/
done
# Loops from 2 to N and copies the respective q_X folders to the local _ph0 directory

cp "${PREFIX}"*/_ph0/*.phsave/dynmat.*.xml _ph0/*.phsave/
# Copies all dynmat XML files from all source phsave directories into the local phsave directory

for (( X=1; X<=N; X++ )); do
    cp "${PREFIX}${X}/_ph0/"*.phsave/patterns.${X}.xml _ph0/*.phsave/
done
# Loops from 1 to N to copy each patterns.X.xml file into the local phsave directory

for (( X=1; X<=N; X++ )); do
    for f in "${PREFIX}${X}"/*.dyn*; do
        if [[ -e "$f" && "$f" != *.dyn0 ]]; then
            cp "$f" .
        fi
    done
done
# Iterates through all folders and copies all .dyn files to the current directory, excluding .dyn0 files

cp "${PREFIX}1"/*.dyn0 .
# Copies the .dyn0 file specifically from the first folder into the current directory
