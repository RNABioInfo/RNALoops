#!/bin/bash
#literally just a bash script to avoid having to type python3 before every run of RNALoops
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python3 "$SCRIPT_DIR/src/RNALoops.py" $@