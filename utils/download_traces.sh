#!/bin/bash

# MIT License
#
# Copyright (c) 2024 David Schall and EASE lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Execute this script using
#   ./eval_all.sh

set -x -e



## Parse the TRACES.list file
## TRACES can be commented out with '#'
TRACES=""


TRACES="${TRACES} mwnginxfpm-wiki"
TRACES="${TRACES} dacapo-kafka"
TRACES="${TRACES} dacapo-tomcat"
TRACES="${TRACES} dacapo-spring"
TRACES="${TRACES} renaissance-finagle-chirper"
TRACES="${TRACES} renaissance-finagle-http"
TRACES="${TRACES} benchbase-tpcc"
TRACES="${TRACES} benchbase-twitter"
TRACES="${TRACES} benchbase-wikipedia"
TRACES="${TRACES} nodeapp-nodeapp"
TRACES="${TRACES} charlie.1006518"
TRACES="${TRACES} delta.507252"
TRACES="${TRACES} merced.467915"
TRACES="${TRACES} whiskey.426708"



TRACE_DIR="./traces"


for fn in $TRACES; do

    wget -O $TRACE_DIR/$fn.champsim.trace.gz https://zenodo.org/records/13133243/files/$fn.champsim.trace.gz?download=1

done
