#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t trinityctat/ctat_mutations_benchmark:${VERSION} .
docker build -t trinityctat/ctat_mutations_benchmark:latest .
