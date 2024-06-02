#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t trinityctat/ctat_mutations_benchmarking:${VERSION} .
docker build -t trinityctat/ctat_mutations_benchmarking:latest .
