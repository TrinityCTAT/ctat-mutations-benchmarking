#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t trinityrnaseq/ctat_mutations_benchmark:${VERSION} .
docker build -t trinityrnaseq/ctat_mutations_benchmark:latest .
