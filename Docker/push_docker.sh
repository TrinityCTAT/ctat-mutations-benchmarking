#!/bin/bash

VERSION=`cat VERSION.txt`

docker push trinityctat/ctat_mutations_benchmarking:${VERSION} 
docker push trinityctat/ctat_mutations_benchmarking:latest 
