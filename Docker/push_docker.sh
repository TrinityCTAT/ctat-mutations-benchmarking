#!/bin/bash

VERSION=`cat VERSION.txt`

docker push trinityrnaseq/ctat_mutations_benchmark:${VERSION} 
docker push trinityrnaseq/ctat_mutations_benchmark:latest 
