#!/bin/bash

VERSION=`cat VERSION.txt`

docker push trinityctat/ctat_mutations_benchmark:${VERSION} 
docker push trinityctat/ctat_mutations_benchmark:latest 
