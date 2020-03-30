#!/usr/bin/env bash
while :
do
	_TMPDIR=$RANDOM
	_DIR=/tmp/dampe_builder/$_TMPDIR
	echo "Current WD: $_DIR"
	if [ ! -d $_DIR ]; then
		mkdir $_DIR
		break
	fi
done
cd $_DIR
git clone --recurse-submodules git@github.com:ecatanzani/eFlux.git
docker run --rm -v $(pwd):/home dampe_builder:centos
