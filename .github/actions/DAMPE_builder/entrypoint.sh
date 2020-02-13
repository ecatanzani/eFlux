#!/usr/bin/env bash

_CLONE="git clone git@github.com:ecatanzani/eFlux.git"
_BUILD="cd eFlux && make rebuild"

${_CLONE}
${_BUILD}