#!/usr/bin/env bash

_CLONE="git clone /root/git@github.com:ecatanzani/eFlux.git"
_BUILD="cd /root/eFlux && make rebuild"

${_CLONE}
${_BUILD}