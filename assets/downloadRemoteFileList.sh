#!/usr/bin/env bash

_SSH_CONFIG="cnaf_dampe"
while read file ; do
    scp $_SSH_CONFIG:$file $2
done <$1