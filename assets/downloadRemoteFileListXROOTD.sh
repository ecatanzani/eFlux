#!/usr/bin/env bash

while read file ; do
    xrdcp $file $2
done <$1
