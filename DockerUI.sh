#!/usr/bin/env bash
ip=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
echo "X11 Docker forwarding IP: $ip"
xhost + $ip
docker run --rm -it -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$ip:0 -v /Users/enrico/Repo:/mnt ecatanzani/dampesw:centos bash
