#!/usr/bin/env bash

docker build -t docker-build:latest docker-build;
#docker build --no-cache -t docker-build:latest docker-build;
docker run -v $PWD:/home/docker/BuddySuite docker-build:latest;
