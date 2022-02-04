#!/bin/bash

wget https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz
tar -xzf v4.2.2.tar.gz 
rm v4.2.2.tar.gz 
cd shapeit4-4.2.2/
make all
cd maps/
tar -xzf genetic_maps.b38.tar.gz
