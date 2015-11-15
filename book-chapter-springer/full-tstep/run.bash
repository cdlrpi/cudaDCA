#!/bin/bash
cd parallel
make
./runSim
cd ..
cd serial
make
./runSim
