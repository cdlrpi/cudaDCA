#!/bin/bash
cd run-sim
make
./runSim
cd ..
cd serial
make
./runSim
