#!/bin/bash 

# Download the data
git clone -b feature/pathwl https://github.com/marek094/wlcpp.git

# Run the program

cd wlcpp

mkdir build

cd build

sudo apt-get install libboost-all-dev
sudo apt-get install libopenmpi-dev
sudo apt-get install cmake curl


for i in {15..10}
do 
    sudo apt-get install g++-$i gcc-$i
    if [ $? -eq 0 ]
    then
        export CXX=g++-$i CC=gcc-$i
        break
    fi
done


cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4 main 
cd ..


name=graph9c.g6
urladdr=https://users.cecs.anu.edu.au/~bdm/data/$name
mkdir graph_data
cd graph_data
curl -O $urladdr
cd ..


./build/src/main graph_data/$name 12