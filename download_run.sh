#!/bin/bash 

# Download the data
git clone --recursive -b feature/pathwl https://github.com/marek094/wlcpp.git

# Run the program

cd wlcpp

mkdir build

cd build

apt install -y libboost-all-dev
apt install -y libopenmpi-dev
apt install -y cmake curl


for i in {15..10}
do 
    apt-get install -y g++-$i gcc-$i
    if [ $? -eq 0 ]
    then
        export cxx=g++-$i cc=gcc-$i
        break
    fi
done


CXX=$cxx CC=$cc cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4 main 
cd ..


name=${name:="graph9c.g6"}
cpus=${cpus:=12}
urladdr=https://users.cecs.anu.edu.au/~bdm/data/$name
mkdir graph_data
cd graph_data
curl -O $urladdr
cd ..


./build/src/main graph_data/$name $cpus