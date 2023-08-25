#!/usr/bin/env bash
#bash runmicromegas.sh && echo SPheno && read && bash runSpheno.sh
if [ ! "$1" ];then
    echo "Usage: $0 N[\\\\/M]"
    exit
fi
d=$1
if [ "$(echo $d | grep '/')" ]; then
    if [ ! "$(echo $d | grep '\\/')" ]; then
	echo "ERROR! use \\\\/ for /"
	exit
    fi
fi
sed -i 's/^\s*dRInput\s*\=\s*.*;/dRInput = '$d';/' SARAH/Models/Zprime/config.m
echo "Compiling the model for dRInput=$d ..."
bash diffkkk.sh
grep -B1 -A1 "{Xq,Xl" kkk
echo micromegas
#read
bash runmicromegas.sh 
echo SPheno 
#read
bash runSpheno.sh
#TODO: Rescale gp to get the correct relic density
mkdir -p ZprimeOut
df=.$(echo $d | sed 's/\\\//_/')
./SPheno/bin/SPhenoZprime LesHouches.in.Zprime
./micromegas/Zprime/CalcOmega_with_DDetection_MOv5 SPheno.spc.Zprime > ZprimeOut/micromegas$df.out
mv SPheno.spc.Zprime ZprimeOut/SPheno.spc.Zprime$df
#cat micromegas.54_43.out | grep -A1000  ^Xf |  grep -E "^\s*\w" | grep -Ev '^\s*[0-9]' | grep -B1000 -A2 "cross sections"
