#!/bin/bash

read -p "Number of qubits: " n
read -p "Pauli observable:[X,Y,Z] " pauli
read -p "Support:[half,full,treshold] " supp
read -p "Number of samples" m
read -p "Number of batches: " l
read -p "Frame ensemble:[BW/LC] " fr

mkdir -p ./${n}_qubits/logs

count=0
j=1
while [ "$count" -lt "$l" ]
do
	if [ ! -f "./${n}_qubits/data/${fr}_${supp}_${pauli}_${m}_${j}.npy" ]
	then
		let "count=count+1"
		nohup python3 ./pauli_estimation.py -n $n --supp $supp --outfile ${m}_${j} --samples $m --pauli $pauli --frame $fr &> ./${n}_qubits/logs/${fr}_${pauli}_${supp}_${m}_log${j}.txt &
	fi
	let "j=j+1"
done
