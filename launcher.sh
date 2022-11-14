#!/bin/bash

InputFile=$1
OutputFile=${InputFile%%.inp}.out
ScratchFile=${InputFile%%.inp}.integrals
QuoteFile=$((1 + $RANDOM % 45))

cat /Users/vcastor/Documents/GitHub/Nika/writer/main_to_write/header > ${OutputFile}

/Users/vcastor/Documents/GitHub/Nika/Nika.exe $1

cat ./tmp/out.out >> ${OutputFile}

echo "*******UNIQUE OVERLAP MATRIX VALUES*******" >> ${ScratchFile}
cat ./tmp/Overlap.int >> ${OutputFile}

echo "*******UNIQUE KINETIC MATRIX VALUES*******" >> ${ScratchFile}
cat ./tmp/Kinetic.int >> ${OutputFile}

echo "******UNIQUE POTENTIAL MATRIX VALUES******" >> ${ScratchFile}
cat ./tmp/Potential.int >> ${OutputFile}

echo "*****UNIQUE TWO ELECTRON MATRIX VALUES****" >> ${ScratchFile}
cat ./tmp/TwoElectron.int >> ${OutputFile}

cat ./writer/succesfull_quotes/${QuoteFile} >> ${OutputFile}

rm tmp.1 ./tmp/out.out ./tmp/Overlap.int ./tmp/Kinetic.int ./tmp/Potential.int ./tmp/TwoElectron.int
