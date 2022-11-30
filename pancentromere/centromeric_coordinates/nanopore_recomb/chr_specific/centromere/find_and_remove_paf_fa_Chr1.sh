#!/bin/bash

for CHROM in $(seq 1);
do
  echo Chr${CHROM};
  echo Removing Col-0 segment PAFs and FASTAs;
  find Chr${CHROM}/segments/Col-0/co/ -mindepth 1 -maxdepth 1 -type f -name "*.*a*" -exec rm {} + ;
  find Chr${CHROM}/segments/Col-0/nco/ -mindepth 1 -maxdepth 1 -type f -name "*.*a*" -exec rm {} + ;
  echo Removing Ler-0 segment PAFs and FASTAs;
  find Chr${CHROM}/segments/Ler-0/co/ -mindepth 1 -maxdepth 1 -type f -name "*.*a*" -exec rm {} + ;
  find Chr${CHROM}/segments/Ler-0/nco/ -mindepth 1 -maxdepth 1 -type f -name "*.*a*" -exec rm {} + ;
done 
