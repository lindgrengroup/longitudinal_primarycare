#!/bin/bash

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb
do
	zcat BOLT_results/${STRATA}_assoc.stats.gz | head -n 1 > replication_results/${STRATA}_replication.txt
	zgrep -wFf snps_to_replicate.txt BOLT_results/${STRATA}_assoc.stats.gz >> replication_results/${STRATA}_replication.txt
done
