#!/bin/bash
rm *~

for prod in PID_0.5_to_10_GeV_gamma_neutron_pi_LUXE
do
    cp CheckInteraction.C $prod/.
    cd $prod
    echo $prod
    for sample in PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_gamma_0.5to10GeV PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_neutron_0.5to10GeV PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_pi-_0.5to10GeV
    do
	root -q CheckInteraction.C\(\"$sample\"\) > ${sample}_events.out
    done
    rm CheckInteraction.C
    cd -
done
