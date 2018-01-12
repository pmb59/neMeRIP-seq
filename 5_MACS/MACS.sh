OUT="../m6A/MACS/"

#Activin
ip1="../m6A/tophat_m6a_IP_A1/m6a_IP_A1.bam"
ip2="../m6A/tophat_m6a_IP_A2/m6a_IP_A2.bam"
ip3="../m6A/tophat_m6a_IP_A3/m6a_IP_A3.bam"
input="../m6A/tophat_m6a_input2_A/m6a_input2_A.bam"

macs2 callpeak -t $ip1 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n A1
macs2 callpeak -t $ip2 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n A2
macs2 callpeak -t $ip3 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n A3

#SB
ip1="../m6A/tophat_m6a_IP_S1/m6a_IP_S1.bam"
ip2="../m6A/tophat_m6a_IP_S2/m6a_IP_S2.bam"
ip3="../m6A/tophat_m6a_IP_S3/m6a_IP_S3.bam"
input="../m6A/tophat_m6a_input2_S/m6a_input2_S.bam"

macs2 callpeak -t $ip1 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n S1
macs2 callpeak -t $ip2 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n S2
macs2 callpeak -t $ip3 -c $input -f BAM --outdir $OUT -q 0.001 --nomodel --keep-dup all -n S3





