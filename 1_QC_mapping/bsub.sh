
for i in m6a_IP_A1 m6a_IP_A2 m6a_IP_A3 m6a_IP_S1 m6a_IP_S2 m6a_IP_S3   m6a_input2_A1 m6a_input2_A2 m6a_input2_A3 m6a_input2_S1 m6a_input2_S2 m6a_input2_S3 ; 
do
bsub -o o.merip_seq_${i} -e e.merip_seq_${i} -sp 100 -G $TEAM -q long -n8 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]" -M16000 sh merip_seq.sh ${i}
done
