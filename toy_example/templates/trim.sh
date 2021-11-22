echo ${reads} > ${sampleId}.txt

trimmomatic \
PE \
/home/guglib/rnaseqs/PE/${reads} \
/home/guglib/rnaseqs/PE/${reads} \
/home/guglib/test/output_forward_paired.fq.gz \
/home/guglib/test/output_forward_unpaired.fq.gz \
/home/guglib/test/output_reverse_paired.fq.gz \
/home/guglib/test/output_reverse_unpaired.fq.gz \
LEADING:${t_lead} \
TRAILING:${t_trail} \
MINLEN:${min_len}
