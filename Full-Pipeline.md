# Full pipeline
Using the CASP13 target T0957s2-D1 as an example.
<hr>  

### Step1. Make alignments
```bash
/ssdA/common-tools/DeepMSA/hhsuite2/scripts/build_MSA_e0.01.py ./T0957s2-D1.fasta -outdir=./T0957s2-D1/ -hhblitsdb=/ssdA/common-tools/uniclust30_2018_08_hhsuite/uniclust30_2018_08 -jackhmmerdb=/ssdA/common-tools/uniref_2019_10_24/uniref90.fasta -hmmsearchdb=/ssdA/common-tools/fasta_188GB_plass_soedinglab/SRC.fasta:/ssdA/common-tools/metaclust-2018-Jun-22/metaclust_nr.fasta
```
