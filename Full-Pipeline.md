# Full pipeline
Using the CASP13 target T0957s2-D1 as an example.
<hr>  

### Step1. Make alignments
```bash
/ssdA/common-tools/DeepMSA/hhsuite2/scripts/build_MSA_e0.001.py ./T0957s2-D1.fasta \
-outdir=./T0957s2-D1/ \
-hhblitsdb=/ssdA/common-tools/uniclust30_2018_08_hhsuite/uniclust30_2018_08 \
-jackhmmerdb=/ssdA/common-tools/uniref_2019_10_24/uniref90.fasta \
-hmmsearchdb=/ssdA/common-tools/fasta_188GB_plass_soedinglab/SRC.fasta\
:/ssdA/common-tools/metaclust-2018-Jun-22/metaclust_nr.fasta
```
Repeat for various e-value thresholds
### Step2. Predict distances and contacts
```bash
python3 /home/badri/casp14/35-PrayogRealDistance/src/predict-distance.py \
-w ../35-PrayogRealDistance/job/dist.hdf5.val.hdf5 \
-a ../38-DeepMSAforCASP13FM/1-e0.001/T0957s2.fasta/T0957s2.aln \
-o e0.001
```
Repeat for various e-value thresholds
### Step3. Evaluate contacts
```bash
/home/badri/PDNET/v3-ICML/scripts/coneva.pl \
-pdb ~/PDNET/v3-ICML/data/casp13/chains/T0957s2.pdb \
-rr ./e0.001/T0957s2.rr
```
