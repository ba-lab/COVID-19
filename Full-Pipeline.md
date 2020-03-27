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
#Input: T0957s2.aln
#Outputs: T0957s2.realdist.distmap.npy and T0957s2.rr
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
### Step4. Visualize distance map and Obtain Rosetta constraints file
```python
import numpy as np

sequence  = 'SNAMINVNSTAKDIEGLESYLANGYVEANSFNDPEDDALECLSNLLVKDSRGGLSFCKKILNSNNIDGVFIKGSALNFLLLSEQWSYAFEYLTSNADNITLAELEKALFYFYCAKNETDPYPVPEGLFKKLMKRYEELKNDPDAKFYHLHETYDDFSKAYPLNN'
predicted_dmap = 'T0957s2.realdist.distmap.npy'
THRESHOLD = 20
CSTWEIGHT = 1.0
MINSEPARATION = 6

cb_map = np.load(predicted_dmap, allow_pickle = True)
cb_map[ cb_map > THRESHOLD ] = THRESHOLD
cb_map[ cb_map < 3.5 ] = 3.5
print(cb_map.shape)
assert len(cb_map[:, 0, 0]) == len(sequence)

#import seaborn as sns
#sns.heatmap(cb_map[:, :, 0], cmap='Spectral')

distances = {}
for j in range(0, len(sequence)):
    for k in range(j, len(sequence)):
        if cb_map[j, k, 0] >= THRESHOLD:
            continue
        if abs(j - k) < MINSEPARATION:
            continue
        distances[str(j+1) + ' ' + str(k+1)] = cb_map[j, k, 0]
print(f"Total constaints: {len(distances)}")

fcst = open('constraints.cst', 'w')
for pair, d in distances.items():
    a, b = pair.split()
    atom1, atom2 = 'CB'
    dev = 0.273 * d - 0.455 # +- 4A deviation
    dmin = d - dev / 2.0
    dmax = d + dev / 2.0
    line = f"AtomPair {atom1} {a} {atom2} {b} SCALARWEIGHTEDFUNC {CSTWEIGHT:.1f} BOUNDED {dmin:.1f} {dmax:.1f} 0.5 NOE"
    fcst.write(line + '\n')
fcst.close()
```
