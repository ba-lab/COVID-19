# Full pipeline (example)
Using the CASP13 target T0957s2-D1 as an example.
<hr>  

### Step1. Make alignments
```console
/ssdA/common-tools/DeepMSA/hhsuite2/scripts/build_MSA_e0.001.py ./T0957s2-D1.fasta \
-outdir=./T0957s2-D1/ \
-hhblitsdb=/ssdA/common-tools/uniclust30_2018_08_hhsuite/uniclust30_2018_08 \
-jackhmmerdb=/ssdA/common-tools/uniref_2019_10_24/uniref90.fasta \
-hmmsearchdb=/ssdA/common-tools/fasta_188GB_plass_soedinglab/SRC.fasta\
:/ssdA/common-tools/metaclust-2018-Jun-22/metaclust_nr.fasta
```
Repeat for various e-value thresholds

### Step2. Predict distances and contacts
```console
#Input: T0957s2.aln
#Outputs: T0957s2.realdist.distmap.npy and T0957s2.rr
python3 /home/badri/casp14/35-PrayogRealDistance/src/predict-distance.py \
-w ../35-PrayogRealDistance/job/dist.hdf5.val.hdf5 \
-a ../38-DeepMSAforCASP13FM/1-e0.001/T0957s2.fasta/T0957s2.aln \
-o e0.001
```
Repeat for various e-value thresholds

### Step3. Evaluate contacts (if native is present)
```console
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
    atom1 = atom2 = 'CB'
    if sequence[int(a)-1] == 'G':
        atom1 = 'CA'
    if sequence[int(b)-1] == 'G':
        atom2 = 'CA'
    dev = 0.273 * d - 0.455 # +- 4A deviation
    dmin = d - dev / 2.0
    dmax = d + dev / 2.0
    line = f"AtomPair {atom1} {a} {atom2} {b} SCALARWEIGHTEDFUNC {CSTWEIGHT:.1f} BOUNDED {dmin:.1f} {dmax:.1f} 0.5 NOE"
    fcst.write(line + '\n')
fcst.close()
```

### Step5. Obtain Rosetta fragment files
* Submit fasta file at http://old.robetta.org/fragmentqueue.jsp and collect `aat000_03_05.200_v1_3` and `aat000_09_05.200_v1_3`
* Otherwise, we will need to get the 'make_fragments.pl' (old) or Fragment picker to work

### Step6. Run Rosetta with restraints
```console
/ssdA/common-tools/rosetta_bin_linux_2019.35.60890_bundle/main/source/bin/AbinitioRelax.static.linuxgccrelease \
-database /ssdA/common-tools/rosetta_bin_linux_2019.35.60890_bundle/main/database/ \
-in:file:fasta ../e0.001/T0957s2.fasta \
-in:file:frag3 ./aat000_03_05.200_v1_3 \
-in:file:frag9 ./aat000_09_05.200_v1_3 \
-nstruct 10 -out:pdb \
-constraints:cst_fa_weight 0.5 -constraints:cst_weight 0.5 \
-cst_fa_file ../constraints.cst \
-abinitio:relax -out:overwrite
```

### Step7. Run Rosetta without any restraints (for comparison)
```console
/ssdA/common-tools/rosetta_bin_linux_2019.35.60890_bundle/main/source/bin/AbinitioRelax.static.linuxgccrelease \
-database /ssdA/common-tools/rosetta_bin_linux_2019.35.60890_bundle/main/database/ \
-in:file:fasta ../e0.001/T0957s2.fasta \
-in:file:frag3 ./aat000_03_05.200_v1_3 \
-in:file:frag9 ./aat000_09_05.200_v1_3 \
-nstruct 10 -out:pdb -abinitio:relax -out:overwrite
```

### Step8. Evaluate predicted structures using TM-score and RMSD (if native is present)
* Submit the native structure and the models with minimum score (in the scores.sfc file) to https://zhanglab.ccmb.med.umich.edu/TM-align/.
* In this case:
  * TM-score of best model without CST ~ 0.38
  * TM-score of best model with CST ~ 0.6


