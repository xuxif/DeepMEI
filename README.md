# DeepMEI
DeepMEI is a convolutional neural network based tool to identify non-reference MEIs in short-read sequencing data
<br/>
![This is an image](https://github.com/xuxif/DeepMEI/blob/main/workflow.png)
<br/>
##
## Software version requirements : <br />
1. samtools 1.15.1<br />
2. bedtools v2.30.0<br />
3. pysam 0.17.0<br />
4. BWA<br />
5. RepeatMasker<br />
6. Tensorflow v2.7.0<br />
7. python 3.8<br />
8. perl v5.32.1<br />
9. Anaconda (optinal)<br />
### Recommended Environment Configuration Steps
1. Insatll Anaconda in your server.
2. Create a new conda environment 
 
 ```
  conda env create -n deepmei 
  conda activate deepmei
 ```
 <br />
 
3. Install Tensorflow 2.7.0

```
  conda install tensorflow=2.7.0
 ```
4. Pip install pysam 

```
  pip install pysam=0.17.0
 ```
5. Conda install samtools, bedtools, BWA and RepeatMasker
 
 ```
  conda activate deepmei
  conda install samtools bedtools  bwa -y
  conda install repeatmasker -y
  ```
#### Configure your server with a .yml file of conda 
The above steps detail the installation process of all dependencies. We also provide a conda environment configuration file (deepmei.yml).Users only need to run the following code to configure the required environment.
 ```
   conda env create -n deepmei -f deepmei.yml
 ```
</br>

##  How to run DeepMEI <br />

### Run DeepMEI
1. Clone the DeepMEI:<br/>

```
  git clone https://github.com/xuxif/DeepMEI.git
```
<br />

2. Move the deepmei folder to the working directory (save all running and final result files of the program) <br />

```
  mv DeepMEI your_workdir
```
<br />

3. input <br/>

-   Bam file aligned with BWA (index is optinal)<br/>
-   Reference genome (include .fai)<br/>


4. Runing DeepMEI <br />

 ```
   conda activate deepmei
   cd your_workdir/DeepMEI
   bash DeepAlu_model/model_test_batch.sh -i test.bam -r reference.fa 
 ```
 <br />
 
#### Parameters:

##### Required

-  -r reference genome (37 or 38 for docker)
-  -i input bam file
-  -v the full path of the DeepMEI folder (except for docker)

##### Optinal
-  -q quick model
-  -b use existed candidate MEIs (skip candidate detection)
-  -m mobile element consensus sequence (fastq format and indexed with BWA)
-  -o output file name
-  -w script directory (for docker)
-  -d sequencing depth (optinal or calculate in chr 1M-2M)
-  -c remove all tempory file (default remove)

### Output

-  1. A vcf file
-  2. Bam file around detected MEI
