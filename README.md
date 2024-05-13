# DeepMEI
DeepMEI is a convolutional neural network based tool to identify non-reference MEIs in short-read sequencing data
<br/>
![This is an image](https://github.com/xuxif/DeepMEI/blob/main/workflow.png)
<br/>
##
We are excited to announce significant updates to DeepMEI, making it more efficient and adaptable for genomic analysis across various computational settings. Our continuous efforts to optimize DeepMEI have led to notable improvements in runtime efficiency, especially in environments with higher computational capacity.

Docker Version Update to v1.6.25:
The latest Docker version of DeepMEI, v1.6.25, has been enhanced to provide even better performance across a range of hardware configurations. Our testing, using the NA12878 37.3X whole-genome sequencing data, has shown remarkable efficiency improvements:

On a system with 12 CPUs and 16 GB of memory, DeepMEI now completes its run in 96 minutes.
With 32 CPUs and 64 GB memory, the runtime significantly decreases to just under 50 minutes.
Further increasing resources to 32 CPUs and 128 GB memory, DeepMEI finishes in 48 minutes.
Maximizing the available resources with 64 CPUs and 256 GB memory, DeepMEI completes the task in approximately 40 minutes.
These timings highlight DeepMEI's enhanced efficiency and scalability, confirming its suitability for high-throughput settings where runtime is critical.

Recommendation:
For ease of use and to avoid server configuration and dependency management challenges, we highly recommend using the Docker version of DeepMEI. Our team frequently updates this version to ensure optimal performance across a wide range of computational environments. For instance, DeepMEI can analyze 30x WGS sequence data in significantly reduced times, demonstrating its adaptability and performance on both personal and server-grade computers.

Should you require faster analysis, please do not hesitate to reach out to us. Our goal is to continuously improve DeepMEI, making genomic analysis more efficient and accessible to researchers worldwide.
## Docker: <br />
Pull docker image from docker hub
```
  docker pull xuxiaofeiscu/deepmei:latest
```
### To reference your bam/cram file (along with its bai/crai index) located in /home/ubuntu/bam/input.bam, please follow these steps:<br /><br />
1. Replace the placeholder 'Bind_mount_a_volume_to_include_input_bam_file' with the directory path '/home/ubuntu/bam/' in your command.
2. Replace the placeholder 'your_bam_file.bam' with the actual file name 'input.bam'.
  GRCh38:
```
  sudo docker run -it  -v /Bind_mount_a_volume_to_include_input_bam_file/:/root/data/ -w /root xuxiaofeiscu/deepmei:latest  /bin/bash -c 'export PATH=/root/miniconda3/bin:$PATH;bash DeepMEI/DeepMEI_model/model_test_batch.sh -i /root/data/you_bam_file.bam  -r 38 -w /root/data/'
```
  hs37d5, hg19 or GRCh37:
```
  sudo docker run -it  -v /Bind_mount_a_volume_to_include_input_bam_file/:/root/data/ -w /root xuxiaofeiscu/deepmei:latest  /bin/bash -c 'export PATH=/root/miniconda3/bin:$PATH;bash DeepMEI/DeepMEI_model/model_test_batch.sh -i /root/data/you_bam_file.bam  -r 19 -w /root/data/'
```
## Software version requirements (without docker): <br />
1. samtools 1.15.1 (Other versions need to test whether the "samtools coverage" function is included)<br />
2. bedtools v2.30.0<br />
3. pysam 0.17.0<br />
4. BWA<br />
5. RepeatMasker<br />
6. Tensorflow v2.7.0<br />
7. python 3.8 or 3.9 (Versions above python3.10 may encounter compatibility issues, such as the inability to install tensorflow 2.7.0 correctly）<br />
8. perl v5.32.1<br />
9. bc <br />
11. Anaconda <br />
12. xarg v4.5 (v4.8 will prompt parameter conflict）
### Recommended Environment Configuration Steps (with conda and pip)
1. Insatll Anaconda in your server (Not Miniconda which met errors).
2. Create a new conda environment 
 
 ```
  conda env create -n deepmei 
  conda activate deepmei
 ```
 <br />
 
3. Install Tensorflow 2.7.0

```
   conda install -c conda-forge tensorflow=2.7.0
 ```
   Or
 ```
   pip install tensorflow==2.7.0
   pip install protobuf==3.20.* 
 ```
4. Pip install pysam 

```
  pip install pysam==0.17.0
 ```
5. Conda install samtools, bedtools, BWA and RepeatMasker
 
 ```
  conda install samtools bedtools  bwa -y
  conda install repeatmasker -y
  ```

#### Configure your server with a .yml file of conda (Optional)
The above steps detail the installation process of all dependencies. We also provide a conda environment configuration file (deepmei.yml).Users only need to run the following code to configure the required environment.
 ```
   conda env create -f deepmei.yml
 ```
</br>

#### Configure your server with bioconda (Optional)
Due to the large size of the DeepMEI model file (over 600 MB), it is not feasible to include it in the conda package. Consequently, the bioconda package for DeepMEI only contains the necessary environment configuration requirements. Users will need to download the DeepMEI code separately from GitHub. To set up the required environment, users should run the following command:
 ```
  conda env create -n deepmei 
  conda activate deepmei
  conda install -c bioconda deepmei -y
 ```
</br>
##  Install DeepMEI code from GitHub<br />

### Run DeepMEI
1. Clone the DeepMEI:<br/>

```
  git clone https://github.com/xuxif/DeepMEI.git
```
<br />

2. Move the deepmei folder to the working directory and do not rename the directory 'DeepMEI'. (Optional） <br />

```
  mv DeepMEI your_workdir
```
<br />

3. Input <br/>

-   Bam file aligned with BWA (index is optinal)<br/>
-   Reference genome (include .fai and .dict)<br/>


4. Runing DeepMEI <br />

 ```
   conda activate deepmei
   cd /path/to/DeepMEI
  /path/to/DeepMEI/DeepMEI -i /path/to/your/bam_file/ -r /path/to/reference_sequence.fa -w  /path/to/DeepMEI output directory/
 ```
 <br />
 
#### Parameters:

##### Required

-  -r reference genome (full path required), reference genome should have reference.fai and reference.dict (Could be generated by 'samtools faidx reference.fa, samtools dict reference.fa >reference.dict')
-  -i input bam file (full path required), 1. bam file shuld have .bai file ;2. If input is cram file a corresponded reference.fa should offered with -r   

##### Optinal
-  -q quick model
-  -b use existed candidate MEIs (skip candidate detection)
-  -m mobile element consensus sequence (fastq format and indexed with BWA)
-  -o output file name
-  -w script directory 
-  -d sequencing depth (optinal or calculate in chr 1M-2M)
-  -c remove all tempory file (default remove)
-  -j Joint calling with existing MEI locus
-  -w the full path of the DeepMEI output folder (except for docker)

### Output

-  1. A VCF file
-  2. Bam file of detected MEIs

### Joint calling for MEIs
   If your ananlysis contains a cohort (or pedigree) a joint calling strategy is necessary. We have test joint calling in 1kGP high coverage WGS and a greatly improvement of sensitivity have been made. In the next version we will provide a joint calling module for DeepMEI. Plesase let us know if you are intersted in this module and we will send you a non-final version, but with joint-calling.
