# DeepMEI
DeepMEI is a convolutional neural network based tool to identify non-reference MEIs in short-read sequencing data
<br/>
![This is an image](https://github.com/xuxif/DeepMEI/blob/main/workflow.png)
<br/>
##
We are excited to announce significant updates to DeepMEI, making it more efficient and adaptable for genomic analysis across various computational settings. Our continuous efforts to optimize DeepMEI have led to notable improvements in runtime efficiency, especially in environments with higher computational capacity.

Docker Version Update to v1.6.26:
The latest Docker version of DeepMEI, v1.6.26, has been enhanced to provide even better performance across a range of hardware configurations. Our testing, using the NA12878 37.3X whole-genome sequencing data, has shown remarkable efficiency improvements:

On a system with 12 CPUs and 16 GB of memory, DeepMEI now completes its run in 83 minutes.
With 32 CPUs and 64 GB memory, the runtime significantly decreases to just under 40 minutes.
Maximizing the available resources with 64 CPUs and 256 GB memory, DeepMEI completes the task in approximately 30 minutes.
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
2. Replace the placeholder 'your_bam_file.bam' with the actual file name 'input.bam'. DeepMEI  requires input BAM files to be aligned using BWA mem. This requirement ensures optimal performance and accuracy in MEI detection.


  GRCh38:
```
  sudo docker run -it  -v /Bind_mount_a_volume_to_include_input_bam_file/:/root/data/ -w /root xuxiaofeiscu/deepmei:latest  /bin/bash -c 'export PATH=/root/miniconda3/bin:$PATH;./DeepMEI/DeepMEI -i /root/data/you_bam_file.bam  -r 38 -w /root/data/'
```
  hs37d5, hg19 or GRCh37:
```
  sudo docker run -it  -v /Bind_mount_a_volume_to_include_input_bam_file/:/root/data/ -w /root xuxiaofeiscu/deepmei:latest  /bin/bash -c 'export PATH=/root/miniconda3/bin:$PATH;./DeepMEI/DeepMEI -i /root/data/you_bam_file.bam  -r 19 -w /root/data/'
```
  chm13
```
  sudo docker pull xuxiaofeiscu/deepmei:chm13
  sudo docker run -it  -v /Bind_mount_a_volume_to_include_input_bam_file/:/root/data/ -w /root xuxiaofeiscu/deepmei:chm13  /bin/bash -c 'export PATH=/root/miniconda3/bin:$PATH;./DeepMEI/DeepMEI -i /root/data/you_bam_file.bam  -r /root/DeepMEI/DeepMEI_model/reference/chm13v2.0.fa -w /root/data/'
```
## Installation and Usage with Apptainer (Singularity) <br />
We have provided an Apptainer (formerly Singularity) version of DeepMEI for easy deployment and use. Follow these steps to get started:
1. Download the DeepMEI Apptainer Image
You can download the DeepMEI Apptainer image from Zenodo:
Download <a href="https://zenodo.org/records/14928972/files/deepmei.sif?download=1">deepmei.sif</a>
2. Convert the Image to a Sandbox (Optional but Recommended)
Converting the .sif image to a sandbox allows for easier modification and access to the container's filesystem.
```
  apptainer build --sandbox deepmei_app deepmei.sif
```
Or
```
  singularity build --sandbox deepmei_app deepmei.sif
```
This will create a directory named deepmei_app containing the sandboxed environment.

3. Run DeepMEI
Once converted the image into a container , you can run DeepMEI with the following command:
```
apptainer exec --cleanenv --no-home --writable --bind /path/to/input/data:/opt/data deepmei_app/ DeepMEI -i /opt/data/you_bam_file.bam -r /opt/data/DeepMEI/DeepMEI_model/reference/hs37d5.fa
```
Or
```
singularity exec --cleanenv --no-home --writable --bind /path/to/input/data:/opt/data deepmei_app/ DeepMEI -i /opt/data/you_bam_file.bam -r /opt/data/DeepMEI/DeepMEI_model/reference/hs37d5.fa
```
Replace /path/to/input/data with the actual path to your input bam file and correct referece file. <br />
The container includes pre-loaded reference genomes:<br />
‌hs37d5 (hg19)‌<br />
  Path: /opt/data/DeepMEI/DeepMEI_model/reference/hs37d5.fa<br />
‌GRCh38 (Homo_sapiens_assembly38)‌<br />
  Path: /opt/data/DeepMEI/DeepMEI_model/reference/Homo_sapiens_assembly38.fasta<br />
The output will be generated in a DeepMEI_output folder within the location of the BAM file on the host system.<br />
## Software version requirements (without container): <br />
1. samtools 1.15.1 (Other versions need to test whether the "samtools coverage and samtools import" function is included)<br />
2. bedtools v2.30.0<br />
3. pysam 0.17.0<br />
4. BWA<br />
5. RepeatMasker 4.1.2<br />
6. Tensorflow v2.7.0 (before v2.12)<br />
7. python 3.8 or 3.9 (Versions above python3.10 may encounter compatibility issues, such as the inability to install tensorflow 2.7.0 correctly）<br />
8. perl v5.32.1<br />
9. bc <br />
11. Anaconda <br />
12. xarg v4.5 (v4.8 will prompt parameter conflict but still work properly）

### Configure your server with Conda (Recommanded)
Due to the large size of the DeepMEI model file (over 600 MB), it is not feasible to include it in the conda package. Users will need to download the DeepMEI code separately from GitHub. To set up the required environment, users should run the following command:
 ```
  conda create -n deepmei -c conda-forge -c bioconda 'tensorflow>=2.7.0,<2.12' 'samtools>=1.15' 'pysam=0.17' 'bedtools=2.30' 'bwa>=0.7.17' 'repeatmasker=4.1.2' 'python>=3.8,<3.10' 'perl=5.32' bc -y
  conda activate deepmei
 ```
  Installing DeepMEI with Conda might take a long time, primarily because the installation of RepeatMasker requires significant time to download the library files.
</br>

### Configure your server with a .yml file of Conda (Recommanded)
The above steps detail the installation process of all dependencies. We also provide a conda environment configuration file (deepmei.yml).Users only need to run the following code to configure the required environment.
 ```
   conda env create -f deepmei.yml
   conda activate deepmei
 ```
</br>


### Manual configure your server with conda and pip
1. Insatll Anaconda in your server (Not Miniconda which met errors).
2. Create a new conda environment 
 
 ```
  conda env create -n deepmei 
  conda activate deepmei
 ```
 <br />
 
3. Install Tensorflow 2.7.0 or later ( Do not support Tensorflow 2.12 or later)

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

##  Install DeepMEI code from GitHub<br />

### Run DeepMEI
1. Clone the DeepMEI:<br/>

```
  git clone https://github.com/xuxif/DeepMEI.git
```
<br />

2. Move the deepmei folder to the working directory. (Optional） <br />

```
  mv DeepMEI your_workdir
```

  Add add the location of DeepMEI to your PATH.
```
  export PATH=/Path/to/DeepMEI/:$PATH
```
To add the path permanently, append PATH=/Path/to/DeepMEI/:$PATH to the ~/.bash_profile file.

<br />

3. Input <br/>

-   Bam file aligned with BWA (index is optinal)<br/>
-   Reference genome (include .fai and .dict)<br/>


4. Runing DeepMEI <br />

 ```
   conda activate deepmei
   cd /path/to/DeepMEI
  /path/to/DeepMEI/DeepMEI -i /path/to/your/bam_file/ -r /path/to/reference_sequence.fa
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
-  -p threads available
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
