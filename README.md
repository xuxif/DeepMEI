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

##  How to run DeepMEI
### Download<br/>
1. Clone the DeepMEI:<br/>
```
git clone https://github.com/xuxif/DeepMEI.git
```
<br />
2. Move the deepmei folder to the working directory (save all running and final result files of the program) <br />
```
mv DeepMEI <your_workdir>
```
3. input
 3.1 Indexed bam (index is optinal)
 3.2 Reference genome (include .fai)
4. Runing DeepMEI <br />
 ```
 cd your_workdir/DeepMEI
 bash DeepAlu/model_test_batch.sh -i <bam_file_path.bam> -r <genome_file_path.fa> 
 ```
 <br />
