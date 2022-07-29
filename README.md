# DeepMEI
DeepMEI is a convolutional neural network based tool to identify non-reference MEIs in short-read sequencing data
<br/>
![This is an image](https://github.com/xuxif/DeepMEI/blob/main/workflow.png)
<br/>
## Download<br/>
To clone the DeepMEI, use the following command:<br/>
```
git clone https://github.com/xuxif/DeepMEI.git
```
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
2. Create a new conda environment using the deepmei.yml file we provided
 ```
  conda env create -n deepmei -f deepmei.yml
  conda activate deepmei
 ```
3. Pip install pysam and 
4. Conda install samtools, bedtools and BWA
  ```
  conda activate deepmei
  conda install samtools bedtools  bwa -y
  conda install repeatmasker -y
  ```

</br>
