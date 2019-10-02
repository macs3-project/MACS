# Official MACS2 v2.1.4 docker

MACS2 is a bioinformatics algorithm to analyze ChIP-seq datasets. 

# Get Started

## Pull the image

Currently, there are two types of MACS2 images, one based on official Python3.7 docker image, and the other based on official Python3.7-slim docker image. Only if space is a concern, and you keep using Python3.7-slim docker for all the Python 3.7 applications, pull the slim one.

To pull the regular MACS2 image:

```docker pull fooliu/macs2```

To pull the MACS2 image based on Python3.7-slim

```docker pull fooliu/macs2:py27-slim```

## Analyze your data

### Example of analyzing ChIP-seq data in the CURRENT directory

Let's assume you pulled the regular MACS2 image. If not, change the image name accordingly. First, ```cd``` to the working directory containing ChIP-seq alignment files such as ```chip-seq-file.bam``` and ```control-seq-file.bam```. Then

```docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/chip-seq-file.bam -c /data/control-seq-file.bam -n test-run --outdir /data/```

The first part ```-v $PWD:/data/``` will mount the CURRENT directory ```$PWD``` to ```/data/``` in the container, so please don't forget to add ```/data/``` to the path of input files with ```-t``` and/or ```-c```, and don't forget to set the ```--outdir``` option to ```/data/```. The final outputs will be directly written into the CURRENT directory. Extra MACS2 options can be modified or added after ```docker run -v $PWD:/data/ fooliu/macs2```. The ENTRYPOINT (or the default command when run the container) has been set as ```macs2```.

# Built with

* Python3.7 docker image 3.7
* pip install numpy version 1.17, cython 0.29, and pytest 4.6 
* git clone MACS codes from github master branch then run ```setup.py```

# Author

This Docker image is maintained and uploaded by Tao Liu <vladimir.liu@gmail.com>. MACS2 is actively maintained on Github, many users contribute codes and suggestions. 
