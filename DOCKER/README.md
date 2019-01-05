# Official MACS2 v2.1.2 docker

Maintained and uploaded by Tao Liu <vladimir.liu@gmail.com>

This image is based on MACS2 v2.1.2. MACS2 was compiled using a python2.7 container with gcc, then the executable and library files were tranfered to a python2.7-slim container.

# How to use this MACS2 docker image

## Pull the image

```docker pull fooliu/macs2:macs2v2.1.2```

## Analyze your data

### Example of analyzing ChIP-seq data in the CURRENT directory

First, ```cd``` to the working directory containing ChIP-seq alignment files such as ```chip-seq-file.bam``` and ```control-seq-file.bam```. Then

```docker run -v $PWD:/data/ fooliu/macs2:macs2v2.1.2 macs2 callpeak -t /data/chip-seq-file.bam -c /data/control-seq-file.bam -n test-run --outdir /data/```

The first part ```-v $PWD:/data/``` will mount the CURRENT directory ```$PWD``` to ```/data/``` in the container, so please don't forget to add ```/data/``` to the path of input files with ```-t``` and/or ```-c```, and don't forget to set the ```--outdir``` option to ```/data/```. The final outputs will be directly written into the CURRENT directory. Extra MACS2 options can be modified or added after ```docker run -v $PWD:/data/ fooliu/macs2:macs2v2.1.2 macs```.
