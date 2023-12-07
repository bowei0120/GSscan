# GSscan: An LSTM-based deep-learning framework capturing genomic instability-related status using low-pass whole genome sequencing (WGS) data

## Introduction
GSscan, a long short-term memory (LSTM)-based deep-learning framework, which utilizes low-pass whole genome sequencing (WGS) data to capture genomic instability-related features. Here, we provide a GSscan model directly surveys homologous recombination deficiency (HRD) status independent of other existing biomarkers.

## Usage
#### 1. Generating bin count from BAM file
We use nextflow pipline to generate the bin counts file
+ step 1 : download the docker image
```
docker pull lengyang/yj_bin
```
+ step 2 : install nextflow
  
```
curl -s https://get.nextflow.io | bash
```

+ step3: modify the nextflow.config

```
params {
    outbin_alias="1m" ## output bin name
    inbin=100 # minimal scale of bin, 100 means 100k
    outbin=1000 # output scale of bin, 1000 means 1000k (1m)
    parallel_number=10 # number of parallel sample
}

process{ 
  container="lengyang/yj_bin:latest"
  containerOptions="-v /mnt:/mnt" # Map local path to docker image
}
```

+ step 4 : generate bin counts from bam 
```
nextflow run ./workflow.nf \
    --input samplesheet.csv \ # path to sample sheet
    --out_dir  ./result \ # path to output dir
    -with-docker \
    -bg \
    -resume
```
+ input :samplesheet.csv
example of samplesheet.csv: column 1 is sample id and columns 2 is path to bam file
```
sample_id,bam
test1,/mnt/results/test1.bam
test2,/mnt/results/test2.bam
```
+ output: bin counts file （./result/yj_bin_merge_1m.txt）

#### 2. Predict HRD scores using GSscan
+ Usage
```
usage: predict.py [-h] -i IN_PATH -w WEIGHT_PATH [-o OUTPUT_FILE]

Call GSscan model to predict HRD risk score from bin count file

optional arguments:
  -h, --help            show this help message and exit
  -i IN_PATH, --in-path IN_PATH
                        Path of bin count file
  -w WEIGHT_PATH, --weight-path WEIGHT_PATH
                        directory where torch model and model weight is stored
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        output file
```


## Cite

