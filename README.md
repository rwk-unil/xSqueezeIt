# VCF / BCF Compressor

## Build

```shell
# Clone
git clone https://github.com/rwk-unil/pbwt_exp.git #TODO Public repo
cd pbwt_exp
# Clone and build htslib (if you already have htslib set Makefile accordingly and skip)
git submodule update --init htslib
cd htslib
autoheader
autoconf
./configure
make
cd ..
# Build application
make
```

## Run

### Compression

```shell
# ./console_app <-c|-x> -f <input file> -o <output file>
mkdir output

# Compression :
./console_app -c -f /path/to/my/data/chr20.bcf -o output/chr20.bin
# This will output two files in output
# output/chr20.bin which is the samples and genotype data in binary encoded format (can still be compressed e.g., with gzip)
# output/chr20.bin_var.bcf which is the variant data, can be opened with bcftools
```

### Extraction

```shell
# Extraction :
./console_app -x -f output/chr20.bin -o output/chr20.bcf
./console_app -x -f output/chr20.bin > output/chr20.bcf # Alternative command
```

### Pipe into bcftools

```shell
# Or pipe directly into bcftools (some examples) :
./console_app -x -f output/chr20.bin | bcftools view | less
./console_app -x -f output/chr20.bin | bcftools view -s HG00111,NA12878 | less
./console_app -x -f output/chr20.bin | bcftools stats
```

## Notes

- Only supports phased data for the moment
- Only supports bi-allelic sites for the moment
    - Convert multi-allelic VCF/BCF to bi-allelic with bcftools :  
      ```shell
      bcftools norm -m any multi_allelic.bcf -o bi_allelic.bcf -O b
      ```