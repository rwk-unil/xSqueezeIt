# VCF / BCF Compressor

## Run

### Compression

```shell
# ./console_app <-c|-x> -f <input file> -o <output file>
```

```shell
# Compression :
./squishit -c -f /path/to/my/data/chr20.bcf -o output/chr20.bin
```

Options :
- `--V3` Compresses with version 3
- `--V2` Compresses with version 2 (deprecated)
- `--zstd` Compresses blocks with an extra zstd compression layer (only for version 3)
- `--maf <value>` Sets the minor allele frequency (MAF) for the minor allele count (MAC) threshold that selects if a variant is encoded as sparse or word aligned hybrid (WAH), typical values are around 0.001 give or take an order of magnitude
- `reset-sort-block-length <value>` Sets the size of the encoded blocks in number of variants. A bigger size can results in better compression, a smaller size can result in faster random access.

### Extraction

```shell
# Extraction :
./console_app -x -f output/chr20.bin -o output/chr20.bcf # To BCF
```

```shell
# Extraction (requires both files generated above) :
./console_app -x -f output/chr20.bin -o output/chr20.bcf # To compressed BCF
./console_app -x -f output/chr20.bin > output/chr20.bcf # Alternative command (uncompressed BCF)
```

#### Region extraction
```shell
# Extraction (requires both files generated above) :
./console_app -x -r "20:200000-200100" -f output/chr20.bin -o output/chr20.bcf # To compressed BCF
./console_app -x -r "20:200000-200100" -f output/chr20.bin | bcftools view # Pipes uncompressed BCF
# The above command is much faster than decompressing and using -r in bcftools
# because only the chosen regions are decompressed, both generate the same result
```

#### Sample extraction
```shell
# Extraction (requires both files generated above) :
./console_app -x -s HG00101,NA12878 -f output/chr20.bin -o output/chr20.bcf # To compressed BCF
./console_app -x -s HG00101,NA12878 -f output/chr20.bin | bcftools view # Pipes uncompressed BCF
```

### Pipe into bcftools

```shell
# Or pipe directly into bcftools (some examples) :
./console_app -x -f output/chr20.bin | bcftools view | less
./console_app -x -f output/chr20.bin | bcftools stats > output/chr20_stats.txt
```