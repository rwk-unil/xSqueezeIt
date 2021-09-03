```shell
# Compression :
./xsqueezeit -c -f /path/to/my/data/chr20.bcf -o output/chr20.bin
# Extraction :
./xsqueezeit -x -f output/chr20.bin -o output/chr20.bcf # To BCF
```

```shell
# Or pipe directly into bcftools (some examples) :
./xsqueezeit -x -f output/chr20.bin | bcftools view | less
./xsqueezeit -x -f output/chr20.bin | bcftools view -s HG00111,NA12878 | less
./xsqueezeit -x -f output/chr20.bin | bcftools stats > output/chr20_stats.txt
```