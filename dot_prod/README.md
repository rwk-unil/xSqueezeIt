# Dot product test app

This is a simple application to test dot products between phenotypes and genotypes from different formats. It will traverse a BCF and extract the genotype data array of each record and do a dot product with phenotypes to finally print the time it took. When an XSI file is loaded the dot product will use the encoded data structures directly without unpacking them (computation on the compressed data structures in memory).

## Build

```shell
make
```

## Run

```shell
./dot_prod -f chr20.bcf
```

```shell
./dot_prod -f chr20.xsi
```