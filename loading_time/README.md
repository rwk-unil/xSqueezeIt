# Loading time test app

This is a simple application to test loading times from different formats. It will traverse a BCF and extract the genotype data array of each record and finally print the time it took.

## Build

```shell
make
```

## Run

```shell
./loading_time -f chr20.bcf
# Loading time test
# Loading gt data from file chr20.bcf
# Time elapsed = 27[s] 27656[ms] 27656318[us] 
```

```shell
./loading_time -f chr20.bin                
# Loading time test
# Loading gt data from file chr20.bin
# Time elapsed = 15[s] 15833[ms] 15833665[us] 
```