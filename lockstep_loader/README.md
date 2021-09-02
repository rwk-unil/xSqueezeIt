# Lockstep loader app

This application takes two files as input and will check if they both contain the same genotype data.

## Build

```shell
make
```

```shell
./lockstep_loader --file1 chr20.bcf --file2 chr20.bin
# Lockstep loader test
# Files have the same GT data
```

```shell
./lockstep_loader --file1 chr20.bin --file2 chr20_zstd.bin
# Lockstep loader test
# Files have the same GT data
```

```shell
./lockstep_loader --file1 chr20.bin --file2 chr20_small.bcf
# Lockstep loader test
# Files don't have the same variants !
```