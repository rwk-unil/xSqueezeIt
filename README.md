# xSqueezeIt (XSI) - VCF / BCF Genotype Compressor

VCF / BCF Genotype data compressor based on sparse representation for rare variants and positional Burrows-Wheeler transform (PBWT) followed by 16-bit Word Aligned Hybrid (WAH) encoding for common variants. (Minor Allele Frequency threshold is selectable for rare/common variants).

Variant information is left in BCF format to remain compatible with HTSLIB / BCFTools, genotype data is custom encoded as described above. The encoded genotype data can then optionnaly be further compressed with zstd https://github.com/facebook/zstd/.

The compressor was realized with haploid/diploid genotype data in mind (human population genetics data). The compressor supports multi-allelic variant sites. Polyploid samples with ploidy > 2 are not supported yet, mixed ploidy samples are supported (e.g., chrX). The main goal is to provide an alternative file format for storing large human genotype datasets, to reduce loading times and speed-up computation (e.g., with computation on the encoded data directly).

*Note :* This is under development and no stable release / version / tag is available yet. The first stable release will be announced soon.

## Example results

### SHAPEIT4

The `dev` branch on the repository https://github.com/rwk-unil/shapeit4/tree/dev illustrates integration of xSqueezeIt file format support in SHAPEIT4. This serves as an example of how to add xSqueezeIt support in existing tools.

### Compression

Compressing chromosomes 1-22 of 1000 Genomes Phase 3 (1KGP3) https://www.internationalgenome.org/category/phase-3/ (2504 samples, 5008 haplotypes, >88M variants) and chromosomes 1-22 of Haplotype Reference Consortium (HRC) https://www.nature.com/articles/ng.3643 (64976 haplotypes, >39M variants).

<img src="images/compression_rates.png" alt="compression rates" width="600"/>

### Loading

Loading of the data from file format (right CHR1 is vcf.gz) :

<img src="images/loading_times.png" alt="loading times" width="600"/>

HTSlib is used to get the genotype data into an array for each variant entry (record). With xSqueezeIt the genotype data is extracted from the binary compressed file.

Normal loading with HTSlib and traditional BCF :
```C
bcf_sr_add_reader(reader, "chr20.bcf");
...
while (bcf_sr_next_line (reader)) { // While there are records in the BCF
    ...
    bcf_get_genotypes(header, line, &genotype_array, &ngt); // HTSlib
    ...
}
```

Loading genotype data from xSqueezeIt binary file and associated BCF variant info file :
```C
Accessor accessor("chr20.xsi"); // xSqueezeIt binary compressed file
bcf_sr_add_reader(reader, accessor.get_variant_filename());
...
while (bcf_sr_next_line (reader)) { // While there are records in the BCF
    ...
    accessor->get_genotypes(header, line, &genotype_array, &ngt);
    ...
}
```

## Build

### Building the xSqueezeIt command line tool

This build requires GCC 8+ because modern C++17 features are used.

```shell
# Clone
git clone https://github.com/rwk-unil/xSqueezeIt.git
cd xSqueezeIt

# Clone and build htslib (if you already have htslib set Makefile accordingly and skip)
git submodule update --init htslib
cd htslib
autoheader
autoconf
./configure
make
sudo make install
sudo ldconfig
cd ..

# Clone and build zstd (if you already have zstd set Makefile accordingly and skip)
git clone https://github.com/facebook/zstd.git
cd zstd
make
cd ..

# Build application
make
```

### Generating the xSqueezeIt support library

In order to add support for xSqueezeIt into other software it can be exported as sources to build a library.

```shell
make package-sources
```

This will generate a directory named `xsqueezeit_export` which provides a Makefile to build the `libxsqueezeit` library. The library is built with `g++` but provides a C API as well. To integrate into a C application include the file `include/c_api.h` in you sources and link with `libxsqueezeit.a` using `g++`.

An example is given in the `c_api_test` directory, where a simple C program (`main.c` file) is compiled with `cc` (`gcc`) and linked with the library using `g++`.

Integration into C++ software allows to access the C++ internal API, through any of the `.hpp` files. Since this is the internal API it is not necessarily the easiest to use, but C++ programs can use the C API as well, which is well suited if they use HTSLIB because it is very similar. For an example see : https://github.com/rwk-unil/shapeit4/blob/dev/src/io/genotype_reader2.cpp (which is the genotype reader source of SHAPEIT4 with optional support for xSqueezeIt (`#ifdef __XSI__`).

### Dependencies
This software depends on :
- `htslib` https://github.com/samtools/htslib
- `zstd` https://github.com/facebook/zstd/


## Run

### Compression
- `-c,--compress`

```shell
# ./xsqueezeit <-c|-x> -f <input file> -o <output file>
mkdir output

# Compression :
./xsqueezeit -c -f /path/to/my/data/chr20.bcf -o output/chr20.xsi
# This will output two files in output
# output/chr20.xsi which is the samples and genotype data in binary encoded format (can still be compressed e.g., with gzip)
# output/chr20.xsi_var.bcf which is the variant data, can be opened with bcftools
```

Options :
- `--zstd` Compresses blocks with an extra zstd compression layer (only for version 3)
- `--maf <value>` Sets the minor allele frequency (MAF) for the minor allele count (MAC) threshold that selects if a variant is encoded as sparse or word aligned hybrid (WAH), typical values are around 0.001 give or take an order of magnitude

### Extraction
- `-x,--extract`

```shell
# Extraction (requires both files generated above) :
./xsqueezeit -x -f output/chr20.xsi -o output/chr20.bcf # To compressed BCF
./xsqueezeit -x -f output/chr20.xsi > output/chr20.bcf # Alternative command (uncompressed BCF)
```

#### Region extraction
- `-r,--regions <regions>`
- `-R,--regions-file <filename>`

```shell
# Extraction (requires both files generated above) :
./xsqueezeit -x -r "20:200000-200100" -f output/chr20.xsi -o output/chr20.bcf # To compressed BCF
./xsqueezeit -x -r "20:200000-200100" -f output/chr20.xsi | bcftools view # Pipes uncompressed BCF
# The above command is much faster than decompressing and using -r in bcftools
# because only the chosen regions are decompressed, both generate the same result
```

#### Sample extraction
- `-s,--samples <samples>`
- `-S,--samples-file <filename>`

```shell
# Extraction (requires both files generated above) :
./xsqueezeit -x -s HG00101,NA12878 -f output/chr20.xsi -o output/chr20.bcf # To compressed BCF
./xsqueezeit -x -s HG00101,NA12878 -f output/chr20.xsi | bcftools view # Pipes uncompressed BCF
# List of samples in file
./xsqueezeit -x -S samples.txt -f output/chr20.xsi -o output/chr20.bcf
```

### Pipe into bcftools

```shell
# Or pipe directly into bcftools (some examples) :
./xsqueezeit -x -f output/chr20.xsi | bcftools view | less
./xsqueezeit -x -f output/chr20.xsi | bcftools view -s HG00111,NA12878 | less
./xsqueezeit -x -f output/chr20.xsi | bcftools stats > output/chr20_stats.txt
```

The default output of `xsqueezeit` on `stdout` is BCF, to output VCF and be human readable add `-Ov`, in order to make the above operations faster the `-p` or `-Ou` option can be passed to `xsqueezeit` in order to output uncompressed BCF which is the fastest to pipe into BCFTools as mentioned in their [documentation](https://samtools.github.io/bcftools/bcftools.html#common_options) :

> "Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion."

Combining this with sample extraction can even speed up analysis for examples running `bcftools roh` (run of homozigosity) :

With BCFTools :
```shell
time bcftools roh -G30 --AF-dflt 0.4 chr1.bcf -s "NA12878,HG00100,HG00101,HG00102,HG00103" > results.txt
# 3:28.22 total
```

xSqueezeIt with `-p` pipe into BCFTools :
```shell
time xsqueezeit -d -f chr1.xsi -s "NA12878,HG00100,HG00101,HG00102,HG00103" -p | bcftools roh -G30 --AF-dflt 0.4 > results.txt
# 1:16.78 total
```

Running without the `-p` option will result in longer execution times because of the unnecessary xsi->VCF (in xSqueezeIt) and VCF->BCF conversions (in BCFTools). Compared to xsi->(uncompressed)BCF that can directly be processed by BCFTools.

xSqueezeIt without `-p` pipe into BCFTools - (Do not do this, better to use `-p` when piping) :
```shell
time xsqueezeit -d -f chr1.xsi -s "NA12878,HG00100,HG00101,HG00102,HG00103" | bcftools roh -G30 --AF-dflt 0.4 > results.txt
# 8:42.54 total
```
### Explore the "variant-only" generated file

The compressor generates a BCF file without the GT data (so variants only) and a binary file with the compressed GT data. This way the BCF file for the variants can still be explored and used.
```
# ./xsqueezeit <-c|-x> -f <input file> -o <output file>

# Compression :
./xsqueezeit -c -f /path/to/my/data/chr20.bcf -o chr20.xsi
# This will output two files
# chr20.xsi which is the samples and genotype data in binary encoded format (can still be compressed e.g., with gzip)
# chr20.xsi_var.bcf which is the variant data, can be opened with bcftools

bcftools view chr20.xsi_var.bcf | less
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20150218
...
##FORMAT=<ID=BM,Number=1,Type=Integer,Description="Position in GT Binary Matrix">
##XSI=filename.xsi
##bcftools_viewCommand=view tmp.xsi_var.bcf; Date=Tue Aug 17 16:06:34 2021
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BIN_MATRIX_POS
20      60343   rs527639301     G       A       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=20377;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0;EUR_AF=0;SAS_AF=0;AA=.|||;VT=SNP   BM      0
20      60419   rs538242240     A       G       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=19865;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0.001;AA=.|||;VT=SNP    BM      1
20      60479   rs149529999     C       T       100     PASS    AC=17;AF=0.00339457;AN=5008;NS=2504;DP=20218;EAS_AF=0;AMR_AF=0.0043;AFR_AF=0.0106;EUR_AF=0;SAS_AF=0;AA=.|||;VT=SNP      BM      2
```

For example here the AC, AF, AB, etc. entries have not been modified, so they can already be used for some computations, since the BCF file is much smaller than the original one it will faster to process.

This also shows that the compressor keeps all the information that is sometimes lost with other compressors.

In this snippet we can also see the single `BIN_MATRIX_POS` sample with a `BM` entry. This corresponds to the position of the variant(s) in the binary block that is encoded/compressed in the binary file. This allows for random access to the GT data from a given variant in this BCF file. This number is the block:offset (details in file format description). In the snippet above the block is 0 and it increases from 0 to 1 to 2 etc. Multi-allelic variant sites will increase this offset by the number of ALT_ALLELES. For example :

```
20      64139   rs186497980     G       A,T     100     PASS    AC=2,2;AF=0.000399361,0.000399361;AN=5008;NS=2504;DP=20791;EAS_AF=0,0;AMR_AF=0,0;AFR_AF=0,0.0008;EUR_AF=0,0.001;SAS_AF=0.002,0;AA=.|||;VT=SNP;MULTI_ALLELIC     BM      103
20      64150   rs7274499       C       A       100     PASS    AC=102;AF=0.0203674;AN=5008;NS=2504;DP=20555;EAS_AF=0;AMR_AF=0.0058;AFR_AF=0.0741;EUR_AF=0;SAS_AF=0;AA=.|||;VT=SNP      BM      105
```

Where BM increases from 103 to 105 because the previous variant has two ALT_ALLELES A and T. Because in the XSI file each alt allele is stored as a single (encoded) binary line.

The `BM` entry allows to extract GT data directly from a region query on the BCF, this is needed to achieve constant time random access. This also helps when because overlapping a region may not be contiguous (e.g., with indels).

## File Format Description (internal version 4)

The compressor takes an input BCF and output two files :
1) A BCF file with the original variant info, with following fields (`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT`). These fields are unaltered. And a single sample `BIN_MATRIX_POS` with `BM` format field which holds an index at each variant entry.
2) A binary file `.xsi` with the compressed/encoded information (e.g., "GT" data).

<img src="images/xsi_arch.png" alt="XSI architecture" width="1000"/>

### BCF File contents

The BCF (named "filename.xsi_var.bcf") has the corresponding XSI file name in the header and the indices allow to query the data.

```
##FORMAT=<ID=BM,Number=1,Type=Integer,Description="Position in GT Binary Matrix">
##XSI=filename.xsi
```

The `BM` index is a position, a 32-bit integer that can be decomposed as follows : The 15 upper bits represent the block where the data is stored, the 17 lower bits represent an index inside that block. This field is passed to an "accessor" to seek the block that holds the data and then seek inside the block itself.

### XSI File contents

The XSI File has the following organization :

| Header                         |
|--------------------------------|
| Binary Data Blocks             |
| Sample ID's                    |
| Block indices for random access|

- The header is a 256-byte header with version information, compression information, endianness info, and location of the following fields.
- The binary data blocks as optionnaly zstd compressed blocks (the header will tell if they are zstd compressed or not).
- The Sample ID's are a list of sample ID strings. These are the names (IDs) of the samples in the original BCF file.
- Finally a list of indices of the (compressed) blocks, this index is queried to get the location of a binary block inside the file from a `BM` index.

#### XSI Binary blocks

If the binary blocks are zstd https://github.com/facebook/zstd compressed each block consists of :

- A 32-bit number indicating the uncompressed block size.
- Followed by the compressed data.

The 32-bit number allows to allocate the required memory before uncompressing the block.

The uncompressed binary blocks consist of :

| Dictionnary               |
|---------------------------|
| GT Block                  |
| Other data specific block |

The dictionnary is a 32-bit associative table (key-value, both 32-bit). The first key-value pair is -1 (key) and the dictionnary size (without this entry). Reading this first key-value pair allows to know the size of the dictionnary. Each subsequent entry represent a data specific block type (key) and location relative to the start of the current block (value).

For the moment only GT (genotype) data is stored in XSI, the key for this data is 256. 32-bit keys (with -1 being reserved) allow for more than 4 billion data specific types. Ranges of keys can be requested and allocated through a github issue / pull request if developers are interested in encoding other specifc blocks or encode existing data types with alternative compression schemes.

The binary blocks each encode a fixed number of BCF lines (8192 per default). This is an option the can be set during compression. Because for compression the `BM` index is used and this contains the block number and offset inside the block, the decompressor does not need to explicitely know how many BCF lines are in each block.

#### Genotype GT block

The genotype GT blocks (`include/gt_block.hpp`) consist of :

| Dictionnary                 |
|-----------------------------|
| Common variants (PBWT+WAH)  |
| Rare variants (Sparse)      |
| Metadata (missing, phasing) |

As for the binary blocks the dictionnary is a 32-bit (key-value associative table). The first key-value pais is -1 (key) and the dictionnary size (without this entry). The following entries can be found in the dictionnary :

Scalar key-values :

- `0x0 KEY_BCF_LINES` : The number of BCF lines in the block.
- `0x1 KEY_BINARY_LINES` : The number of binary lines in the block (sum of number of alt-alleles for all BCF lines in block for example, a BCF line with multiple alt alleles will be encoded on multiple binary lines (one alt allele per binary line)).
- `0x2 KEY_MAX_LINE_PLOIDY` : The may ploidy encountered in all the lines in the block.
- `0x3 KEY_DEFAULT_PHASING` : The default phase state of the majority of samples, either phased `0|1` or unphased `0/1`.
- `0x4 KEY_WEIRDNESS_STRATEGY` : (For internal use), this desribes the strategy used to compress/encode the missing / end of vector / non default phase information.

Vector key-values : (values are location of the vector)

- `0x10 KEY_LINE_SORT` : Binary vector that indicates if the corresponding binary line is used to sort the PBWT arrangement of subsequent common variant genotypes.
- `0x11 KEY_LINE_SELECT` : Binary vector that indicates if the corresponding binary line is PBWT WAH encoded or sparse encoded.
- `0x12 KEY_LINE_HAPLOID` : Optional binary vector that indicates if the corresponding binary line is haploid. (E.g., for handling chromosome X/Y where some lines can be totally haploid).
- `0x15 KEY_VECTOR_LENGTH` : Unused vector, for future support of mixed ploidy or polyploidy > 2.
- `0x16 KEY_LINE_MISSING` : Optional binary vector that indicates if the corresponding binary line has missing values.
- `0x17 KEY_LINE_NON_UNIFORM_PHASING` : Optional binary vector that indicates if the corresponding binary line has non uniform phasing (mixed phased and unphased samples).
- `0x18 KEY_LINE_END_OF_VECTORS` : Optional binary vector that indicates if the corresponding binary line has "end of vectors", BCF values the represent lower ploidy or lack of information compared to the biggest vector in the BCF line. (This allows support for mixed ploidy).

Matrix key-values : (values are location of the matrix)

- `0x20 KEY_MATRIX_WAH` : Matrix of WAH encoded binary lines
- `0x21 KEY_MATRIX_SPARSE` : Matrix of sparse encoded binary lines
- `0x26 KEY_MATRIX_MISSING` : Matrix of `strategy` encoded binary lines of missing values. Where `strategy` is given by `0x4 KEY_WEIRDNESS_STRATEGY`. Which consists of either (WAH, PBWT+WAH, or Sparse (in matrix below)). For performance the strategy is set to sparse for the moment (not possible to change this from command line). The choice intended for development and testing.
- `0x27 KEY_MATRIX_NON_UNIFORM_PHASING`: Matrix of non uniformly phased positions.
- `0x28 KEY_MATRIX_END_OF_VECTORS` : Same as missing but for "end of vector" BCF entries.
- `0x36 KEY_MATRIX_MISSING_SPARSE` : Sparse matrix for missing values, (both missing matrices are used if the strategy is mixed...)
- `0x38 KEY_MATRIX_END_OF_VECTOR_SPARSE` : Sparse matrix for "end of vector" values, (both end of vector matrices are used if the strategy is mixed...)

Decoding the genotype block uses the vectors (lines) and matrices above to either recreate the original BCF data line or simply give access to the internal data structures.

Internal data structures :

- WAH lines : The WAH (Word Aligned Hybrid) lines are encoded using a 16-bit WAH scheme (see `include/wah.hpp`). If a line updates the PBWT sort all subsequent WAH lines are reordered given the PBWT order.
- Sparse lines : Depending on the number of samples the sparse lines are either encoded on 16-bits (if less than 65536 samples) or 32-bits. Each sparse binary line starts by the number of sparse values followed by the indices of the samples (corresponding samples IDs strings are stored in the XSI file once). The sparse lines are not PBWT reordered because if makes no difference in compression.
- `a` : The PBWT arrangement at the current posision as in the Durbin 2014 PBWT paper. This is not stored in the XSI file but created on the fly on traversal of the file. Binary lines that update this structure are marked by a set bit in the `LINE_SORT` vector.
- Other data structures include the missing, end of vector, and phasing info.

The binary vector lines are themselves WAH encoded to save some space.

#### Sample ID's
A sequence of null terminated strings corresponding to the sample ID's of the input BCF file.

## C API

The C API `include/c_api.h` provides a C compatible API for reading XSI files in a C program. The libxsqueezit library and C API can be exported with the command :

```shell
make package-sources
```

which will create a directory `xsqueezeit_export` with a makefile that allows to build `libxsqueezeit.a`, the xsqueezeit library. Including `include/c_api.h` from the same folder and linking with `libxsqueezeit.a` allow for mixed XSI/VCF/BCF file reading from C programs.

An example application can be found in `c_api_test`.

Another example is the addition of XSI support in SHAPEIT4 https://github.com/rwk-unil/shapeit4/tree/dev. Specifically see : https://github.com/rwk-unil/shapeit4/blob/dev/src/io/genotype_reader2.cpp (`#ifdef __XSI__` for reference).

## Allele count and Allele number computation

See directory `af_stats` which provides a program that recomputes allele counts and allele numbers from an XSI file.

## Dot products

See directory `dot_prod` which provides a benchmark that computes dot products between all bi-allelic sites and a random phenotype vector. The genotypes can be loaded from either a BCF or XSI files. Internal data structures and "compressive acceleration" is used when an XSI file is provided.

## Loading time

Loads all genotypes of a file (either XSI or BCF) into memory, one line at the time, once every line has been loaded once the time is shown. This allows to benchmark data loading from either format. The HTSLIB is used in both cases, only the method `bcf_get_genotypes()` is replaced by our own decompression when an XSI file is used.

## Lockstep loader

Load genotypes from two files and checks that the two files have the same genotypes for all variants they contain. This allows to check if an XSI file and a BCF file hold the same data. This also allows to compare XSI against XSI or BCF against BCF.

Example :

```
lockstep_loader % ./lockstep_loader --file1 ../../Data/1kgp3/chrX.bcf --file2 ../chrX/chrX.xsi
Lockstep loader test
ext1 bcf ext2 xsi
Checked 17368209744 GT entries
Files have the same GT data
```

## Validation through integration testing

The directory `test` provides cukinia scripts to run integration tests that check that the features of XSI work.

- Check if missing data works
- Check if end of vector works (samples with smaller ploidy mixed in)
- Check if non uniform phasing works (samples that are unphased at some variants)
- Check if missing and non uniform phasing work together
- Check if some variants can have a global different ploidy than other (does not yet work)
- Check if files that span multiple blocks are recovered
- Check if zstd compression works
- Check if region extraction works
- Check if sample extraction works
- Check combinations of the above...

The tests rely on comparing the output of `bcftools view` of the XSI compressed-decompressed files through the utility `diff` against the reference input files, diff may run out of memory on big files (because the uncompressed textual result of bcftools is checked, i.e. plain VCF). These test allow to check that the files contain the same data.

For testing bigger files we recommend either using the lockstep loader above (for checking GT data) or bcftools to compare the rest of the data.