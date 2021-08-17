# VCF / BCF Compressor

VCF / BCF Genotype data compressor based on positional Burrows-Wheeler transform (PBWT) and 16-bit Word Aligned Hybrid (WAH) encoding.

Variants are left in BCF format, genotype data is custom encoded. The genotype data can then be further compressed with standard tools such as gzip.

## Build

This build requires GCC 8+ because modern C++17 features are used.

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
sudo make install
sudo ldconfig
cd ..
# Build application
make
```

## Run

### Compression

```shell
# ./console_app <-c|-x> -f <input file> -o <output file>
mkdir output

# Compression :
./console_app -c -f /path/to/my/data/chr20.bcf -o output/chr20.bin
# This will output two files in output
# output/chr20.bin which is the samples and genotype data in binary encoded format (can still be compressed e.g., with gzip)
# output/chr20.bin_var.bcf which is the variant data, can be opened with bcftools
```

Options :
- `--iota` allows to use natural ordering for checkpoints instead of saving the permutation arrays in the binary file, this results in smaller binary file for collection with many samples. This has no noticeable impact on speed.

### Extraction

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
./console_app -x -f output/chr20.bin | bcftools view -s HG00111,NA12878 | less
./console_app -x -f output/chr20.bin | bcftools stats > output/chr20_stats.txt
```

## File Format Description

### Version 3

#### File organization

This compressor takes an input BCF file and outputs two files :

1) A BCF file with the original variant info, with following fields (`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT`). These field are unaltered (other compressors often replace fields with `.` (missing))
2) A binary file with the "GT" information encoded

Technical details :

the BCF file has a single "sample" `BIN_MATRIX_POS` with a single field `BM` which contains the position of the variant in the "binary matrix" of all variants when encoded as "one-hot" (this is so that multi-allelic variant sites can be encoded efficiently).

The encoded file format has the following organization

| Header                         |
|--------------------------------|
| GT Data Blocks                 |
| Block indices for random access|
| Sample ID's                    |
| Missing positions as sparse    |
| Non standard phase as sparse   |

The header is the main entry point to the file and other data can be accessed from there (offsets to other entries given in header).

##### Blocks

GT Data blocks encode a given number of variants, the offset of each block in the file is given by the block indices, this allows for constant time random access to a block.

The blocks can be compressed with zstd https://github.com/facebook/zstd, a bit in the header indicates if this is the case. If the blocks are compressed two 32-bit words precede each block in memory to indicate the compressed and uncompressed sizes.

Each (uncompressed) block has the following hierarchy :

| Dictionnary |
|-------------|
| Data        |

###### Data
where data are filled with :

| WAH encoded genotypes    |
| Sparse encoded genotypes |
| Binary encoding track    |
| Binary sort track        |

- WAH encoded genotypes are Word Aligned Hybrid encoded genotypes ordered given the sample order until the binary sort track position is "true" where the haplotypes are then reordered given the current variant (Positional Burrows Wheeler Transform).
- Sparse encoded genotypes are the id's of the haplotypes that have the minor allele. (Can be ref or alt, if ref, the number of sparse encodings has MSB bit set, this bit is not used to encode the minor allele, because minor allele means that less than 50% of all haplotypes carry that allele).
- The binary encoding track maps each variant to either a WAH encoding or a sparse encoding.
- The sort track indicates if the following WAH encoded genotypes ordering must be updated by the current variant (PBWT).
Note : the encoding and sort track can be the same track.

###### Dictionnary
The dictionnary is a dictionnary with 32-bit keys and 32-bit values. Terminated by a {-1,-1} key pair ({0xFFFF'FFFF, 0xFFFF'FFFF}).

The dictionnary allows to store the following information :
- Block size, if the block has a smaller size (needed for the last block, because the number of variants is not necesarrily a multiple of block length) this information is put in the dictionnary.
- Offset of the WAH encoded data.
- Offset of the Sparse encoded data.
- Offset of the encoding track.
- Offset of the sort track.
- Other information, if needed by later extension to the format.

##### Sample ID's
A sequence of null terminated strings corresponding to the sample ID's of the input BCF file.

##### Missing positions as sparse
Contains a map of variants positions that have samples with a missing entry "." in BCF/VCF. The samples that have missing entries are encoded as sparse (i.e., list of sample haplotypes that show a missing allele).

##### Non standard phase
The compressor will check out a portion of the input BCF file to see if most samples are phased or unphased. Based on the result the format assumes that most other entries will have samples that are also phased (or unphased).

For each variant where there are samples that are not phased in a mostly phased file (respectively phased in a mostly unphased file). The positions of the samples will be encoded as sparse.

These encoded variant site are then stored in a map as for the missing positions above.

#### Notes

- The missing and non standard phase maps are there to support files that are not entirely phased or unphased as well as files with missing data. Theses files are however not the main focus of this compressor, therefore not much care has been put in the encoding of this. Therefore files with lots of missing data or mixed phase will not benefit from good compression ratios.
- - Improving on this encoding or adding a second compression layer as e.g., zstd could improve this.

### Version 2

**DEPRECATED**

But still supported by compressor / decompressor

#### File organization

This compressor takes an input BCF file and outputs two files :

1) A BCF file with the original variant info, with following fields (`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT`). These field are unaltered (other compressors often replace fields with `.` (missing))
2) A binary file with the "GT" information encoded

Technical details : 

the BCF file has a single "sample" `BIN_MATRIX_POS` with a single field `BM` which contains the position of the variant in the "binary matrix" of all variants when encoded as "one-hot" (this is so that multi-allelic variant sites can be encoded efficiently).

The encoded file format has the following organization

| Header          |
|-----------------|
| Indices WAH     |
| Indices Sparse  |
| WAH data        |
| Sample ID's     |
| Condition track |
| Sparse data     |

The header is the main entry point to the file and other data can be accessed from there

##### Header

The header is 256 bytes long and serves as an entry point for the binary file.

It contains information on the data representation for example the number of bytes used for different representations, e.g., 2 bytes for `uint16_t`.

It contains the offsets (addresses) to the different parts of the file.

The exact C representation can be found below :

```C
struct header_s {
    // "rsvd" fields are "reserved", unused for the moment and kept for future additions

    // 32 bytes
    uint32_t endianness = ENDIANNESS; // Endianness sanity check
    uint32_t first_magic = MAGIC;     // File format sanity check
    uint32_t version = VERSION;       // File format version
    uint8_t  ploidy = PLOIDY_DEFAULT; // Ploidy of samples encoded
    uint8_t  ind_bytes = 0;           // Number of bytes used to save indices
    uint8_t  aet_bytes = 0;           // Number of bytes used to save positions
    uint8_t  wah_bytes = 0;           // Number of bytes used for WAH
    union {
        uint8_t special_bitset = 0;
        struct {
            bool has_missing : 1;     // The input file had missing data
            bool non_uniform_phasing : 1; // The input file has phased / unphased mixed data
            uint8_t rsvd__1 : 6;
        };
    };
    union {
        uint8_t specific_bitset = 0;
        struct {
            bool iota_ppa : 1;        // Reset sort instead of saving permutation arrays
            bool no_sort : 1;         // Data is not permutated
            uint8_t rsvd__2 : 6;
        };
    };
    uint8_t  rsvd_bs[2] = {0,};
    uint32_t rsvd_1[3] = {0,};

    // 64 bytes
    uint64_t hap_samples = 0;         // Number of haplotypes
    uint64_t num_variants = 0;        // Number of variants (total number of ALTs)
    uint32_t block_size = 0;          // DEPRECATED
    uint32_t number_of_blocks = 0;    // DEPRECATED
    uint32_t ss_rate = 0;             // Sub Sample rate of permutation arrays / reset sort rate
    // Offsets, positions of data in the binary file
    uint32_t number_of_ssas = 0;      // Number of sampled loci for random access = ceil(num_variants/ss_rate)
    uint32_t indices_offset = 0;      // Position in the binary file of WAH indices
    uint32_t ssas_offset = 0;         // Position in the binary file of sub sampled permutation arrays (if any)
    uint32_t wahs_offset = 0;         // Position in the binary file of WAH data
    uint32_t samples_offset = 0;      // Position in the binary file of samples (e.g., "NA12878", "HG00101")
    uint32_t indices_sparse_offset = 0; // Position in the binary file of indices for the sparse data
    uint32_t rsvd_offset = 0;
    uint32_t rearrangement_track_offset = 0; // Position in the binary file of the rearrangement track
    uint32_t sparse_offset = 0;       // Position in the binary file of the sparse data

    // 128 bytes
    uint32_t rare_threshold = 0;      // Threshold for the rearrangement track / sorting / wah vs sparse
    uint64_t xcf_entries = 0;         // Num entries in the BCF file (may be less than num_variants if multi-allelic)
    uint8_t rsvd_3[116] = {0,};

    // 32 bytes
    uint32_t rsvd_4[3] = {0,};
    uint32_t sample_name_chksum = 0;; // Checksum unused for now
    uint32_t bcf_file_chksum = 0;     // Checksum unused for now
    uint32_t data_chksum = 0;         // Checksum unused for now
    uint32_t header_chksum = 0;       // Checksum unused for now
    uint32_t last_magic = MAGIC;      // Sanity check magic
} __attribute__((__packed__));
```

##### Indices

Indices allow for random access, they are sampled every `header.ss_rate` variant, they reflect the position of the pointers to the encoded data at the sampled position.

Because sparse data as well as WAH data (similar to RLE) encode a fixed number of bits with a variable number of bits, random access is not directly possible. The indices allow for random access at sampled positions from there access to a given position is possible by traversing the data.

##### Data

Variants are represented internally as bi-allelic 0-1 values (0 is REF) (1 is ALT). The BCF is unaltered in this regard.

Bi-allelic sites with minor allele frequency (MAF) above a given threshold are encoded as WAH, those below are encoded as Sparse (indices of samples with minor allele).

##### Condition track

The condition track represents if a given bi-allelic variant is saved as WAH or as Sparse, as well as, if the variant is used to sort the data.

### Version 1

**DEPRECATED**

- Uses reverse prefix sort, as in the Positional Burrows Wheeler Transform (PBWT), to group related haplotypes together.
- Only uses variants with MAF > 0.01 for sorting
- Can store permutation arrays for random access (costly if many samples) every 8192 variants
- Can reset the sort for faster random access every 8192 variants
  - The sort being reset reduces compression of genotype data but this is offset by the fact that     permutation arrays are not stored
- Encodes the genotype data with Word Aligned Hybrid (WAH) on 16-bits

## Notes / TODO

- Rename the compressor (currently named console_app ...)
- ~~Only outputs data as phased for the moment~~ Done !
- ~~Handle missing~~ Done !
- ~~Only supports bi-allelic sites for the moment~~ Done !
    - (not needed anymore) Convert multi-allelic VCF/BCF to bi-allelic with bcftools :  
      ```shell
      bcftools norm -m any multi_allelic.bcf -o bi_allelic.bcf -O b
      ```

## Further works

- ~~Extraction~~Done !
- ~~Filtering~~ Done !
- ~~Based on the block compression scheme for faster access~~ Done !
- ~~Based on permutation sub sampling for faster / parallel access~~ Done !