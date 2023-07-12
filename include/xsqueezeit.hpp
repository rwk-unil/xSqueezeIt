/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef __XSQUEEZEIT_HPP__
#define __XSQUEEZEIT_HPP__

#include "CLI11.hpp"

#define XSI_BCF_VAR_EXTENSION "_var.bcf"

class GlobalAppOptions
{
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name, default is stdio");
        app.add_option("-o,--output", ofname, "Output file name, default is stdio");
        app.add_option("-O, --output-type", output_type, "Output type b (BCF) u (uncompressed BCF) z (VCF.gz) v (VCF) x (XSI)");
        app.add_flag("-p,--fast-pipe", fast_pipe, "Outputs uncompressed BCF (-Ou) when writing to stdout (default)");
        app.add_flag("-c,--compress", compress, "Compress");
        app.add_flag("-v,--verbose", verbose, "Verbose, prints progress");
        app.add_flag("-d,--decompress", decompress, "Decompress");
        app.add_flag("-x,--extract", decompress, "Extract (Decompress)");
        //app.add_flag("--wait", wait, "DEBUG - wait for int input");
        //app.add_flag("--verify", verify, "DEBUG - verify");

        //app.add_flag("--iota", iota, "DEBUG - Seed PWBT with natural order");
        //app.add_flag("--no-sort", no_sort, "DEBUG - No PBWT sort");
        //app.add_flag("--v2", v2, "DEPRECATED - Version 2");
        //app.add_flag("--v3", v3, "DEBUG - New version");
        //app.add_flag("--v33", v33, "DEBUG - Experimental Streaming Version 3");
        //app.add_flag("--exp", v4, "DEBUG - Experimental Version");
        app.add_flag("--zstd", zstd, "Compress blocks with zstd");
        app.add_option("--zstd-level,--zl", zstd_compression_level, "zstd compression level");
        app.add_option("--maf", maf, "Minor Allele Frequency threshold");
        app.add_flag("-i,--info", info, "Get info on file");
        //app.add_option("--variant-block-length", reset_sort_block_length, "V3 - variant block length");

        //app.add_flag("--sandbox", sandbox, "DEBUG - ...");
        //app.add_flag("--inject-phase-switches", inject_phase_switches, "DEBUG injects phase switches");
        //app.add_option("--phase-switch-prob", phase_switch_prob, "DEBUG probability for phase switch injection");
        //app.add_flag("--compare-matrices", compare_matrices, "DEBUG compare matrices from two bi-allelic bcf files");
        //app.add_flag("--copy-bcf", copy_bcf, "DEBUG copies a bcf (not useful)");
        //app.add_flag("--phase", phase, "DEBUG - phases the file according to PBWT sort");
        //app.add_flag("--phase2", phase_2, "DEBUG - phases the file given samples");
        //app.add_flag("--het-info", het_info, "DEBUG - displays number of het sites per sample");
        //app.add_flag("--compute-phase-switch-errors", compute_phase_switch_errors, "Compute phase switch errors between input and output file");
        //app.add_flag("--count-xcf", count_xcf, "DEBUG - counts number of variant entries in VCF/BCF");
        //app.add_flag("--create-map", create_map, "DEBUG - create map");
        //app.add_flag("--unphase", unphase, "Removes phasing and reorders alleles in natural order e.g., 1|0 => 0/1");
        //app.add_flag("--unphase-random", unphase_random, "Removes phasing and reorders alleles randomly");
        //app.add_flag("--sprinkle_missing", sprinkle_missing, "Sprinkles missing values randomly with 1% chance");
        //app.add_flag("--bitmap", bitmap, "DEBUG - creates bitmap");
        //app.add_flag("--bitmap_pbwt", bitmap_pbwt, "BEBUG - apply PBWT in bitmap");
        //app.add_flag("--bitmap_het", het_bitmap, "DEBUG - HET Bitmap, possible PBWT option");
        //app.add_flag("--color_bitmap16", color_bitmap16, "DEBUG - creates bitmap");
        //app.add_flag("--sorted_bitmap", sorted_bitmap, "DEBUG - creates sorted bitmap");
        //app.add_flag("--block_sorted_bitmap", block_sorted_bitmap, "DEBUG - creates block sorted bitmap");
        //app.add_option("--bblock", block_size, "DEBUG - block size for block sorted bitmap");
        //app.add_flag("--histogram_info", histogram_info, "DEBUG - some histogram info");
        //app.add_flag("--partial_pbwt", partial_pbwt, "DEBUG - creates partial tree like pbwt bitmap");
        /*
         * Options and descriptions below are copied from : https://samtools.github.io/bcftools/bcftools.html#common_options
         * They were copied to remain coherent with the bcftools interface (uses htslib underneath)
         */
        app.add_option("-s,--samples", samples, "Comma-separated list of samples to include or exclude if prefixed with \"^\".");
        app.add_option("-S,--samples-file", samples_file, "File of sample names to include or exclude if prefixed with \"^\". One sample per line.");
        app.add_option("-r,--regions", regions, "chr|chr:pos|chr:beg-end|chr:beg-[,...]"); //"Comma-separated list of regions, see also -R, --regions-file. Overlapping records are matched even when the starting coordinate is outside of the region, unlike the -t/-T options where only the POS coordinate is checked. Note that -r cannot be used in combination with -R.");
        app.add_option("-R,--regions-file", regions_file, "Region file (same as bcftools)"); //"Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file can contain either positions (two-column format) or intervals (three-column format): CHROM, POS, and, optionally, END, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and END (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, "chr20" is not the same as "20". Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.");
        app.add_option("-t,--targets", targets, "[^]chr|chr:pos|chr:from-to|chr:from-[,...]"); //"Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with "^" to request logical complement. For example, "^X,Y,MT" indicates that sequences X, Y and MT should be skipped. Yet another difference between the -t/-T and -r/-R is that -r/-R checks for proper overlaps and considers both POS and the end position of an indel, while -t/-T considers the POS coordinate only. Note that -t cannot be used in combination with -T.");
        //app.add_option("-T,--targets-file", targets_file, ""); // "Same -t, --targets, but reads regions from a file. Note that -T cannot be used in combination with -t.\nWith the call -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele. Note that the file must be compressed and indexed.");
        app.add_flag("-H,--no-header", no_header, "Suppress the header in VCF output (-Ov/-Oz)");
        app.add_option("--nthreads", n_threads, "Number of threads for compression");

        app.add_flag("--shapeit5", shapeit5_format, "EXPERIMENTAL - Generate an XSI file with SHAPEIT5 internal format");
        app.add_flag("--GP", compress_gp, "EXPERIMENTAL - Compress GP array");
    }
    CLI::App app{"xSqueezeIt - VCF/BCF Compressor"};

    std::string filename = "-";
    std::string ofname = "-";
    bool fast_pipe = false;
    //char O = 'u';
    std::string output_type = "b"; // CL11 doesn't convert to char (runtime e.g., with -O v)
    bool compress = false;
    bool verbose = false;
    bool decompress = false;
    bool info = false;
    bool wait = false;
    bool verify = false;
    bool iota = false;
    bool v4 = true;
    bool zstd = false;
    int  zstd_compression_level = 7; // Some acceptable default value
    double maf = 0.001;
    size_t reset_sort_block_length = 8192;
    bool no_sort = false;
    bool count_xcf = false;
    bool sandbox = false;
    bool compare_matrices = false;
    bool copy_bcf = false;
    bool compute_phase_switch_errors = false;
    bool inject_phase_switches = false;
    double phase_switch_prob = 0.1;
    bool create_map = false;
    bool phase = false;
    bool phase_2 = false;
    bool het_info = false;
    bool unphase = false;
    bool unphase_random = false;
    bool sprinkle_missing = false;
    bool bitmap = false;
    bool bitmap_pbwt = false;
    bool het_bitmap = false;
    bool color_bitmap16 = false;
    bool sorted_bitmap = false;
    bool block_sorted_bitmap = false;
    uint32_t block_size = 32;
    bool histogram_info = false;
    bool partial_pbwt = false;
    bool replace_pseudo = false;
    std::string samples = "";
    std::string samples_file = "";
    std::string regions = "";
    std::string regions_file = "";
    std::string targets = "";
    std::string targets_file = "";
    bool no_header = false;

    size_t n_threads = 0;

    bool shapeit5_format = false;
    bool compress_gp = false;
};

#endif /* __XSQUEEZEIT_HPP__ */