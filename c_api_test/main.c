/*
 * Very small example to show the modifications to a C application to add xSqueezeIt support through the C API
 *
 * */

#include <stdio.h>

#include "xsqueezeit_export/include/c_api.h"
#include "vcf.h"
#include "synced_bcf_reader.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Please pass file as argument\n");
        exit(-1);
    }

    // Add and allocate the mixed xSqueezeIt + BCF/VCF helper class
    c_xcf *c_xcf_p = c_xcf_new();

    bcf_srs_t *sr = bcf_sr_init();
    if (!bcf_sr_add_reader(sr, argv[1])) {
        printf("Could not load file : %s\n", argv[1]);
        exit(-1);
    }

    // Add the readers to the helper class
    c_xcf_add_readers(c_xcf_p, sr);

    int* genotypes = NULL;
    int ngt_arr = 0;
    int ngt = 0;
    int nset = 0;
    bcf1_t* line = NULL;
    int records = 0;
    const int FIRST_READER = 0;

    // The number of samples for xSqueezeIt is not in the variant file header, so a specific method is needed
    printf("The number of samples in %s is %d\n", argv[1], c_xcf_nsamples(argv[1]));

    while ((nset = bcf_sr_next_line(sr))) {
        line = bcf_sr_get_line(sr, FIRST_READER);

        // original HTSLIB bcf_get_genotypes call, to be replaced by the function with support for xSqueezeIt
        //ngt = bcf_get_genotypes(sr->readers[FIRST_READER].header, line, (void**)&genotypes, &ngt_arr);

        // Call the mixed "xcf" xSqueezeIt + VCF/BCF "get_genotypes" function instead of the bcf one
        ngt = c_xcf_get_genotypes(c_xcf_p, FIRST_READER, sr->readers[FIRST_READER].header, line, (void**)&genotypes, &ngt_arr);

        records++;
    }

    // Deallocate the helper class
    c_xcf_delete(c_xcf_p);

    printf("Exctacted %d records\n", records);

    return 0;
}