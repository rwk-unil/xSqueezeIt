# xSqueezeIt integration and unit tests

## Integration tests

The integration tests are run through `cukinia` (https://github.com/savoirfairelinux/cukinia). They run bash scripts to validate the compression / decompression at the command line level.

### Integration tests done

The integration tests compress and decompress and check the output against the bcftools view command. The files should be exactly the same (except maybe the bcftools command in the header). This is checked by running both results through "diff" (diff does not handle gigantic files because it loads them in memory before doing the comparison...).

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

### Running the integration tests

```shell
# Clone cukinia
git submodule update --init cukinia
# Run the tests
./cukinia/cukinia cukinia_v4.conf
```

## Unit tests

Not yet uploaded.

But they test :

- That WAH encoding / decoding works
- That small functions work
- etc.

If the integration tests are passing, there is a high probability the unit tests are OK, because if some core components are broken the integration tests would fail...