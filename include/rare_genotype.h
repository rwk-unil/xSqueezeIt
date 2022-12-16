/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
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

#ifndef _RARE_GENOTYPE_H
#define _RARE_GENOTYPE_H

#define SETBIT(n,i)	(n)|=(1UL<<(i));
#define CLRBIT(n,i)	(n)&=~(1UL<<(i));
#define GETBIT(n,i)	(((n)>>(i))&1U);


class rare_genotype {
public:
	/** @note bit order is implementation dependent, so not the best idea to write a file */
	unsigned int idx : 27;	//Individual index of the non Major/Major genotype
	unsigned int het : 1;	//Is it het?
	unsigned int mis : 1;	//Is it missing?	If none of the two, it is then Minor/Minor
	unsigned int al0 : 1;	//What's the allele on haplotype 0?
	unsigned int al1 : 1;	//What's the allele on haplotype 1?
	unsigned int pha : 1;	//Is the genotype phased?

	rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
	}

	rare_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1, bool _pha) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1; pha = _pha;
	}

	~rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
	}

	unsigned int get() {
		unsigned int value = (idx << 5);
		if (het) SETBIT(value, 4);
		if (mis) SETBIT(value, 3);
		if (al0) SETBIT(value, 2);
		if (al1) SETBIT(value, 1);
		if (pha) SETBIT(value, 0);
		return value;
	}

	void set(unsigned int value) {
		idx = (value >> 5);
		het = GETBIT(value, 4);
		mis = GETBIT(value, 3);
		al0 = GETBIT(value, 2);
		al1 = GETBIT(value, 1);
		pha = GETBIT(value, 0);
	}
};

class rare_genotype2 {
public:
	uint32_t idx : 27;
	uint32_t al0 : 2; /* In BCF format (without phase bit) = bits 1 and 2 of BCF "int" */
	uint32_t al1 : 2; /* In BCF format (without phase bit) = bits 1 and 2 of BCF "int" */
	uint32_t pha : 1;

	/**
	 * @brief Convert BCF "int" allele to 2 bit representation (bit 1 and 2 only of bits 0 to 31)
	 *        Therefore can only handle "missing", "ref", "alt1", "alt2" alleles types
	 *        Only "missing", "ref", and "alt" ("alt1") are used here
	 *        Note: The 32-bit BCF format has bit 0 for phasing and can represent more alts
	 *
	 * @param bcf_binary_allele
	 * @return uint32_t
	 */
	inline uint32_t bcf_binary_allele_to_2bit(int32_t bcf_binary_allele) {
		return (bcf_binary_allele >> 1) & 0x3;
	}

	rare_genotype2() :
		idx(0), al0(0), al1(0), pha(0)
	{}

	rare_genotype2(uint32_t idx, int32_t al0, int32_t al1) :
		idx(idx), al0(bcf_binary_allele_to_2bit(al0)), al1(bcf_binary_allele_to_2bit(al1)),
		pha(bcf_gt_is_phased(al0) | bcf_gt_is_phased(al1)) /* usually the bit is set only on allele 1 (see VCF spec) */
	{}

	/**
	 * @brief Serializes to uint32_t (careful endianness sensitive ! (note that XSI checks endianness))
	 *
	 * @return uint32_t
	 */
	uint32_t get() {
		/* Serialize */
		uint32_t val = idx;
		val = (val << 2) | (al0 & 0x3);
		val = (val << 2) | (al1 & 0x3);
		val = (val << 1) | (pha & 0x1);

		return val;
	}

	/**
	 * @brief Deserializes from uint32_t (careful endianness sensitive ! (note that XSI checks endianness))
	 *
	 * @param val
	 */
	void set(uint32_t val) {
		/* Deserialize (mirror of the above) */
		pha = val & 0x1;
		al1 = (val >> 1) & 0x3;
		al0 = (val >> 3) & 0x3;
		idx = (val >> 5) & 0x07FFFFFF;
	}

	/**
	 * @brief Get allele 0 as BCF "int" representation (handles missing, and has phase bit)
	 *
	 * @return int32_t
	 */
	inline int32_t get_allele0_bcf() {
		return int32_t(al0) << 1;
	}

	/**
	 * @brief Get allele 1 as BCF "int" representation (handles missing, and has phase bit)
	 *
	 * @return int32_t
	 */
	inline int32_t get_allele1_bcf() {
		return (int32_t(al1) << 1) | pha;
	}
};

#endif
