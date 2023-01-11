#include "dfc.h"

//Public functions
dfc_value_t dfc_cvtps_df(const float value, 
						 const float maxErr, 
						 const float zeroThresh) {
	//Zero encoding
	if(dfc_abs_ps(value) <= zeroThresh) {
		return (dfc_value_t){0, DFC_ZERO_LITERAL};
	}

	dfc_floatMem_t binding = (dfc_floatMem_t){.single = value};
	//NaN and inf encoding
	if(binding.exp == 255) {
		uint8_t tag = binding.frac ? DFC_NAN_LITERAL : DFC_INF_LITERAL;
		return (dfc_value_t){0, tag};
	}

	//Find best fraction/mantissa
	dfc_value_t encoded = dfc_nomm_encodeFrac(binding, maxErr);
	
	//Encoded as boundary
	if(encoded.tag.value == 0xFF) {
		encoded.tag.value = 0;
		binding.exp += encoded.memory;
	}

	//Pack exponent
	int exponentLength;
	uint8_t packedExponent = dfc_packExponent(binding.exp, &exponentLength);
	encoded.tag.exp = (uint8_t)exponentLength;

	//Save exponent with its signature
	encoded.memory <<= exponentLength + 1;
	encoded.memory |= packedExponent;

	//Save negative sign
	encoded.memory <<= 1;
	encoded.memory |= binding.sign;
	return encoded;
}
float dfc_cvtdf_ps(const dfc_value_t value) {
	dfc_floatMem_t buffer;
	switch(value.tag.value) {
		case DFC_ZERO_LITERAL:
			return 0.0f;
		case DFC_INF_LITERAL:
			buffer.memory = 0x7F800000;
			return buffer.single;
		case DFC_NAN_LITERAL:
			buffer.memory = 0x7FFFFFFF;
			return buffer.single;
		default: //2^0
			buffer.memory = 0x3F800000;
			break;
	}

	//Decode the sign of the float
	uint32_t memory = value.memory;
	buffer.sign = memory & 1;
	memory >>= 1;

	//Decode the exponent
	if(value.tag.exp) {
		buffer.exp = dfc_unpackExponent(memory, value.tag.exp);
		memory >>= value.tag.exp + 1;
	}

	//value has 23 real bits as the float frac
	if((value.tag.value & DFC_F32_LITERAL) == DFC_F32_LITERAL) {
		buffer.frac = memory;
		return buffer.single;
	}

	//Construct real part
	uint8_t realOffset = 23 - value.tag.real;
	uint32_t realMask = dfc_bitMask(value.tag.real);
	buffer.frac = (memory & realMask) << realOffset;
	memory >>= value.tag.real;

	//Construct run part
	memory &= dfc_bitMask(value.tag.run * 4);
	memory *= dfc_mulTable[value.tag.run];
	memory >>= (value.tag.real + 1);
	buffer.frac |= memory & (~realMask);
	
	return buffer.single;
}
int dfc_getEncodedLength(const dfc_tag_t tag) {
	switch (tag.value) {
		case DFC_ZERO_LITERAL:
		case DFC_INF_LITERAL:
		case DFC_NAN_LITERAL:
			return 0;
	}
	//Init with exponent length if its not zero
	int size = (tag.exp + 1) * (tag.exp != 0);
	//if encoding is a float, 23 real else encoded size
	if((tag.value & DFC_F32_LITERAL) == DFC_F32_LITERAL) {
		size += 23;
	} else {
		size += tag.run * 4 + tag.real;
	}
	return size;
}

//Compiler specific implementations
inline void dfc_bitScanReverse(int* const pos, const int value, const int stride) {
#if defined(DFC_COMPILE_MSVC)
	_BitScanReverse((unsigned long*)pos, value);
#elif defined(DFC_COMPILE_GCC)
	*pos = stride - __builtin_clz(value);
#else
	int c = 1<<stride;
	*pos = stride;
	while(1) {
		if(value & c)
			break;
		c >>= 1;
		(*pos)--;
	}
#endif
}

//Private functions
uint8_t dfc_packExponent(uint8_t exponent, int* const length) {
	exponent -= 127;
	if(exponent == 0) {
		*length = 0;
		return 0;
	}
	uint8_t sign = exponent >> 7;
	if(sign) exponent = -exponent;
	dfc_bitScanReverse(length, exponent, 7);
	(*length)++;
	exponent |= sign << (*length);
	return exponent;
}
uint8_t dfc_unpackExponent(uint8_t value, const uint8_t length) {
	if(length == 0) return 127u;
	uint8_t mask = dfc_bitMask(length);
	uint8_t sign = (value >> length) & 1;
	value &= mask;
	if(sign) value = -value;
	return 127+value;
}

inline float dfc_sad_ps(const float a, const float b) {
	//float diff = a - b;
	//uint32_t mem = *(uint32_t*)&diff;
	//mem &= 0x7FFFFFFF;
	//return *(float*)&mem;
	return fabsf(a - b);
}
inline float dfc_abs_ps(const float value) {
	//uint32_t mem = *(uint32_t*)&value;
	//mem &= 0x7FFFFFFF;
	//return *(float*)&mem;
	return fabsf(value);
}

dfc_value_t dfc_nomm_encodeFrac(const dfc_floatMem_t value, const float maxErr) {
	uint32_t baseFrac = value.frac;

	uint32_t mask = 0xFF800000;
	uint32_t valueTemplate = value.memory & mask;
	
	//Test the boundary exponents
	dfc_floatMem_t bound;
	bound.memory = valueTemplate;
	bound.exp++;
	if(dfc_sad_ps(bound.single, value.single) <= maxErr)
		return (dfc_value_t){1, 0xFF};
	bound.exp -= 2;
	if(dfc_sad_ps(bound.single, value.single) <= maxErr)
		return (dfc_value_t){(-1), 0xFF};

	for(int i = 0; i < dfc_length(dfc_tagTable) - 1; i++) {
		dfc_tag_t tag = (dfc_tag_t){.value = dfc_tagTable[i]};

		//Construct the real part
		int realOffset = 23 - tag.real;
		uint32_t realMask = (mask >> tag.real);

		uint32_t evalBuffer = baseFrac & realMask;
		uint32_t encodeBuffer = baseFrac >> realOffset;
	
		//Constructing the run can be expensive
		if(tag.run) {
			int runBits = tag.run * 4;
			int runOffset = (23 - (runBits + tag.real));
			uint32_t runValue = baseFrac >> runOffset;
			runValue &= dfc_bitMask(runBits);
			encodeBuffer |= (runValue << tag.real);

			runValue *= dfc_mulTable[tag.run];
			runValue >>= (tag.real + 1);
			evalBuffer |= runValue & (~realMask);
		}
		//Evaluate the encoded values error
		evalBuffer |= valueTemplate;
		dfc_floatMem_t encoded = (dfc_floatMem_t){.memory = evalBuffer};
		float error = dfc_sad_ps(encoded.single, value.single);

		//Flip the last real bit (if used) and redo evaluation
		if(tag.real) {
			encoded.memory ^= (1 << realOffset);
			float flip = dfc_sad_ps(encoded.single, value.single);

			//update buffer if the flip was better
			if(flip < error) {
				encodeBuffer ^= 1;
				error = flip;
			}
		}

		//If encoding was suitable within maxErr
		if(error <= maxErr) {
			dfc_value_t df;
			df.memory = encodeBuffer;
			df.tag = tag;
			return df;
		}
	}
	//Fallback F32
	dfc_value_t df;
	df.memory = value.frac;
	df.tag.value = DFC_F32_LITERAL;
	return df;
}