#pragma once
#include <stdint.h>
#include <math.h> //fabsf only

//Windows header removal
#if defined(WIN32)||defined(_WIN32)||defined(__WIN32)&&!defined(__CYGWIN__)
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#endif

#if !defined(DFC_BUILD_STATIC)
//DLL export defines
	#if defined(DFC_BUILD_DLL)
	#define DFC_API __declspec(dllexport)
	#else
	#define DFC_API __declspec(dllimport)
	#endif
#else
//Static export defines
	#if defined(__cplusplus)
	#define DFC_API
	#else
	#define DFC_API extern
	#endif
#endif

#if defined(_MSC_VER)
#define DFC_COMPILE_MSVC
#elif (defined(__GNUC__) && !defined(__clang__))
#define DFC_COMPILE_GCC
#endif

#if defined(DFC_COMPILE_MSVC)
#include <intrin.h>
#endif

#if defined(__cplusplus)
namespace dfc {
#endif

#define dfc_length(arr) (sizeof(arr)/sizeof(*arr))
#define dfc_bitMask(bits) ((1<<(bits))-1)

#define DFC_F32_LITERAL 0x1F
#define DFC_ZERO_LITERAL 0xFC
#define DFC_INF_LITERAL 0xFD
#define DFC_NAN_LITERAL 0xFE

#if defined(__cplusplus)
extern "C" {
#endif

static const uint16_t dfc_tagTable[32] = {
	0,  4,  8,  12, 16, 1,  20, 5,
	24, 9,  28, 13, 2,  17, 21, 6,
	25, 10, 29, 14, 3,  18, 22, 7,
	26, 11, 30, 15, 19, 23, 27, 31
};
static const int dfc_mulTable[4] = {
	0, 0x111111, 0x10101, 0x1001
};

typedef union dfc_floatMem {
	float single;
	uint32_t memory;
	struct {
		uint32_t frac : 23;
		uint32_t exp : 8;
		uint32_t sign : 1;
	};
} dfc_floatMem_t;

typedef union dfc_tag {
	uint8_t value;
	struct {
		uint8_t run : 2;
		uint8_t real : 3;
		uint8_t exp : 3;
	};
} dfc_tag_t;

typedef struct dfc_value {
	uint32_t memory;
	dfc_tag_t tag;
} dfc_value_t;

//Public Functions

/// <summary>
/// Encodes a float to a dynamic float.
/// </summary>
/// <param name="value">: The input float</param>
/// <param name="maxErr">: Max encoding precision error +-</param>
/// <param name="zeroThresh">: The threshold to truncate encoding to zero</param>
/// <returns>The encoded DFC value</returns>
DFC_API dfc_value_t dfc_cvtps_df(const float value, 
								 const float maxErr, 
								 const float zeroThresh);
/// <summary>
/// Decodes a dynamic float to a single precision float.
/// </summary>
/// <param name="value">: The encoded dfc value</param>
/// <returns>The decoded float</returns>
DFC_API float dfc_cvtdf_ps(const dfc_value_t value);
/// <summary>
/// Gets the length of the encoded dfc buffer.
/// </summary>
/// <param name="tag">: The encoder tag</param>
/// <returns>the length in bits</returns>
DFC_API int dfc_getEncodedLength(const dfc_tag_t tag);


#if defined(__cplusplus)
//Private namespace
namespace impl {
#endif

//Compiler specific implementations
//Stride 8 bit = 7, 16 bit = 15, 32 bit = 31
inline void dfc_bitScanReverse(int* const pos, const int value, const int stride);

//Private functions

/// <summary>
/// Packs an exponent value into its smallest binary form.
/// </summary>
/// <param name="exponent">: The float exponent</param>
/// <param name="length">: The length in bits (output)</param>
/// <returns>The encoded buffer</returns>
uint8_t dfc_packExponent(uint8_t exponent, int* const length);
/// <summary>
/// Unpacks an encoded exponent to a float compliant exponent.
/// </summary>
/// <param name="value">: The encoded buffer</param>
/// <param name="length">: The bit length of the buffer</param>
/// <returns>A float exponent</returns>
uint8_t dfc_unpackExponent(uint8_t value, const uint8_t length);

/// <summary>
/// Returns the absolute difference between <paramref name="a"/> and <paramref name="b"/>
/// <remarks>sad : "Sum of Absolute Differences"</remarks>
/// </summary>
/// <param name="a">: First float argument</param>
/// <param name="b">: Second float argument</param>
/// <returns>Absolute difference</returns>
inline float dfc_sad_ps(const float a, const float b);
/// <summary>
/// The absolute of a float.
/// </summary>
/// <param name="value"></param>
/// <returns></returns>
inline float dfc_abs_ps(const float value);

/// <summary>
/// Evaluates the best float fraction storage method
/// </summary>
/// <remarks> nomm: no SIMD </remarks>
/// <param name="value">: The input float binding</param>
/// <param name="maxErr">: The max error +-</param>
/// <returns>The fraction encoding</returns>
dfc_value_t dfc_nomm_encodeFrac(const dfc_floatMem_t value, const float maxErr);

#if defined(__cplusplus)
} //private namespace
} //extern "C"
} //namespace
#endif