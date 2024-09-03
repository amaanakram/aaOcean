#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623158e+308
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.402823466e+38F
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define CLAMP_TO_RANGE 1

static const float aa_GRAVITY       = 9.80665f;
static const float aa_GRAVITYSQ		= aa_GRAVITY * aa_GRAVITY;
static const float aa_EPSILON       = 0.0000001f;
static const float aa_PI            = 3.14159265358979323846f;
static const float aa_TWOPI         = aa_PI * 2.0f;
static const float aa_PIBYTWO       = 1.57079632679489661923f;
static const float aa_180BYPI       = 57.295779513082320876798154814105f;
static const float aa_PIBY180       = 0.01745329251994329576922222222222f;
static const float aa_INV_PIBYTWO   = 0.63661977236758134307607071493546f;
static const float aa_INV_SQRTTWO   = 0.70710678118654752440084436210485f;
static const float aa_FLT_EPSILON   = 1.0e-6f;

#define BOUNDARY 16 // alignment boundary

#ifdef _MSC_VER
#define ALIGN(x) __declspec(align(x))
#else // gcc
#define ALIGN(x) __attribute__((aligned(x)))
#endif

#ifdef _MSC_VER
#define f_inline __forceinline
#else // gcc
#define f_inline inline
#endif

#endif /*CONSTANTS_H*/
