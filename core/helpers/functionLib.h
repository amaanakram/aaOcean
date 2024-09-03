#ifndef FUNCTION_LIB_H
#define FUNCTION_LIB_H

#include <assert.h>
#include "constants.h"
#ifdef __GNUC__ 
#include <stdint.h>
#endif
#include <stdlib.h>

inline float DegsToRads(float degrees)  // Degrees to radians conversion...
{ 
    return(degrees * aa_PIBY180);  
}
inline float RadsToDegs(float rads)     // Radians to degrees conversion...
{ 
    return(rads * aa_180BYPI); 
}

template <class T> const T& maximum ( const T& a, const T& b ) 
{
  return (b<a)?a:b;     
}

template <class T> const T& minimum ( const T& a, const T& b ) 
{
  return (b>a)?a:b;     
}

inline float clamp(float x, float a, float b)
{
    return x < a ? a : (x > b ? b : x);
}

inline bool isfEqual(float x, float y, const float epsilon)
{
   return fabs(x - y) <= epsilon;  // very basic float comparator
}

bool isFlEqual(float x, float y, const float epsilon = aa_FLT_EPSILON)
{
    if (fabs(x - y) <= epsilon) {
        return true;
    }
    // Check relative error
    float relativeError = fabs(x - y) / std::max(fabs(x), fabs(y));
    return relativeError <= epsilon;
}

inline float isEven(int x)
{
    // used in FFT
    // (-1)^(x+y), shift back from centre representation for zero-frequencies to left
    // x, y being coordinates of 2d array holding fft data
    if((x % 2)==0)
        return 1.0f;
    else
        return -1.0f;
}

inline long int_sqrt(long r) // paul bourke
{
   long t,b,c = 0;

   for (b = 0x10000000; b != 0; b >>= 2) 
   {
      t = c + b;
      c >>= 1;
      if (t <= r) 
      {
         r -= t;
         c += b;
      }
   }
   return(c);
}

#if defined(_MSC_VER)
#if _MSC_VER<1900
inline int round(float x)
{
   return  (int)(x > 0.0f ? x + 0.5f : x - 0.5f);
}
#endif
#endif

inline bool isInt(float a)
{
    if(static_cast<int>(a) == a)
        return true;
    else
        return false;
}

inline float rescale(const float& value, const float& oldMin, const float& oldMax, const float& newMin, const float& newMax)
{
  const float oldDistance = oldMax - oldMin;
  const float newDistance = newMax - newMin;

  const float distance = (value - oldMin) / maximum(oldDistance, 0.0001f);
  const float newValue = newMin + (distance * newDistance);

  return newValue;
}

float intpow( float base, int exponent )
{
    float out = 1.f, curpwr = base;
    for( ; exponent > 0; exponent = exponent >> 1)
    {
        if ((exponent & 1) > 0)
            base *= curpwr;
        curpwr *= curpwr;
    }
    return out;
}

float fastPow(float a, float b) 
{
  union 
  {
    float d;
    int x[2];
  } 
  u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline float lerp(float t, float a, float b)
{ 
    // t = 0 retunes a
    // t = 1 returns b
    return a*(1.0f - t) + b * t;
}

bool isAligned(void* data, int alignment = 16)
{
    // check that the alignment is a power of two
    assert((alignment & (alignment-1)) == 0); 
    return ((uintptr_t)data & (alignment-1)) == 0;
}

void* aligned_malloc(unsigned int size, unsigned int alignment = 16)
{
#ifdef __GNUC__ 
    // only supporting aligned malloc for windows for now
    void* mem = malloc(size);
    return mem;
#else
    void* aligned =_aligned_malloc( size, alignment);
    return aligned;
#endif

}

void aligned_free(void* aligned)
{
#ifdef __GNUC__
    // only supporting aligned malloc for windows for now
    free(aligned);
#else
    _aligned_free(aligned);
#endif

}

#if defined(_MSC_VER)
    #define DATA_ALIGN(declaration, alignment) __declspec(align(alignment)) declaration
#elif defined(GCC)
    #define DATA_ALIGN(declaration, alignment) declaration __attribute__ ((aligned (alignment)))
#else
    #define DATA_ALIGN(declaration, alignment)
#endif

#endif /*FUNCTION_LIB_H*/
