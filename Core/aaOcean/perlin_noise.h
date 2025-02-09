// Perlin Noise Generator
// Author: Amaan Akram 
// https://linkedin.com/in/amaan
// Outputs scalar noise for use in tiling
//
// LICENSE: 
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html
//
// A "New BSD" License for aaOcean can be obtained by contacting the author
// For more details on aaOcean and associated 3rd Party licenses, please see
// license.txt file that is part of the aaOcean repository:
// https://github.com/amaanakram/aaOcean

#ifndef PERLIN_NOISE_H
#define PERLIN_NOISE_H

#include <cmath>
#include <array>
#include <tuple>
#define HAVE_SSE2 1
#include <immintrin.h>  // AVX2 for SIMD optimization
#define SFMT_MEXP 19937
#include "dSFMT/dSFMT.c"

// Rest of the implementation remains exactly the same as before...
class PerlinNoise {
private:
    // Permutation table
    std::array<int, 512> p;
    dsfmt_t sfmt;  // SFMT state
    u_int32_t m_seed;
    
    // 3D gradient vectors (12 edges of a cube)
    static constexpr std::array<std::tuple<float, float, float>, 16> grads3D = {{
        {1,1,0}, {-1,1,0}, {1,-1,0}, {-1,-1,0},
        {1,0,1}, {-1,0,1}, {1,0,-1}, {-1,0,-1},
        {0,1,1}, {0,-1,1}, {0,1,-1}, {0,-1,-1},
        {1,1,0}, {0,-1,1}, {-1,1,0}, {0,-1,-1}
    }};
    
    // 2D gradient vectors for backward compatibility
    static constexpr std::array<std::pair<float, float>, 8> grads2D = {{
        {1.0f, 1.0f}, {-1.0f, 1.0f}, {1.0f, -1.0f}, {-1.0f, -1.0f},
        {1.0f, 0.0f}, {-1.0f, 0.0f}, {0.0f, 1.0f}, {0.0f, -1.0f}
    }};

    // Fade function (improved with 6t^5 - 15t^4 + 10t^3)
    static float fade(float t) 
    {
        return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    }
    
    // Linear interpolation
    static float lerp(float t, float a, float b) 
    {
        return a + t * (b - a);
    }
    
    // 2D gradient function
    float grad2D(int hash, float u, float v) const 
    {
        const auto& gradient = grads2D[hash & 7];
        return gradient.first * u + gradient.second * v;
    }
    
    // 3D gradient function
    float grad3D(int hash, float x, float y, float z) const 
    {
        const auto& gradient = grads3D[hash & 15];
        return std::get<0>(gradient) * x + 
               std::get<1>(gradient) * y + 
               std::get<2>(gradient) * z;
    }

public:
    // Constructor with optional seed
    explicit PerlinNoise(uint32_t seed = 42) 
    {
        m_seed = seed;

        // Initialize permutation table with values 0-255
        for (int i = 0; i < 256; ++i) 
            p[i] = i;
        
        
        // Create random engine with seed
        dsfmt_init_gen_rand(&sfmt, m_seed);
        
        // Shuffle permutation table
        for (int i = 255; i > 0; --i) 
        {
            int j = dsfmt_genrand_uint32(&sfmt) % (i + 1);
            std::swap(p[i], p[j]);
        }
        
        // Duplicate permutation table to avoid overflow
        for (int i = 0; i < 256; ++i) 
            p[256 + i] = p[i];
    }
    
    // 2D noise function
    float noise(float u, float v) const 
    {
        // Find unit square that contains point
        int u_int = static_cast<int>(std::floor(u));
        int v_int = static_cast<int>(std::floor(v));
        
        // Get relative coordinates within unit square
        u -= std::floor(u);
        v -= std::floor(v);
        
        // Compute fade curves
        float u_fade = fade(u);
        float v_fade = fade(v);
        
        // Hash coordinates of cube corners
        int a = p[p[u_int & 255] + (v_int & 255)];
        int b = p[p[(u_int + 1) & 255] + (v_int & 255)];
        int c = p[p[u_int & 255] + ((v_int + 1) & 255)];
        int d = p[p[(u_int + 1) & 255] + ((v_int + 1) & 255)];
        
        // Calculate dot products with gradients
        float x1 = lerp(u_fade,
            grad2D(a, u, v),
            grad2D(b, u - 1, v));
        float x2 = lerp(u_fade,
            grad2D(c, u, v - 1),
            grad2D(d, u - 1, v - 1));
        
        // Interpolate along y
        return lerp(v_fade, x1, x2);
    }
    
    // 3D noise function
    float noise(float x, float y, float z) const 
    {
        // Find unit cube that contains point
        int X = static_cast<int>(std::floor(x)) & 255;
        int Y = static_cast<int>(std::floor(y)) & 255;
        int Z = static_cast<int>(std::floor(z)) & 255;
        
        // Get relative coordinates within unit cube
        x -= std::floor(x);
        y -= std::floor(y);
        z -= std::floor(z);
        
        // Compute fade curves
        float u = fade(x);
        float v = fade(y);
        float w = fade(z);
        
        // Hash coordinates of cube corners
        int A  = p[X] + Y;
        int AA = p[A] + Z;
        int AB = p[A + 1] + Z;
        int B  = p[X + 1] + Y;
        int BA = p[B] + Z;
        int BB = p[B + 1] + Z;
        
        // Calculate dot products with gradients
        float gradAA  = grad3D(p[AA], x, y, z);
        float gradBA  = grad3D(p[BA], x-1, y, z);
        float gradAB  = grad3D(p[AB], x, y-1, z);
        float gradBB  = grad3D(p[BB], x-1, y-1, z);
        float gradAA1 = grad3D(p[AA+1], x, y, z-1);
        float gradBA1 = grad3D(p[BA+1], x-1, y, z-1);
        float gradAB1 = grad3D(p[AB+1], x, y-1, z-1);
        float gradBB1 = grad3D(p[BB+1], x-1, y-1, z-1);
        
        // Trilinear interpolation
        return  lerp(w,
                lerp(v,
                    lerp(u, gradAA, gradBA),
                    lerp(u, gradAB, gradBB)
                ),
                lerp(v,
                    lerp(u, gradAA1, gradBA1),
                    lerp(u, gradAB1, gradBB1)
                )
        );
    }
    
    // Normalized noise functions (output in [0, 1] range)
    float normalizedNoise(float u, float v) const 
    {
        return (noise(u, v) + 1.0f) * 0.5f;
    }
    
    float normalizedNoise(float x, float y, float z) const 
    {
        return (noise(x, y, z) + 1.0f) * 0.5f;
    }

    u_int32_t getSeed()
    {
        return m_seed;
    }
    
    // 2D octave noise
    /*float octaveNoise(
        float u, 
        float v, 
        int octaves, 
        float frequency = 1.0f,
        float persistence = 0.7f,
        float lacunarity = 2.0f
    ) const {
        float total = 0.0f;
        float amplitude = 1.0f;
        float maxValue = 0.0f;
        
        for (int i = 0; i < octaves; ++i) {
            total += noise(u * frequency, v * frequency) * amplitude;
            maxValue += amplitude;
            amplitude *= persistence;
            frequency *= lacunarity;  // Use lacunarity parameter instead of hardcoded 2.0f
        }
        
        return total / maxValue;
    }*/
    
    // 3D octave noise
    float octaveNoise(
        float x, 
        float y, 
        float z, 
        int octaves, 
        float frequency = 1.0f,
        float persistence = 0.7f,
        float lacunarity = 2.0f
    ) const 
    {
        float total = 0.0f;
        float amplitude = 1.0f;
        float maxValue = 0.0f;
        
        for (int i = 0; i < octaves; ++i) 
        {
            total += noise(x * frequency, y * frequency, z * frequency) * amplitude;
            maxValue += amplitude;
            amplitude *= persistence;
            frequency *= lacunarity;  
        }
        
        return total / maxValue;
    }
};

#endif //PERLIN_NOISE_H