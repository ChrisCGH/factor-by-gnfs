#ifndef BITOPERATIONS_H
#define BITOPERATIONS_H

#include <stdint.h>

struct BitOperations
{
    enum { BITS_IN_WORD = 32 };
    static inline uint32_t colMask(size_t cols)
    {
        if (cols >= BITS_IN_WORD)
        {
            return 0xFFFFFFFF;
        }
        return (1U << cols) - 1U;
    }

    static inline bool bitSet(size_t i, uint32_t w)
    {
        // returns zero if bit i is not set in w
        // returns non-zero otherwise
        // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
        // where 0x00000001 has bit 0 set
        if (i >= BITS_IN_WORD)
            return false;
        return (w & (1U << i)) != 0;
    }
    static inline bool bitClear(size_t i, uint32_t w)
    {
        // returns true if bit i is not set in w
        // returns false otherwise
        // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
        // where 0x00000001 has bit 0 set
        if (i >= BITS_IN_WORD)
            return true;
        return (w & (1U << i)) == 0;
    }
    static inline void setBit(size_t i, uint32_t& w)
    {
        // sets bit i (counting from 0) in w
        if (i >= BITS_IN_WORD)
            return;
        w |= (1U << i);
    }
    static inline bool clearBit(size_t i, uint32_t& w)
    {
        // clears bit i (counting from 0) in w
        if (i >= BITS_IN_WORD)
            return false;
        uint32_t ww(w);
        w &= ~(1U << i);
        return (ww != w);
    }
    static inline void copyBit(int bit, size_t i, uint32_t& w)
    {
        if (i >= BITS_IN_WORD)
            return;
        if (!bit)
            w &= ~(1U << i);
        else
            w |= (1U << i);
    }
    static inline void toggleBit(size_t i, uint32_t& w)
    {
        if (i >= BITS_IN_WORD)
            return;
        w ^= (1U << i);
    }
    static int highestBitSetInByte(unsigned char b)
    {
        static int highestBit [256] =
        {   -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
        return highestBit[b];
    }
    static inline int highestBitSet(uint32_t w)
    {
#ifdef __GNUC__
        if (!w) return -1;
        return 31 - __builtin_clz(w);
#else
        if (!w)
        {
            return -1;
        }
        unsigned char byte = static_cast<unsigned char>((w >> 24) & 0xFF);
        if (byte)
        {
            return highestBitSetInByte(byte) + 24;
        }
        byte = static_cast<unsigned char>((w >> 16) & 0xFF);
        if (byte)
        {
            return highestBitSetInByte(byte) + 16;
        }
        byte = static_cast<unsigned char>((w >> 8) & 0xFF);
        if (byte)
        {
            return highestBitSetInByte(byte) + 8;
        }
        byte = static_cast<unsigned char>(w & 0xFF);
        return highestBitSetInByte(byte);
#endif
    }
    static int bitCountInByte(unsigned char b)
    {
        static int bits[256] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
                               };
        return bits[b];
    }
    static inline int bitCount(uint32_t w)
    {
#ifdef __GNUC__
        return __builtin_popcount(w);
#else
        int bit_count = 0;
        for (int i = 0; i < 4; ++i)
        {
            bit_count += bitCountInByte(static_cast<unsigned char>(w & 0xFF));
            w >>= 8;
        }
        return bit_count;
#endif
    }
};

struct BitOperations64
{
    enum { BITS_IN_WORD = 64 };
    static inline unsigned long long int colMask(size_t cols)
    {
        if (cols >= BITS_IN_WORD)
        {
            return 0xFFFFFFFFFFFFFFFFULL;
        }
        return (1ULL << cols) - 1ULL;
    }

    static inline bool bitSet(size_t i, unsigned long long int w)
    {
        // returns zero if bit i is not set in w
        // returns non-zero otherwise
        // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
        // where 0x00000001 has bit 0 set
        if (i >= BITS_IN_WORD)
            return false;
        return (w & (1ULL << i)) != 0;
    }
    static inline bool bitClear(size_t i, unsigned long long int w)
    {
        // returns true if bit i is not set in w
        // returns false otherwise
        // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
        // where 0x00000001 has bit 0 set
        if (i >= BITS_IN_WORD)
            return true;
        return (w & (1ULL << i)) == 0;
    }

    static inline void setBit(size_t i, unsigned long long int& w)
    {
        // sets bit i (counting from 0) in w
        if (i >= BITS_IN_WORD)
            return;
        w |= (1ULL << i);
    }
    static inline void clearBit(size_t i, unsigned long long int& w)
    {
        // clears bit i (counting from 0) in w
        if (i >= BITS_IN_WORD)
            return;
        w &= ~(1ULL << i);
    }
    static inline void copyBit(int bit, size_t i, unsigned long long int& w)
    {
        if (i >= BITS_IN_WORD)
            return;
        if (!bit)
            w &= ~(1ULL << i);
        else
            w |= (1ULL << i);
    }
    static inline void toggleBit(size_t i, unsigned long long int& w)
    {
        if (i >= BITS_IN_WORD)
            return;
        w ^= (1ULL << i);
    }
    static inline int highestBitSet(unsigned long long int w)
    {
#ifdef __GNUC__
        if (!w) return -1;
        return 63 - __builtin_clzll(w);
#else
        if (!w)
        {
            return -1;
        }
        unsigned char byte = static_cast<unsigned char>((w >> 56) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 56;
        }
        byte = static_cast<unsigned char>((w >> 48) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 48;
        }
        byte = static_cast<unsigned char>((w >> 40) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 40;
        }
        byte = static_cast<unsigned char>((w >> 32) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 32;
        }
        byte = static_cast<unsigned char>((w >> 24) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 24;
        }
        byte = static_cast<unsigned char>((w >> 16) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 16;
        }
        byte = static_cast<unsigned char>((w >> 8) & 0xFF);
        if (byte)
        {
            return BitOperations::highestBitSetInByte(byte) + 8;
        }
        byte = static_cast<unsigned char>(w & 0xFF);
        return BitOperations::highestBitSetInByte(byte);
#endif
    }
    static inline int bitCount(unsigned long long int w)
    {
#ifdef __GNUC__
        return __builtin_popcountll(w);
#else
        int bit_count = 0;
        for (int i = 0; i < 8; ++i)
        {
            bit_count += BitOperations::bitCountInByte(static_cast<unsigned char>(w & 0xFF));
            w >>= 8;
        }
        return bit_count;
#endif
    }
};

template <int bit_count>
struct BitArray64
{
    BitArray64() {}
    void set(size_t bit)
    {
        // BITS_IN_WORD = 64 = 2^6, so bit % 64 = bit & 63, bit / 64 = bit >> 6
        BitOperations64::setBit(bit & 63, word[bit >> 6]);
    }
    void clear(size_t bit)
    {
        BitOperations64::clearBit(bit & 63, word[bit >> 6]);
    }
    bool isSet(size_t bit) const
    {
        return BitOperations64::bitSet(bit & 63, word[bit >> 6]);
    }
    void clear()
    {
        for (size_t i = 0; i < bit_count / BitOperations64::BITS_IN_WORD; ++i)
        {
            word[i] = 0ULL;
        }
    }
    void set()
    {
        for (size_t i = 0; i < bit_count / BitOperations64::BITS_IN_WORD; ++i)
        {
            word[i] = ~0ULL;
        }
    }
    unsigned long long int word[bit_count / BitOperations64::BITS_IN_WORD];
};

#endif
