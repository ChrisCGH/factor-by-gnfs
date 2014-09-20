#ifndef BITOPERATIONS_H
#define BITOPERATIONS_H

#include <stdint.h>

struct BitOperations
{
   enum { BITS_IN_WORD = 32 };
   static uint32_t colMask(size_t cols)
   {
      static uint32_t mask[BITS_IN_WORD + 1] =
      {
         0x00000000,
         0x00000001, 0x00000003, 0x00000007, 0x0000000F,
         0x0000001F, 0x0000003F, 0x0000007F, 0x000000FF,
         0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF,
         0x00001FFF, 0x00003FFF, 0x00007FFF, 0x0000FFFF,
         0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF,
         0x001FFFFF, 0x003FFFFF, 0x007FFFFF, 0x00FFFFFF,
         0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF,
         0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF
      };
      if (cols > BITS_IN_WORD)
      {
         return 0xFFFFFFFF;
      }
      return mask[cols];
   }

   static bool bitSet(size_t i, uint32_t w)
   {
      // returns zero if bit i is not set in w
      // returns non-zero otherwise
      // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
      // where 0x00000001 has bit 0 set
      if (i >= BITS_IN_WORD)
         return false;
      return (w & (1 << i)) != 0;
   }
   static bool bitClear(size_t i, uint32_t w)
   {
      // returns true if bit i is not set in w
      // returns false otherwise
      // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
      // where 0x00000001 has bit 0 set
      if (i >= BITS_IN_WORD)
         return true;
      return (w & (1 << i)) == 0;
   }
   static void setBit(size_t i, uint32_t& w)
   {
      // sets bit i (counting from 0) in w
      if (i >= BITS_IN_WORD)
         return;
      w |= (1 << i);
   }
   static bool clearBit(size_t i, uint32_t& w)
   {
      // clears bit i (counting from 0) in w
      if (i >= BITS_IN_WORD)
         return false;
      uint32_t ww(w);
      w &= ~(1 << i);
      return (ww != w);
   }
   static void copyBit(int bit, size_t i, uint32_t& w)
   {
#if 0
      if (bit) setBit(i, w);
      else clearBit(i, w);
#else
      if (i >= BITS_IN_WORD)
         return;
      if (!bit)
         // clearBit(i, w);
         w &= ~(1 << i);
      else
         // setBit(i, w);
         w |= (1 << i);
#endif
   }
   static void toggleBit(size_t i, uint32_t& w)
   {
      if (bitSet(i, w))
      {
         clearBit(i, w);
      }
      else
      {
         setBit(i, w);
      }
   }
   static int highestBitSetInByte(unsigned char b)
   {
      static int highestBit [256] =
         { -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
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
   static int highestBitSet(uint32_t w)
   {
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
   static int bitCount(uint32_t w)
   {
      int bit_count = 0;
      for (int i = 0; i < 4; ++i)
      {
         bit_count += bitCountInByte(static_cast<unsigned char>(w & 0xFF));
         w >>= 8;
      }
      return bit_count;
   }
};

struct BitOperations64
{
   enum { BITS_IN_WORD = 64 };
   static unsigned long long int colMask(size_t cols)
   {
      static unsigned long long int mask[BITS_IN_WORD + 1] =
      {
         0x00000000,
         0x00000001, 0x00000003, 0x00000007, 0x0000000F,
         0x0000001F, 0x0000003F, 0x0000007F, 0x000000FF,
         0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF,
         0x00001FFF, 0x00003FFF, 0x00007FFF, 0x0000FFFF,
         0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF,
         0x001FFFFF, 0x003FFFFF, 0x007FFFFF, 0x00FFFFFF,
         0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF,
         0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF,
         0x00000001FFFFFFFFULL, 0x00000003FFFFFFFFULL, 0x00000007FFFFFFFFULL, 0x0000000FFFFFFFFFULL,
         0x0000001FFFFFFFFFULL, 0x0000003FFFFFFFFFULL, 0x0000007FFFFFFFFFULL, 0x000000FFFFFFFFFFULL,
         0x000001FFFFFFFFFFULL, 0x000003FFFFFFFFFFULL, 0x000007FFFFFFFFFFULL, 0x00000FFFFFFFFFFFULL,
         0x00001FFFFFFFFFFFULL, 0x00003FFFFFFFFFFFULL, 0x00007FFFFFFFFFFFULL, 0x0000FFFFFFFFFFFFULL,
         0x0001FFFFFFFFFFFFULL, 0x0003FFFFFFFFFFFFULL, 0x0007FFFFFFFFFFFFULL, 0x000FFFFFFFFFFFFFULL,
         0x001FFFFFFFFFFFFFULL, 0x003FFFFFFFFFFFFFULL, 0x007FFFFFFFFFFFFFULL, 0x00FFFFFFFFFFFFFFULL,
         0x01FFFFFFFFFFFFFFULL, 0x03FFFFFFFFFFFFFFULL, 0x07FFFFFFFFFFFFFFULL, 0x0FFFFFFFFFFFFFFFULL,
         0x1FFFFFFFFFFFFFFFULL, 0x3FFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL
      };
      if (cols > BITS_IN_WORD)
      {
         return 0xFFFFFFFFFFFFFFFFULL;
      }
      return mask[cols];
   }

   static bool bitSet(size_t i, unsigned long long int w)
   {
      // returns zero if bit i is not set in w
      // returns non-zero otherwise
      // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
      // where 0x00000001 has bit 0 set
      if (i >= BITS_IN_WORD)
         return false;
      return (w & (1ULL << i)) != 0;
   }
   static bool bitClear(size_t i, unsigned long long int w)
   {
      // returns true if bit i is not set in w
      // returns false otherwise
      // assumes bits are numbered from 0 up to BITS_IN_WORD - 1
      // where 0x00000001 has bit 0 set
      if (i >= BITS_IN_WORD)
         return true;
      return (w & (1ULL << i)) == 0;
   }

   static void setBit(size_t i, unsigned long long int& w)
   {
      // sets bit i (counting from 0) in w
      if (i >= BITS_IN_WORD)
         return;
      w |= (1ULL << i);
   }
   static void clearBit(size_t i, unsigned long long int& w)
   {
      // clears bit i (counting from 0) in w
      if (i >= BITS_IN_WORD)
         return;
      w &= ~(1ULL << i);
   }
   static void copyBit(int bit, size_t i, unsigned long long int& w)
   {
      if (bit) setBit(i, w);
      else clearBit(i, w);
   }
   static void toggleBit(size_t i, unsigned long long int& w)
   {
      if (bitSet(i, w))
      {
         clearBit(i, w);
      }
      else
      {
         setBit(i, w);
      }
   }
   static int highestBitSet(unsigned long long int w)
   {
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
   }
   static int bitCount(unsigned long long int w)
   {
      int bit_count = 0;
      for (int i = 0; i < 8; ++i)
      {
         bit_count += BitOperations::bitCountInByte(static_cast<unsigned char>(w & 0xFF));
         w >>= 8;
      }
      return bit_count;
   }
};

template <int bit_count>
struct BitArray64
{
    BitArray64() {}
    void set(size_t bit)
    {
        BitOperations64::setBit(bit % BitOperations64::BITS_IN_WORD, word[bit / BitOperations64::BITS_IN_WORD]);
    }
    void clear(size_t bit)
    {
        BitOperations64::clearBit(bit % BitOperations64::BITS_IN_WORD, word[bit / BitOperations64::BITS_IN_WORD]);
    }
    bool isSet(size_t bit) const
    {
        //static bool bitSet(size_t i, unsigned long long int w)
        //return (w & (1ULL << i)) != 0;
        //const unsigned long long int& w = word[bit / BitOperations64::BITS_IN_WORD];
        //size_t i = bit % BitOperations64::BITS_IN_WORD;
        //return (w & (1ULL << i)) != 0;
        //return (word[bit / BitOperations64::BITS_IN_WORD] & (1ULL << (bit % BitOperations64::BITS_IN_WORD))) != 0;
        return BitOperations64::bitSet(bit % BitOperations64::BITS_IN_WORD, word[bit / BitOperations64::BITS_IN_WORD]);
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
