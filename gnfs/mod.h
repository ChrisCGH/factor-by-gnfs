#ifndef _MOD_H
#define _MOD_H
#ifndef WIN32
#define DEFINEDWIN32
#define WIN32
#endif
#ifndef WIN32
/*
#define modasm1(ll, p, r) \
{ \
   r = static_cast<long int>(ll % p); \
   if (r < 0) r += p; \
}
*/
#define modasm1(ll, p, r) \
{ \
   long int a = ll >> 32; \
   unsigned long int b = (unsigned long int)(ll & 0xffffffffLL); \
   asm( \
       "xorl	%%eax, %%eax\n\t" \
       "movl	$1, %%edx\n\t" \
       "divl	%[" #p "]\n\t" \
       "movl	%%edx, %%eax\n\t" \
       "imull	%[a]\n\t" \
       "addl	%[b], %%eax\n\t" \
       "adcl	$0, %%edx\n\t" \
       "idivl	%[" #p "]\n\t" \
       "testl	%%edx, %%edx\n\t" \
       "jns	1f\n\t" \
       "addl	%[" #p "], %%edx\n" \
       "1:\n\t" \
       "movl	%%edx, %[" #r "]" \
   : [r] "=&r" (r) \
   : [p] "r" (p), [a] "r" (a), [b] "r" (b) \
   : "dx", "ax"); \
}
#else
#define modasm1(ll, p, r) \
{ \
   r = static_cast<long int>(ll % p); \
   if (r < 0) r += p; \
}
#endif
#ifndef WIN32
#define modasm2(ll, p, r) \
long int r; \
{ \
   long int a = ll >> 32; \
   unsigned long int b = (unsigned long int)(ll & 0xffffffffLL); \
   if (a >= p || p < 100L) \
   { \
   asm( \
       "xorl	%%eax, %%eax\n\t" \
       "movl	$1, %%edx\n\t" \
       "divl	%[" #p "]\n\t" \
       "movl	%%edx, %%eax\n\t" \
       "imull	%[a]\n\t" \
       "addl	%[b], %%eax\n\t" \
       "adcl	$0, %%edx\n\t" \
       "divl	%[" #p "]\n\t" \
       "movl	%%edx, %[" #r "]" \
   : [r] "=&r" (r) \
   : [p] "r" (p), [a] "r" (a), [b] "r" (b) \
   : "dx", "ax"); \
   } \
   else \
   { \
   asm( \
       "divl	%[" #p "]\n\t" \
   : [r] "=&d" (r) \
   : [p] "r" (p), [a] "0" (a), [b] "a" (b) ); \
   } \
}
#else
#define modasm2(ll, p, r) \
long int r; \
{ \
   r = static_cast<long int>(ll % p); \
   if (r < 0) r += p; \
}
#endif

inline long int modasm(long long int ll, long int p)
{
   long int r;
   modasm1(ll, p, r);
   return r;
}

#ifndef WIN32
#define mulmodasm1(l1, l2, p, r) \
long int r; \
{ \
   asm( \
       "imul	%[" #l1 "]\n\t" \
       "idivl	%[" #p "]\n\t" \
   : [r] "=&d" (r) \
   : [l1] "r" (l1), [l2] "a" (l2), [p] "r" (p) ); \
}
#else
#define mulmodasm1(l1, l2, p, r) \
long int r; \
{ \
   r = ((long long)l1 * l2) % p; \
}
#endif

#ifndef WIN32
#define mulmodasm2(l1, l2, p, r) \
{ \
   asm( \
       "imul	%[" #l1 "]\n\t" \
       "idivl	%[" #p "]\n\t" \
       "testl	%%edx, %%edx\n\t" \
       "jns	1f\n\t" \
       "addl	%[" #p "], %%edx\n" \
       "1:\n\t" \
   : [r] "=&d" (r) \
   : [l1] "r" (l1), [l2] "a" (l2), [p] "r" (p) ); \
}
#else
#define mulmodasm2(l1, l2, p, r) \
{ \
   r = static_cast<long int>(((long long)l1 * l2) % p); \
   if (r < 0) r += p; \
}
#endif

#ifndef WIN32
#define mulmodasm3(l1, l2, p, r) \
long int r; \
{ \
   asm( \
       "imul	%[" #l1 "]\n\t" \
       "idivl	%[" #p "]\n\t" \
   : [r] "=&d" (r), [l2] "=&a" (l2) \
   : [l1] "r" (l1), [p] "r" (p) ); \
}
#else
#define mulmodasm3(l1, l2, p, r) \
long int r; \
{ \
   r = ((long long)l1 * l2) % p; \
   l2 = ((long long)l1 * l2) / p; \
}
#endif

inline long int mulmodasm(long int l1, long int l2, long int p)
{
   long int r;
   mulmodasm2(l1, l2, p, r);
   return r;
}
#ifdef DEFINEDWIN32
#undef WIN32
#endif
#endif
