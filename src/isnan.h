#ifndef isnan_h
#define isnan_h
#include <stdint.h>

inline int isnan_safe(double d)
{
    union { double d; uint64_t x; } u = { d };
    return (u.x << 1) > 0xff70000000000000ull;
}
inline int isinf_safe(double d)
{
    union { double d; uint64_t x; } u = { d };
    return (u.x << 1) == 0xff70000000000000ull;
}
inline int isnan_ofast(double d)
{
  return (isnan_safe(d) || isinf_safe(d));
}

#endif
