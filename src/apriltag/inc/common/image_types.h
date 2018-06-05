#ifndef _IMAGE_TYPES_H
#define _IMAGE_TYPES_H

#include <stdint.h>

// to support conversions between different types, we define all image
// types at once. Type-specific implementations can then #include this
// file, assured that the basic types of each image are known.

typedef struct image_u8 image_u8_t;
struct image_u8
{
    const int32_t width;
    const int32_t height;
    const int32_t stride;

    uint8_t *buf;
};

typedef struct image_u8x3 image_u8x3_t;
struct image_u8x3
{
    const int32_t width;
    const int32_t height;
    const int32_t stride; // bytes per line

    uint8_t *buf;
};

#endif
