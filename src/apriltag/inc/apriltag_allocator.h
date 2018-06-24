/**
 * Copyright (c) 2018, The Akatsuki(jacob.lsx). All rights reserved.
 * 
 * This software was developed of Akatsuki(jacob.lsx). This software may be
 * available under alternative licensing terms; contact the address above.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of Akatsuki(jacob.lsx).
 */

#ifndef __APRILTAG_ALLOCATOR_H__
#define __APRILTAG_ALLOCATOR_H__

#include <stdlib.h>

namespace apriltag {

static inline void *apriltagMalloc(size_t _size)
{
    return malloc(_size);
}

static inline void *apriltagCalloc(size_t _nmemb, size_t _size)
{
    return calloc(_nmemb, _size);
}

static inline void *apriltagRealloc(void *_ptr, size_t _new_size)
{
    return realloc(_ptr, _new_size);
}

static inline void apriltagFree(void *_ptr)
{
    free(_ptr);
}

} /* namespace apriltag */

#endif /* __APRILTAG_ALLOCATOR_H__ */
