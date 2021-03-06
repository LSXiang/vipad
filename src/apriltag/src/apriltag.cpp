/* Copyright (C) 2013-2016, The Regents of The University of Michigan.
All rights reserved.

This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the Regents of The University of Michigan.
*/

#include "apriltag.h"

#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

// #include "image_u8.h"
#include "zarray.h"
#include "matd.h"
#include "homography.h"
#include "math_util.h"
#include "g2d.h"
#include "apriltag.h"
#include "apriltag_math.h"

namespace apriltag {

#ifndef M_PI
# define M_PI 3.141592653589793238462643383279502884196
#endif

extern zarray_t *apriltag_quad_gradient(apriltag_detector_t *td, image_u8_t *im);
extern zarray_t *apriltag_quad_thresh(apriltag_detector_t *td, image_u8_t *im);

// Regresses a model of the form:
// intensity(x,y) = C0*x + C1*y + CC2
// The J matrix is the:
//    J = [ x1 y1 1 ]
//        [ x2 y2 1 ]
//        [ ...     ]
//  The A matrix is J'J

struct graymodel
{
    float A[3][3];
    float B[3];
    float C[3];
};

void graymodel_init(struct graymodel *gm)
{
    memset(gm, 0, sizeof(struct graymodel));
}

void graymodel_add(struct graymodel *gm, float x, float y, float gray)
{
    // update upper right entries of A = J'J
    gm->A[0][0] += x*x;
    gm->A[0][1] += x*y;
    gm->A[0][2] += x;
    gm->A[1][1] += y*y;
    gm->A[1][2] += y;
    gm->A[2][2] += 1;

    // update B = J'gray
    gm->B[0] += x * gray;
    gm->B[1] += y * gray;
    gm->B[2] += gray;
}

void graymodel_solve(struct graymodel *gm)
{
    mat33_sym_solve((float*) gm->A, gm->B, gm->C);
}

float graymodel_interpolate(struct graymodel *gm, float x, float y)
{
    return gm->C[0]*x + gm->C[1]*y + gm->C[2];
}

struct quick_decode_entry
{
    uint64_t rcode;   // the queried code
    uint16_t id;      // the tag ID (a small integer)
    uint8_t hamming;  // how many errors corrected?
    uint8_t rotation; // number of rotations [0, 3]
};

struct quick_decode
{
    int nentries;
    struct quick_decode_entry *entries;
};

/** if the bits in w were arranged in a d*d grid and that grid was
 * rotated, what would the new bits in w be?
 * The bits are organized like this (for d = 3):
 *
 *  8 7 6       2 5 8      0 1 2
 *  5 4 3  ==>  1 4 7 ==>  3 4 5    (rotate90 applied twice)
 *  2 1 0       0 3 6      6 7 8
 **/
static uint64_t rotate90(uint64_t w, uint32_t d)
{
    uint64_t wr = 0;

    for (int32_t r = d-1; r >=0; r--) {
        for (int32_t c = 0; c < d; c++) {
            int32_t b = r + d*c;

            wr = wr << 1;

            if ((w & (((uint64_t) 1) << b))!=0)
                wr |= 1;
        }
    }

    return wr;
}

void quad_destroy(struct quad *quad)
{
    if (!quad)
        return;

    matd_destroy(quad->H);
    apriltagFree(quad);
}

struct quad *quad_copy(struct quad *quad)
{
    struct quad *q = (struct quad *)apriltagCalloc(1, sizeof(struct quad));
    memcpy(q, quad, sizeof(struct quad));
    if (quad->H)
        q->H = matd_copy(quad->H);
    return q;
}

// void quick_decode_add(struct quick_decode *qd, uint64_t code, int id, int hamming)
// {
//     uint32_t bucket = code % qd->nentries;
// 
//     while (qd->entries[bucket].rcode != UINT64_MAX) {
//         bucket = (bucket + 1) % qd->nentries;
//     }
// 
//     qd->entries[bucket].rcode = code;
//     qd->entries[bucket].id = id;
//     qd->entries[bucket].hamming = hamming;
// }

// void quick_decode_uninit(apriltag_family_t *fam)
// {
//     if (!fam->impl)
//         return;
// 
//     struct quick_decode *qd = (struct quick_decode*) fam->impl;
//     apriltagFree(qd->entries);
//     apriltagFree(qd);
//     fam->impl = NULL;
// }

// void quick_decode_init(apriltag_family_t *family, int maxhamming)
// {
//     assert(family->impl == NULL);
//     assert(family->ncodes < 65535);
// 
//     struct quick_decode *qd = (struct quick_decode *)apriltagCalloc(1, sizeof(struct quick_decode));
//     int capacity = family->ncodes;
// 
//     int nbits = family->d * family->d;
// 
//     if (maxhamming >= 1)
//         capacity += family->ncodes * nbits;
// 
//     if (maxhamming >= 2)
//         capacity += family->ncodes * nbits * (nbits-1);
// 
//     if (maxhamming >= 3)
//         capacity += family->ncodes * nbits * (nbits-1) * (nbits-2);
// 
//     qd->nentries = capacity * 3;
// 
// //    printf("capacity %d, size: %.0f kB\n",
// //           capacity, qd->nentries * sizeof(struct quick_decode_entry) / 1024.0);
// 
//     qd->entries = (quick_decode_entry *)apriltagCalloc(qd->nentries, sizeof(struct quick_decode_entry));
//     if (qd->entries == NULL) {
//         printf("apriltag.c: failed to allocate hamming decode table. Reduce max hamming size.\n");
//         exit(-1);
//     }
// 
//     for (int i = 0; i < qd->nentries; i++)
//         qd->entries[i].rcode = UINT64_MAX;
// 
//     for (int i = 0; i < family->ncodes; i++) {
//         uint64_t code = family->codes[i];
// 
//         // add exact code (hamming = 0)
//         quick_decode_add(qd, code, i, 0);
// 
//         if (maxhamming >= 1) {
//             // add hamming 1
//             for (int j = 0; j < nbits; j++)
//                 quick_decode_add(qd, code ^ (1L << j), i, 1);
//         }
// 
//         if (maxhamming >= 2) {
//             // add hamming 2
//             for (int j = 0; j < nbits; j++)
//                 for (int k = 0; k < j; k++)
//                     quick_decode_add(qd, code ^ (1L << j) ^ (1L << k), i, 2);
//         }
// 
//         if (maxhamming >= 3) {
//             // add hamming 3
//             for (int j = 0; j < nbits; j++)
//                 for (int k = 0; k < j; k++)
//                     for (int m = 0; m < k; m++)
//                         quick_decode_add(qd, code ^ (1L << j) ^ (1L << k) ^ (1L << m), i, 3);
//         }
// 
//         if (maxhamming > 3) {
//             printf("apriltag.c: maxhamming beyond 3 not supported\n");
//         }
//     }
// 
//     family->impl = qd;
// 
//     if (0) {
//         int longest_run = 0;
//         int run = 0;
//         int run_sum = 0;
//         int run_count = 0;
// 
//         // This accounting code doesn't check the last possible run that
//         // occurs at the wrap-around. That's pretty insignificant.
//         for (int i = 0; i < qd->nentries; i++) {
//             if (qd->entries[i].rcode == UINT64_MAX) {
//                 if (run > 0) {
//                     run_sum += run;
//                     run_count ++;
//                 }
//                 run = 0;
//             } else {
//                 run ++;
//                 longest_run = imax(longest_run, run);
//             }
//         }
// 
//         printf("quick decode: longest run: %d, average run %.3f\n", longest_run, 1.0 * run_sum / run_count);
//     }
// }

// http://en.wikipedia.org/wiki/Hamming_weight

//types and constants used in the functions below
//uint64_t is an unsigned 64-bit integer variable type (defined in C99 version of C language)
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

//This is a naive implementation, shown for comparison,
//and to help in understanding the better functions.
//This algorithm uses 24 arithmetic operations (shift, add, and).
int popcount64a(uint64_t x)
{
    x = (x & m1 ) + ((x >>  1) & m1 ); //put count of each  2 bits into those  2 bits
    x = (x & m2 ) + ((x >>  2) & m2 ); //put count of each  4 bits into those  4 bits
    x = (x & m4 ) + ((x >>  4) & m4 ); //put count of each  8 bits into those  8 bits
    x = (x & m8 ) + ((x >>  8) & m8 ); //put count of each 16 bits into those 16 bits
    x = (x & m16) + ((x >> 16) & m16); //put count of each 32 bits into those 32 bits
    x = (x & m32) + ((x >> 32) & m32); //put count of each 64 bits into those 64 bits
    return x;
}

//This uses fewer arithmetic operations than any other known
//implementation on machines with slow multiplication.
//This algorithm uses 17 arithmetic operations.
int popcount64b(uint64_t x)
{
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    x += x >>  8;  //put count of each 16 bits into their lowest 8 bits
    x += x >> 16;  //put count of each 32 bits into their lowest 8 bits
    x += x >> 32;  //put count of each 64 bits into their lowest 8 bits
    return x & 0x7f;
}

//This uses fewer arithmetic operations than any other known
//implementation on machines with fast multiplication.
//This algorithm uses 12 arithmetic operations, one of which is a multiply.
int popcount64c(uint64_t x)
{
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}


// returns an entry with hamming set to 255 if no decode was found.
static void quick_decode_codeword(apriltag_family_t *tf, uint64_t rcode,
                                  struct quick_decode_entry *entry)
{
//     struct quick_decode *qd = (struct quick_decode*) tf->impl;

    int threshold = imax(tf->h - tf->d - 1, 0);

    for (int ridx = 0; ridx < 4; ridx++) {

        for (int i = 0, j = tf->ncodes; i < j; i++) {
            int hamming = popcount64c(tf->codes[i] ^ rcode);
            if(hamming <= threshold) {
                entry->rcode = rcode;
                entry->id = i;
                entry->hamming = hamming;
                entry->rotation = ridx;
                return;
            }
        }

        rcode = rotate90(rcode, tf->d);
    }

    entry->rcode = 0;
    entry->id = 65535;
    entry->hamming = 255;
    entry->rotation = 0;
}

static inline int detection_compare_function(const void *_a, const void *_b)
{
    apriltag_detection_t *a = *(apriltag_detection_t**) _a;
    apriltag_detection_t *b = *(apriltag_detection_t**) _b;

    return a->id - b->id;
}

void apriltag_detector_remove_family(apriltag_detector_t *td, apriltag_family_t *fam)
{
//     quick_decode_uninit(fam);
    zarray_remove_value(td->tag_families, &fam, 0);
}

void apriltag_detector_add_family_bits(apriltag_detector_t *td, apriltag_family_t *fam, int bits_corrected)
{
    zarray_add(td->tag_families, &fam);

//     if (!fam->impl)
//         quick_decode_init(fam, bits_corrected);
}

void apriltag_detector_clear_families(apriltag_detector_t *td)
{
//     for (int i = 0; i < zarray_size(td->tag_families); i++) {
//         apriltag_family_t *fam;
//         zarray_get(td->tag_families, i, &fam);
//         quick_decode_uninit(fam);
//     }
    zarray_clear(td->tag_families);
}

apriltag_detector_t *apriltag_detector_create()
{
    apriltag_detector_t *td = (apriltag_detector_t*) apriltagCalloc(1, sizeof(apriltag_detector_t));

//     td->quad_decimate = 1.0;
//     td->quad_sigma = 0.0;

    td->qtp.max_nmaxima = 10;
    td->qtp.min_cluster_pixels = 5;

    td->qtp.max_line_fit_mse = 10.0;
    td->qtp.critical_rad = 10 * M_PI / 180;
    td->qtp.deglitch = 0;
    td->qtp.min_white_black_diff = 5;

    td->tag_families = zarray_create(sizeof(apriltag_family_t*));

    td->refine_edges = 1;

    return td;
}

void apriltag_detector_destroy(apriltag_detector_t *td)
{
    apriltag_detector_clear_families(td);

    zarray_destroy(td->tag_families);
    apriltagFree(td);
}

// returns non-zero if an error occurs (i.e., H has no inverse)
int quad_update_homographies(struct quad *quad)
{
    zarray_t *correspondences = zarray_create(sizeof(float[4]));

    for (int i = 0; i < 4; i++) {
        float corr[4];

        // At this stage of the pipeline, we have not attempted to decode the
        // quad into an oriented tag. Thus, just act as if the quad is facing
        // "up" with respect to our desired corners. We'll fix the rotation
        // later.
        // [-1, -1], [1, -1], [1, 1], [-1, 1]
        corr[0] = (i==0 || i==3) ? -1 : 1;
        corr[1] = (i==0 || i==1) ? -1 : 1;

        corr[2] = quad->p[i][0];
        corr[3] = quad->p[i][1];

        zarray_add(correspondences, &corr);
    }

    if (quad->H)
        matd_destroy(quad->H);

    // XXX Tunable
    quad->H = homography_compute(correspondences);
    
    zarray_destroy(correspondences);
    
    if (quad->H)
        return 0;

    return -1;
}

// returns the decision margin. Return < 0 if the detection should be rejected.
float quad_decode(apriltag_family_t *family, image_u8_t *im, struct quad *quad, struct quick_decode_entry *entry, image_u8_t *im_samples)
{
    // decode the tag binary contents by sampling the pixel
    // closest to the center of each bit cell.

    int64_t rcode = 0;

    // how wide do we assume the white border is?
    float white_border = 1.0;

    // We will compute a threshold by sampling known white/black cells around this tag.
    // This sampling is achieved by considering a set of samples along lines.
    //
    // coordinates are given in bit coordinates. ([0, fam->d]).
    //
    // { initial x, initial y, delta x, delta y, WHITE=1 }
    float patterns[] = {
        // left white column
        static_cast<float>(0 - white_border / 2.0), 0.5,
        0, 1,
        1,

        // left black column
        static_cast<float>(0 + family->black_border / 2.0), 0.5,
        0, 1,
        0,

        // right white column
        static_cast<float>(2*family->black_border + family->d + white_border / 2.0), .5,
        0, 1,
        1,

        // right black column
        static_cast<float>(2*family->black_border + family->d - family->black_border / 2.0), .5,
        0, 1,
        0,

        // top white row
        0.5, static_cast<float>(-white_border / 2.0),
        1, 0,
        1,

        // top black row
        0.5, static_cast<float>(family->black_border / 2.0),
        1, 0,
        0,

        // bottom white row
        0.5, static_cast<float>(2*family->black_border + family->d + white_border / 2.0),
        1, 0,
        1,

        // bottom black row
        0.5, static_cast<float>(2*family->black_border + family->d - family->black_border / 2.0),
        1, 0,
        0

        // XXX float-counts the corners.
    };

    struct graymodel whitemodel, blackmodel;
    graymodel_init(&whitemodel);
    graymodel_init(&blackmodel);

    for (int pattern_idx = 0; pattern_idx < sizeof(patterns)/(5*sizeof(float)); pattern_idx ++) {
        float *pattern = &patterns[pattern_idx * 5];

        int is_white = pattern[4];

        for (int i = 0; i < 2*family->black_border + family->d; i++) {
            float tagx01 = (pattern[0] + i*pattern[2]) / (2*family->black_border + family->d);
            float tagy01 = (pattern[1] + i*pattern[3]) / (2*family->black_border + family->d);

            float tagx = 2*(tagx01-0.5);
            float tagy = 2*(tagy01-0.5);

            float px, py;
            homography_project(quad->H, tagx, tagy, &px, &py);

            // don't round
            int ix = px;
            int iy = py;
            if (ix < 0 || iy < 0 || ix >= im->width || iy >= im->height)
                continue;

            int v = im->buf[iy*im->stride + ix];

            if (im_samples) {
                im_samples->buf[iy*im_samples->stride + ix] = (1-is_white)*255;
            }

            if (is_white)
                graymodel_add(&whitemodel, tagx, tagy, v);
            else
                graymodel_add(&blackmodel, tagx, tagy, v);
        }
    }

    graymodel_solve(&whitemodel);
    graymodel_solve(&blackmodel);

    // XXX Tunable
    if (graymodel_interpolate(&whitemodel, 0, 0) - graymodel_interpolate(&blackmodel, 0, 0) < 0)
        return -1;

    // compute the average decision margin (how far was each bit from
    // the decision boundary?
    //
    // we score this separately for white and black pixels and return
    // the minimum average threshold for black/white pixels. This is
    // to penalize thresholds that are too close to an extreme.
    float black_score = 0, white_score = 0;
    float black_score_count = 1, white_score_count = 1;

    for (int bitidx = 0; bitidx < family->d * family->d; bitidx++) {
        int bitx = bitidx % family->d;
        int bity = bitidx / family->d;

        float tagx01 = (family->black_border + bitx + 0.5) / (2*family->black_border + family->d);
        float tagy01 = (family->black_border + bity + 0.5) / (2*family->black_border + family->d);

        // scale to [-1, 1]
        float tagx = 2*(tagx01-0.5);
        float tagy = 2*(tagy01-0.5);

        float px, py;
        homography_project(quad->H, tagx, tagy, &px, &py);

        rcode = (rcode << 1);

        // don't round.
        int ix = px;
        int iy = py;

        if (ix < 0 || iy < 0 || ix >= im->width || iy >= im->height)
            continue;

        int v = im->buf[iy*im->stride + ix];

        float thresh = (graymodel_interpolate(&blackmodel, tagx, tagy) + graymodel_interpolate(&whitemodel, tagx, tagy)) / 2.0;
        if (v > thresh) {
            white_score += (v - thresh);
            white_score_count ++;
            rcode |= 1;
        } else {
            black_score += (thresh - v);
            black_score_count ++;
        }

        if (im_samples)
            im_samples->buf[iy*im_samples->stride + ix] = (1 - (rcode & 1)) * 255;
    }

    quick_decode_codeword(family, rcode, entry);

    return fmin(white_score / white_score_count, black_score / black_score_count);
}

// float score_decodability(apriltag_family_t *family, image_u8_t *im, struct quad *quad, void *user)
// {
//     struct quick_decode_entry entry;
// 
//     float decision_margin = quad_decode(family, im, quad, &entry, NULL);
// 
//     // hamming trumps decision margin; maximum value for decision_margin is 255.
//     return decision_margin - entry.hamming*1000;
// }

static void refine_edges(apriltag_detector_t *td, image_u8_t *im_orig, struct quad *quad)
{
    float lines[4][4]; // for each line, [Ex Ey nx ny]

    for (int edge = 0; edge < 4; edge++) {
        int a = edge, b = (edge + 1) & 3; // indices of the end points.

        // compute the normal to the current line estimate
        float nx = quad->p[b][1] - quad->p[a][1];
        float ny = -quad->p[b][0] + quad->p[a][0];
        float mag = sqrt(nx*nx + ny*ny);
        nx /= mag;
        ny /= mag;

        // we will now fit a NEW line by sampling points near
        // our original line that have large gradients. On really big tags,
        // we're willing to sample more to get an even better estimate.
        int nsamples = imax(16, mag / 8); // XXX tunable

        // stats for fitting a line...
        float Mx = 0, My = 0, Mxx = 0, Mxy = 0, Myy = 0, N = 0;

        for (int s = 0; s < nsamples; s++) {
            // compute a point along the line... Note, we're avoiding
            // sampling *right* at the corners, since those points are
            // the least reliable.
            float alpha = (1.0 + s) / (nsamples + 1);
            float x0 = alpha*quad->p[a][0] + (1-alpha)*quad->p[b][0];
            float y0 = alpha*quad->p[a][1] + (1-alpha)*quad->p[b][1];

            // search along the normal to this line, looking at the
            // gradients along the way. We're looking for a strong
            // response.
            float Mn = 0;
            float Mcount = 0;

            // XXX tunable: how far to search?  We want to search far
            // enough that we find the best edge, but not so far that
            // we hit other edges that aren't part of the tag. We
            // shouldn't ever have to search more than quad_decimate,
            // since otherwise we would (ideally) have started our
            // search on another pixel in the first place. Likewise,
            // for very small tags, we don't want the range to be too
            // big.
//             float range = td->quad_decimate + 1;
            float range = 1.0 + 1;

            // XXX tunable step size.
            for (float n = -range; n <= range; n +=  0.25) {
                // Because of the guaranteed winding order of the
                // points in the quad, we will start inside the white
                // portion of the quad and work our way outward.
                //
                // sample to points (x1,y1) and (x2,y2) XXX tunable:
                // how far +/- to look? Small values compute the
                // gradient more precisely, but are more sensitive to
                // noise.
                float grange = 1;
                int x1 = x0 + (n + grange)*nx;
                int y1 = y0 + (n + grange)*ny;
                if (x1 < 0 || x1 >= im_orig->width || y1 < 0 || y1 >= im_orig->height)
                    continue;

                int x2 = x0 + (n - grange)*nx;
                int y2 = y0 + (n - grange)*ny;
                if (x2 < 0 || x2 >= im_orig->width || y2 < 0 || y2 >= im_orig->height)
                    continue;

                int g1 = im_orig->buf[y1*im_orig->stride + x1];
                int g2 = im_orig->buf[y2*im_orig->stride + x2];

                if (g1 < g2) // reject points whose gradient is "backwards". They can only hurt us.
                    continue;

                float weight = (g2 - g1)*(g2 - g1); // XXX tunable. What shape for weight=f(g2-g1)?

                // compute weighted average of the gradient at this point.
                Mn += weight*n;
                Mcount += weight;
            }

            // what was the average point along the line?
            if (Mcount == 0)
                continue;

            float n0 = Mn / Mcount;

            // where is the point along the line?
            float bestx = x0 + n0*nx;
            float besty = y0 + n0*ny;

            // update our line fit statistics
            Mx += bestx;
            My += besty;
            Mxx += bestx*bestx;
            Mxy += bestx*besty;
            Myy += besty*besty;
            N++;
        }

        // fit a line
        float Ex = Mx / N, Ey = My / N;
        float Cxx = Mxx / N - Ex*Ex;
        float Cxy = Mxy / N - Ex*Ey;
        float Cyy = Myy / N - Ey*Ey;

        float normal_theta = .5 * atan2f(-2*Cxy, (Cyy - Cxx));
        nx = cosf(normal_theta);
        ny = sinf(normal_theta);
        lines[edge][0] = Ex;
        lines[edge][1] = Ey;
        lines[edge][2] = nx;
        lines[edge][3] = ny;
    }

    // now refit the corners of the quad
    for (int i = 0; i < 4; i++) {

        // solve for the intersection of lines (i) and (i+1)&3.
        float A00 =  lines[i][3],  A01 = -lines[(i+1)&3][3];
        float A10 =  -lines[i][2],  A11 = lines[(i+1)&3][2];
        float B0 = -lines[i][0] + lines[(i+1)&3][0];
        float B1 = -lines[i][1] + lines[(i+1)&3][1];

        float det = A00 * A11 - A10 * A01;

        // inverse.
        if (fabs(det) > 0.001) {
            // solve
            float W00 = A11 / det, W01 = -A01 / det;

            float L0 = W00*B0 + W01*B1;

            // compute intersection
            quad->p[i][0] = lines[i][0] + L0*A00;
            quad->p[i][1] = lines[i][1] + L0*A10;
        } else {
            // this is a bad sign. We'll just keep the corner we had.
//            printf("bad det: %15f %15f %15f %15f %15f\n", A00, A11, A10, A01, det);
        }
    }
}


void apriltag_detection_destroy(apriltag_detection_t *det)
{
    if (det == NULL)
        return;

    matd_destroy(det->H);
    apriltagFree(det);
}

int prefer_smaller(int pref, float q0, float q1)
{
    if (pref)     // already prefer something? exit.
        return pref;

    if (q0 < q1)
        return -1; // we now prefer q0
    if (q1 < q0)
        return 1; // we now prefer q1

    // no preference
    return 0;
}

zarray_t *apriltag_detector_detect(apriltag_detector_t *td, image_u8_t *im_orig)
{
    if (zarray_size(td->tag_families) == 0) {
        zarray_t *s = zarray_create(sizeof(apriltag_detection_t*));
        printf("apriltag.c: No tag families enabled.");
        return s;
    }

    ///////////////////////////////////////////////////////////
    // Step 1. Detect quads according to requested image decimation
    // and blurring parameters.
//     image_u8_t *quad_im = im_orig;
//     if (td->quad_decimate > 1) {
//         quad_im = image_u8_decimate(im_orig, td->quad_decimate);
//     }
// 
//     if (td->quad_sigma != 0) {
//         // compute a reasonable kernel width by figuring that the
//         // kernel should go out 2 std devs.
//         //
//         // max sigma          ksz
//         // 0.499              1  (disabled)
//         // 0.999              3
//         // 1.499              5
//         // 1.999              7
// 
//         float sigma = fabsf((float) td->quad_sigma);
// 
//         int ksz = 4 * sigma; // 2 std devs in each direction
//         if ((ksz & 1) == 0)
//             ksz++;
// 
//         if (ksz > 1) {
// 
//             if (td->quad_sigma > 0) {
//                 // Apply a blur
//                 image_u8_gaussian_blur(quad_im, sigma, ksz);
//             } else {
//                 // SHARPEN the image by subtracting the low frequency components.
//                 image_u8_t *orig = image_u8_copy(quad_im);
//                 image_u8_gaussian_blur(quad_im, sigma, ksz);
// 
//                 for (int y = 0; y < orig->height; y++) {
//                     for (int x = 0; x < orig->width; x++) {
//                         int vorig = orig->buf[y*orig->stride + x];
//                         int vblur = quad_im->buf[y*quad_im->stride + x];
// 
//                         int v = 2*vorig - vblur;
//                         if (v < 0)
//                             v = 0;
//                         if (v > 255)
//                             v = 255;
// 
//                         quad_im->buf[y*quad_im->stride + x] = (uint8_t) v;
//                     }
//                 }
//                 image_u8_destroy(orig);
//             }
//         }
//     }

//    zarray_t *quads = apriltag_quad_gradient(td, im_orig);
    zarray_t *quads = apriltag_quad_thresh(td, im_orig);

//     // adjust centers of pixels so that they correspond to the
//     // original full-resolution image.
//     if (td->quad_decimate > 1) {
//         for (int i = 0; i < zarray_size(quads); i++) {
//             struct quad *q;
//             zarray_get_volatile(quads, i, &q);
// 
//             for (int i = 0; i < 4; i++) {
//                 q->p[i][0] *= td->quad_decimate;
//                 q->p[i][1] *= td->quad_decimate;
//             }
//         }
//     }
// 
//     if (quad_im != im_orig)
//         image_u8_destroy(quad_im);

    zarray_t *detections = zarray_create(sizeof(apriltag_detection_t*));

    td->nquads = zarray_size(quads);

    ////////////////////////////////////////////////////////////////
    // Step 2. Decode tags from each quad.
    if (1) {
        for (int i = 0; i < zarray_size(quads); i++) {
            struct quad *quad_original;
            zarray_get_volatile(quads, i, &quad_original);
            
            // refine edges is not dependent upon the tag family, thus
            // apply this optimazation BEFORE the other work.
            //if (td->quad_decimate > 1 && td->refine_edges) {
            if (td->refine_edges) {
                refine_edges(td, im_orig, quad_original);
            }
            
            // make sure the homographies are computed...
            if (quad_update_homographies(quad_original))
                continue;
            
            for (int famidx = 0; famidx < zarray_size(td->tag_families); famidx ++) {
                apriltag_family_t *family;
                zarray_get(td->tag_families, famidx, &family);
                
                // since the geometry of tag families can vary, start any
                // optimization process over with the original quad.
                struct quad *quad = quad_copy(quad_original);
                
                struct quick_decode_entry entry;
                
                float decision_margin = quad_decode(family, im_orig, quad, &entry, NULL);
                
                if (entry.hamming < 255 && decision_margin >= 0) {
                    apriltag_detection_t *det = (apriltag_detection_t *)apriltagCalloc(1, sizeof(apriltag_detection_t));
                    
                    det->family = family;
                    det->id = entry.id;
                    det->hamming = entry.hamming;
                    det->decision_margin = decision_margin;
                    
                    float theta = -entry.rotation * M_PI / 2.0;
                    float c = cos(theta), s = sin(theta);
                    
                    // Fix the rotation of our homography to properly orient the tag_families
                    matd_t *R = matd_create(3, 3);
                    MATD_EL(R, 0, 0) = c;
                    MATD_EL(R, 0, 1) = -s;
                    MATD_EL(R, 1, 0) = s;
                    MATD_EL(R, 1, 1) = c;
                    MATD_EL(R, 2, 2) = 1;
                    
                    det->H = matd_op("M*M", quad->H, R);
                    
                    matd_destroy(R);
                    
                    homography_project(det->H, 0, 0, &det->c[0], &det->c[1]);
                    
                    // [-1, -1], [1, -1], [1, 1], [-1, 1], Desired points
                    // [-1, 1], [1, 1], [1, -1], [-1, -1], FLIP Y
                    // adjust the points in det->p so that they correspond to
                    // counter-clockwise around the quad, starting at -1,-1.
                    for (int i = 0; i < 4; i++) {
                        int tcx = (i == 1 || i == 2) ? 1 : -1;
                        int tcy = (i < 2) ? 1 : -1;

                        float p[2];

                        homography_project(det->H, tcx, tcy, &p[0], &p[1]);

                        det->p[i][0] = p[0];
                        det->p[i][1] = p[1];
                    }

                    zarray_add(detections, &det);
                }
                
                quad_destroy(quad);
            }
        }
    }

    ////////////////////////////////////////////////////////////////
    // Step 3. Reconcile detections--- don't report the same tag more
    // than once. (Allow non-overlapping duplicate detections.)
    if (1) {
        zarray_t *poly0 = g2d_polygon_create_zeros(4);
        zarray_t *poly1 = g2d_polygon_create_zeros(4);

        for (int i0 = 0; i0 < zarray_size(detections); i0++) {

            apriltag_detection_t *det0;
            zarray_get(detections, i0, &det0);

            for (int k = 0; k < 4; k++)
                zarray_set(poly0, k, det0->p[k], NULL);

            for (int i1 = i0+1; i1 < zarray_size(detections); i1++) {

                apriltag_detection_t *det1;
                zarray_get(detections, i1, &det1);

                if (det0->id != det1->id || det0->family != det1->family)
                    continue;

                for (int k = 0; k < 4; k++)
                    zarray_set(poly1, k, det1->p[k], NULL);

                if (g2d_polygon_overlaps_polygon(poly0, poly1)) {
                    // the tags overlap. Delete one, keep the other.

                    int pref = 0; // 0 means undecided which one we'll keep.
                    pref = prefer_smaller(pref, det0->hamming, det1->hamming);     // want small hamming
                    pref = prefer_smaller(pref, -det0->decision_margin, -det1->decision_margin);      // want bigger margins

                    // if we STILL don't prefer one detection over the other, then pick
                    // any deterministic criterion.
                    for (int i = 0; i < 4; i++) {
                        pref = prefer_smaller(pref, det0->p[i][0], det1->p[i][0]);
                        pref = prefer_smaller(pref, det0->p[i][1], det1->p[i][1]);
                    }

                    if (pref == 0) {
                        // at this point, we should only be undecided if the tag detections
                        // are *exactly* the same. How would that happen?
                        printf("uh oh, no preference for overlappingdetection\n");
                    }

                    if (pref < 0) {
                        // keep det0, destroy det1
                        apriltag_detection_destroy(det1);
                        zarray_remove_index(detections, i1, 1);
                        i1--; // retry the same index
                        goto retry1;
                    } else {
                        // keep det1, destroy det0
                        apriltag_detection_destroy(det0);
                        zarray_remove_index(detections, i0, 1);
                        i0--; // retry the same index.
                        goto retry0;
                    }
                }

              retry1: ;
            }

          retry0: ;
        }

        zarray_destroy(poly0);
        zarray_destroy(poly1);
    }

    for (int i = 0; i < zarray_size(quads); i++) {
        struct quad *quad;
        zarray_get_volatile(quads, i, &quad);
        matd_destroy(quad->H);
    }

    zarray_destroy(quads);

    zarray_sort(detections, detection_compare_function);

    return detections;
}


// Call this method on each of the tags returned by apriltag_detector_detect
void apriltag_detections_destroy(zarray_t *detections)
{
    for (int i = 0; i < zarray_size(detections); i++) {
        apriltag_detection_t *det;
        zarray_get(detections, i, &det);

        apriltag_detection_destroy(det);
    }

    zarray_destroy(detections);
}

} /* namespace apriltag */


