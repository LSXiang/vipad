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

#include <math.h>
#include <stdio.h>

#include "matd.h"
#include "zarray.h"
#include "homography.h"
#include "math_util.h"

namespace apriltag {

// correspondences is a list of float[4]s, consisting of the points x
// and y concatenated. We will compute a homography such that y = Hx
matd_t *homography_compute(zarray_t *correspondences)
{
    // compute centroids of both sets of points (yields a better
    // conditioned information matrix)
    double x_cx = 0, x_cy = 0;
    double y_cx = 0, y_cy = 0;

    for (int i = 0; i < zarray_size(correspondences); i++) {
        float *c;
        zarray_get_volatile(correspondences, i, &c);

        x_cx += c[0];
        x_cy += c[1];
        y_cx += c[2];
        y_cy += c[3];
    }

    int sz = zarray_size(correspondences);
    x_cx /= sz;
    x_cy /= sz;
    y_cx /= sz;
    y_cy /= sz;

    // NB We don't normalize scale; it seems implausible that it could
    // possibly make any difference given the dynamic range of IEEE
    // doubles.

    matd_t *A = matd_create(9,9);
    for (int i = 0; i < zarray_size(correspondences); i++) {
        float *c;
        zarray_get_volatile(correspondences, i, &c);

        // (below world is "x", and image is "y")
        double worldx = c[0] - x_cx;
        double worldy = c[1] - x_cy;
        double imagex = c[2] - y_cx;
        double imagey = c[3] - y_cy;

        double a03 = -worldx;
        double a04 = -worldy;
        double a05 = -1;
        double a06 = worldx*imagey;
        double a07 = worldy*imagey;
        double a08 = imagey;

        MATD_EL(A, 3, 3) += a03*a03;
        MATD_EL(A, 3, 4) += a03*a04;
        MATD_EL(A, 3, 5) += a03*a05;
        MATD_EL(A, 3, 6) += a03*a06;
        MATD_EL(A, 3, 7) += a03*a07;
        MATD_EL(A, 3, 8) += a03*a08;
        MATD_EL(A, 4, 4) += a04*a04;
        MATD_EL(A, 4, 5) += a04*a05;
        MATD_EL(A, 4, 6) += a04*a06;
        MATD_EL(A, 4, 7) += a04*a07;
        MATD_EL(A, 4, 8) += a04*a08;
        MATD_EL(A, 5, 5) += a05*a05;
        MATD_EL(A, 5, 6) += a05*a06;
        MATD_EL(A, 5, 7) += a05*a07;
        MATD_EL(A, 5, 8) += a05*a08;
        MATD_EL(A, 6, 6) += a06*a06;
        MATD_EL(A, 6, 7) += a06*a07;
        MATD_EL(A, 6, 8) += a06*a08;
        MATD_EL(A, 7, 7) += a07*a07;
        MATD_EL(A, 7, 8) += a07*a08;
        MATD_EL(A, 8, 8) += a08*a08;

        double a10 = worldx;
        double a11 = worldy;
        double a12 = 1;
        double a16 = -worldx*imagex;
        double a17 = -worldy*imagex;
        double a18 = -imagex;

        MATD_EL(A, 0, 0) += a10*a10;
        MATD_EL(A, 0, 1) += a10*a11;
        MATD_EL(A, 0, 2) += a10*a12;
        MATD_EL(A, 0, 6) += a10*a16;
        MATD_EL(A, 0, 7) += a10*a17;
        MATD_EL(A, 0, 8) += a10*a18;
        MATD_EL(A, 1, 1) += a11*a11;
        MATD_EL(A, 1, 2) += a11*a12;
        MATD_EL(A, 1, 6) += a11*a16;
        MATD_EL(A, 1, 7) += a11*a17;
        MATD_EL(A, 1, 8) += a11*a18;
        MATD_EL(A, 2, 2) += a12*a12;
        MATD_EL(A, 2, 6) += a12*a16;
        MATD_EL(A, 2, 7) += a12*a17;
        MATD_EL(A, 2, 8) += a12*a18;
        MATD_EL(A, 6, 6) += a16*a16;
        MATD_EL(A, 6, 7) += a16*a17;
        MATD_EL(A, 6, 8) += a16*a18;
        MATD_EL(A, 7, 7) += a17*a17;
        MATD_EL(A, 7, 8) += a17*a18;
        MATD_EL(A, 8, 8) += a18*a18;

        double a20 = -worldx*imagey;
        double a21 = -worldy*imagey;
        double a22 = -imagey;
        double a23 = worldx*imagex;
        double a24 = worldy*imagex;
        double a25 = imagex;

        MATD_EL(A, 0, 0) += a20*a20;
        MATD_EL(A, 0, 1) += a20*a21;
        MATD_EL(A, 0, 2) += a20*a22;
        MATD_EL(A, 0, 3) += a20*a23;
        MATD_EL(A, 0, 4) += a20*a24;
        MATD_EL(A, 0, 5) += a20*a25;
        MATD_EL(A, 1, 1) += a21*a21;
        MATD_EL(A, 1, 2) += a21*a22;
        MATD_EL(A, 1, 3) += a21*a23;
        MATD_EL(A, 1, 4) += a21*a24;
        MATD_EL(A, 1, 5) += a21*a25;
        MATD_EL(A, 2, 2) += a22*a22;
        MATD_EL(A, 2, 3) += a22*a23;
        MATD_EL(A, 2, 4) += a22*a24;
        MATD_EL(A, 2, 5) += a22*a25;
        MATD_EL(A, 3, 3) += a23*a23;
        MATD_EL(A, 3, 4) += a23*a24;
        MATD_EL(A, 3, 5) += a23*a25;
        MATD_EL(A, 4, 4) += a24*a24;
        MATD_EL(A, 4, 5) += a24*a25;
        MATD_EL(A, 5, 5) += a25*a25;
    }

    // make symmetric
    for (int i = 0; i < 9; i++)
        for (int j = i+1; j < 9; j++)
            MATD_EL(A, j, i) = MATD_EL(A, i, j);

    matd_t *H = matd_create(3,3);

    // compute singular vector by (carefully) inverting the rank-deficient matrix.
    double data[] = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    matd_t *b = matd_create_data(9, 1, data);
    matd_t *Ainv = NULL;

    if (0) {
        matd_plu_t *lu = matd_plu(A);
        Ainv = matd_plu_solve(lu, b);
        matd_plu_destroy(lu);
    } else {
        matd_chol_t *chol = matd_chol(A);
        Ainv = matd_chol_solve(chol, b);
        matd_chol_destroy(chol);
    }

    double scale = 0;

    for (int i = 0; i < 9; i++)
        scale += sq(MATD_EL(Ainv, i, 0));
    scale = sqrt(scale);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            MATD_EL(H, i, j) = MATD_EL(Ainv, 3*i+j, 0) / scale;

    matd_destroy(b);
    matd_destroy(Ainv);


    matd_t *Tx = matd_identity(3);
    MATD_EL(Tx,0,2) = -x_cx;
    MATD_EL(Tx,1,2) = -x_cy;

    matd_t *Ty = matd_identity(3);
    MATD_EL(Ty,0,2) = y_cx;
    MATD_EL(Ty,1,2) = y_cy;

    matd_t *H2 = matd_op("M*M*M", Ty, H, Tx);

    matd_destroy(A);
    matd_destroy(Tx);
    matd_destroy(Ty);
    matd_destroy(H);

    return H2;
}

} /* namespace apriltag */

