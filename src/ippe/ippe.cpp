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

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "ippe.h"

namespace ippe {


void homographyFromSquarePoints(matd_t *_targetPts, float halfLength, matd_t *_H)
{
    assert(_targetPts != nullptr);
    assert(_targetPts->nrows == 4 && _targetPts->ncols == 2);
    
    if (_H != nullptr) matd_destroy(_H);
    _H = matd_create(3, 3);
    
    double p1x = -MATD_EL(_targetPts, 0, 0);
    double p1y = -MATD_EL(_targetPts, 0, 1);
    
    double p2x = -MATD_EL(_targetPts, 1, 0);
    double p2y = -MATD_EL(_targetPts, 1, 1);
    
    double p3x = -MATD_EL(_targetPts, 2, 0);
    double p3y = -MATD_EL(_targetPts, 2, 1);
    
    double p4x = -MATD_EL(_targetPts, 3, 0);
    double p4y = -MATD_EL(_targetPts, 3, 1);
    
    /* analytic solution */
    double detsInv = -1 / (halfLength * (p1x * p2y - p2x * p1y - p1x * p4y + p2x * p3y - p3x * p2y + p4x * p1y
                                         + p3x * p4y - p4x * p3y));
    
    MATD_EL(_H, 0, 0) = detsInv * (p1x * p3x * p2y - p2x * p3x * p1y - p1x * p4x * p2y + p2x * p4x * p1y 
                                 - p1x * p3x * p4y + p1x * p4x * p3y + p2x * p3x * p4y - p2x * p4x * p3y);
    
    MATD_EL(_H, 0, 1) = detsInv * (p1x * p2x * p3y - p1x * p3x * p2y - p1x * p2x * p4y + p2x * p4x * p1y
                                 + p1x * p3x * p4y - p3x * p4x * p1y - p2x * p4x * p3y + p3x * p4x * p2y);
    
    MATD_EL(_H, 0, 2) = detsInv * halfLength * (p1x * p2x * p3y - p2x * p3x * p1y - p1x * p2x * p4y + p1x * p4x * p2y
                                              - p1x * p4x * p3y + p3x * p4x * p1y + p2x * p3x * p4y - p3x * p4x * p2y);
    
    MATD_EL(_H, 1, 0) = detsInv * (p1x * p2y * p3y - p2x * p1y * p3y - p1x * p2y * p4y + p2x * p1y * p4y
                                 - p3x * p1y * p4y + p4x * p1y * p3y + p3x * p2y * p4y - p4x * p2y * p3y);
    
    MATD_EL(_H, 1, 1) = detsInv * (p2x * p1y * p3y - p3x * p1y * p2y - p1x * p2y * p4y + p4x * p1y * p2y
                                 + p1x * p3y * p4y - p4x * p1y * p3y - p2x * p3y * p4y + p3x * p2y * p4y);
    
    MATD_EL(_H, 1, 2) = detsInv * halfLength * (p1x * p2y * p3y - p3x * p1y * p2y - p2x * p1y * p4y + p4x * p1y * p2y
                                              - p1x * p3y * p4y + p3x * p1y * p4y + p2x * p3y * p4y - p4x * p2y * p3y);
    
    MATD_EL(_H, 2, 0) = -detsInv * (p1x * p3y - p3x * p1y - p1x * p4y - p2x * p3y + p3x * p2y + p4x * p1y + p2x * p4y - p4x * p2y);
    
    MATD_EL(_H, 2, 1) =  detsInv * (p1x * p2y - p2x * p1y - p1x * p3y + p3x * p1y + p2x * p4y - p4x * p2y - p3x * p4y + p4x * p3y);
    
    MATD_EL(_H, 2, 2) = 1.0;
}

void ippeComputeRotations(double j00, double j01, double j10, double j11, double p, double q, matd_t *_R1, matd_t *_R2)
{
    /**
     * Note that it is very hard to understand what is going on here from the code, so if you want to have a clear explanation
     * then please refer to the IPPE paper (Algorithm 1 and its description). Or, you can read IPPE matlab code that is certainly 
     * easier to read.
     */
    if (_R1 != nullptr) matd_destroy(_R1);
    _R1 = matd_create(3, 3);
    
    if (_R2 != nullptr) matd_destroy(_R2);
    _R2 = matd_create(3, 3);
    
    double a00, a01, a10, a11, ata00, ata01, ata11, b00, b01, b10, b11, binv00, binv01, binv10, binv11;
    double rtilde00, rtilde01, rtilde10, rtilde11;
    double rtilde00_2, rtilde01_2, rtilde10_2, rtilde11_2;
    double b0, b1, gamma, dtinv;
    double sp;
    
    /* Finds the rotation Rv that rotates a vector (p, q, 1) to the z axis (0, 0, 1) */
    double s, t, krs0, krs1, krs0_2, krs1_2, costh, sinth;
    
    s = sqrt(p * p + q * q + 1);
    t = sqrt(p * p + q * q);
    costh = 1 / s;
    sinth = sqrt(1 - 1 / (s * s));

    krs0 = p / t;
    krs1 = q / t;
    krs0_2 = krs0 * krs0;
    krs1_2 = krs1 * krs1;
    
    double rv00, rv01, rv02;
    double rv10, rv11, rv12;
    double rv20, rv21, rv22;
    
    rv00 = (costh - 1) * krs0_2 + 1;
    rv01 = krs0 * krs1 * (costh - 1);
    rv02 = krs0 * sinth;
    
    rv10 = krs0 * krs1 * (costh - 1);
    rv11 = (costh - 1) * krs1_2 + 1;
    rv12 = krs1 * sinth;
    
    rv20 = -krs0 * sinth;
    rv21 = -krs1 * sinth;
    rv22 = (costh - 1) * (krs0_2 + krs1_2) + 1;
    
    /* setup the 2x2 SVD decomposition */
    b00 = rv00 - p * rv20;
    b01 = rv01 - p * rv21;
    b10 = rv10 - q * rv20;
    b11 = rv11 - q * rv21;

    dtinv = 1.0 / ((b00 * b11 - b01 * b10));

    binv00 = dtinv * b11;
    binv01 = -dtinv * b01;
    binv10 = -dtinv * b10;
    binv11 = dtinv * b00;

    a00 = binv00 * j00 + binv01 * j10;
    a01 = binv00 * j01 + binv01 * j11;
    a10 = binv10 * j00 + binv11 * j10;
    a11 = binv10 * j01 + binv11 * j11;

    /* compute the largest singular value of A */
    ata00 = a00 * a00 + a01 * a01;
    ata01 = a00 * a10 + a01 * a11;
    ata11 = a10 * a10 + a11 * a11;

    gamma = sqrt(0.5 * (ata00 + ata11 + sqrt((ata00 - ata11) * (ata00 - ata11) + 4.0 * ata01 * ata01)));

    /* reconstruct the full rotation matrices */
    rtilde00 = a00 / gamma;
    rtilde01 = a01 / gamma;
    rtilde10 = a10 / gamma;
    rtilde11 = a11 / gamma;

    rtilde00_2 = rtilde00 * rtilde00;
    rtilde01_2 = rtilde01 * rtilde01;
    rtilde10_2 = rtilde10 * rtilde10;
    rtilde11_2 = rtilde11 * rtilde11;

    b0 = sqrt(-rtilde00_2 - rtilde10_2 + 1);
    b1 = sqrt(-rtilde01_2 - rtilde11_2 + 1);
    sp = (-rtilde00 * rtilde01 - rtilde10 * rtilde11);

    if (sp < 0) {
        b1 = -b1;
    }
    
    /* save results */
    MATD_EL(_R1, 0, 0) = (rtilde00)*rv00 + (rtilde10)*rv01 + (b0)*rv02;
    MATD_EL(_R1, 0, 1) = (rtilde01)*rv00 + (rtilde11)*rv01 + (b1)*rv02;
    MATD_EL(_R1, 0, 2) = (b1 * rtilde10 - b0 * rtilde11) * rv00 + (b0 * rtilde01 - b1 * rtilde00) * rv01
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv02;
    MATD_EL(_R1, 1, 0) = (rtilde00)*rv10 + (rtilde10)*rv11 + (b0)*rv12;
    MATD_EL(_R1, 1, 1) = (rtilde01)*rv10 + (rtilde11)*rv11 + (b1)*rv12;
    MATD_EL(_R1, 1, 2) = (b1 * rtilde10 - b0 * rtilde11) * rv10 + (b0 * rtilde01 - b1 * rtilde00) * rv11
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv12;
    MATD_EL(_R1, 2, 0) = (rtilde00)*rv20 + (rtilde10)*rv21 + (b0)*rv22;
    MATD_EL(_R1, 2, 1) = (rtilde01)*rv20 + (rtilde11)*rv21 + (b1)*rv22;
    MATD_EL(_R1, 2, 2) = (b1 * rtilde10 - b0 * rtilde11) * rv20 + (b0 * rtilde01 - b1 * rtilde00) * rv21
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv22;

    MATD_EL(_R2, 0, 0) = (rtilde00)*rv00 + (rtilde10)*rv01 + (-b0) * rv02;
    MATD_EL(_R2, 0, 1) = (rtilde01)*rv00 + (rtilde11)*rv01 + (-b1) * rv02;
    MATD_EL(_R2, 0, 2) = (b0 * rtilde11 - b1 * rtilde10) * rv00 + (b1 * rtilde00 - b0 * rtilde01) * rv01
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv02;
    MATD_EL(_R2, 1, 0) = (rtilde00)*rv10 + (rtilde10)*rv11 + (-b0) * rv12;
    MATD_EL(_R2, 1, 1) = (rtilde01)*rv10 + (rtilde11)*rv11 + (-b1) * rv12;
    MATD_EL(_R2, 1, 2) = (b0 * rtilde11 - b1 * rtilde10) * rv10 + (b1 * rtilde00 - b0 * rtilde01) * rv11
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv12;
    MATD_EL(_R2, 2, 0) = (rtilde00)*rv20 + (rtilde10)*rv21 + (-b0) * rv22;
    MATD_EL(_R2, 2, 1) = (rtilde01)*rv20 + (rtilde11)*rv21 + (-b1) * rv22;
    MATD_EL(_R2, 2, 2) = (b0 * rtilde11 - b1 * rtilde10) * rv20 + (b1 * rtilde00 - b0 * rtilde01) * rv21
                          + (rtilde00 * rtilde11 - rtilde01 * rtilde10) * rv22;
}

void ippeComputeTranslation(matd_t *_objectPoints, matd_t *_imgPoints, matd_t *_R, matd_t *_t)
{
    /**
     * This is solved by building the linear system At = b, where t corresponds to the (ALL_DICTS) translation.
     * This is then inverted with the associated normal equations to give t = inv(transpose(A)*A)*transpose(A)*b
     * For efficiency we only store the coefficients of (transpose(A)*A) and (transpose(A)*b)
     */
    assert(_objectPoints != nullptr);
    assert(_imgPoints != nullptr);
    assert(_R != nullptr);
    assert(_R->nrows == 3 && _R->ncols == 3);
    assert(_objectPoints->ncols == 3);
    assert(_imgPoints->ncols == 2);
    assert(_objectPoints->nrows == _imgPoints->nrows);
    
    int numPts = _objectPoints->nrows;
    
    if (_t != nullptr) matd_destroy(_t);
    _t = matd_create(3, 1);
    
    /* coefficients of (transpose(A)*A) */
    double ATA00 = numPts;
    double ATA02 = 0;
    double ATA11 = numPts;
    double ATA12 = 0;
    double ATA20 = 0;
    double ATA21 = 0;
    double ATA22 = 0;

    /* coefficients of (transpose(A)*b) */
    double ATb0 = 0;
    double ATb1 = 0;
    double ATb2 = 0;

    /* S  gives inv(transpose(A)*A)/det(A)^2 */
    double S00, S01, S02;
    double S10, S11, S12;
    double S20, S21, S22;

    double rx, ry, rz;
    double a2;
    double b2;
    double bx, by;
    
    /* now loop through each point and increment the coefficients */
    for (int i = 0; i < numPts; i++) {
        /* rotation(3x3) * vector(3 * 1) */
        rx = MATD_EL(_R, 0, 0) * MATD_EL(_objectPoints, i, 0) + MATD_EL(_R, 0, 1) * MATD_EL(_objectPoints, i, 1)
            + MATD_EL(_R, 0, 2) * MATD_EL(_objectPoints, i, 2);
        ry = MATD_EL(_R, 1, 0) * MATD_EL(_objectPoints, i, 0) + MATD_EL(_R, 1, 1) * MATD_EL(_objectPoints, i, 1)
            + MATD_EL(_R, 1, 2) * MATD_EL(_objectPoints, i, 2);
        rz = MATD_EL(_R, 2, 0) * MATD_EL(_objectPoints, i, 0) + MATD_EL(_R, 2, 1) * MATD_EL(_objectPoints, i, 1)
            + MATD_EL(_R, 2, 2) * MATD_EL(_objectPoints, i, 2);
            
        a2 = -MATD_EL(_imgPoints, i, 0);
        b2 = -MATD_EL(_imgPoints, i, 1);
        
        ATA02 = ATA02 + a2;
        ATA12 = ATA12 + b2;
        ATA20 = ATA20 + a2;
        ATA21 = ATA21 + b2;
        ATA22 = ATA22 + a2 * a2 + b2 * b2;
        
        bx = MATD_EL(_imgPoints, i, 0) * rz - rx;
        by = MATD_EL(_imgPoints, i, 1) * rz - ry;
        
        ATb0 = ATb0 + bx;
        ATb1 = ATb1 + by;
        ATb2 = ATb2 + a2 * bx + b2 * by;
    }
    
    double detAInv = 1.0 / (ATA00 * ATA11 * ATA22 - ATA00 * ATA12 * ATA21 - ATA02 * ATA11 * ATA20);

    /* construct S */
    S00 =  ATA11 * ATA22 - ATA12 * ATA21;
    S01 =  ATA02 * ATA21;
    S02 = -ATA02 * ATA11;
    S10 =  ATA12 * ATA20;
    S11 =  ATA00 * ATA22 - ATA02 * ATA20;
    S12 = -ATA00 * ATA12;
    S20 = -ATA11 * ATA20;
    S21 = -ATA00 * ATA21;
    S22 =  ATA00 * ATA11;
    
    /* solve t */
    MATD_EL(_t, 0, 0) = detAInv * (S00 * ATb0 + S01 * ATb1 + S02 * ATb2);
    MATD_EL(_t, 1, 0) = detAInv * (S10 * ATb0 + S11 * ATb1 + S12 * ATb2);
    MATD_EL(_t, 2, 0) = detAInv * (S20 * ATb0 + S21 * ATb1 + S22 * ATb2);
}

float ippeEvalReprojError(matd_t *_R, matd_t *_t, matd_t *_objectPoints, matd_t *_undistortedPoints)
{
    assert(_R != nullptr);
    assert(_t != nullptr);
    assert(_R->nrows == 3 && _R->ncols == 3);
    assert(_t->nrows == 3 && _t->ncols == 1);
    assert(_objectPoints->ncols == 3 && _undistortedPoints->ncols == 2);
    assert(_objectPoints->nrows == _undistortedPoints->nrows);
    
    int numPts = _objectPoints->nrows;
    float px, py, pz;
    float reprojError = 0;
    float dx, dy;   // residual reprojection error with respect to x and y coordinates
    
    /* now loop over each correspondence and compute the reprojection error */
    for (int i = 0; i < numPts; i++) {
        /* TransformationMatrix[R t](4x4) * vector[x, y, z, 1](4x1) */
        px = static_cast<float>(MATD_EL(_R, 0, 0) * MATD_EL(_objectPoints, i, 0)) +
             static_cast<float>(MATD_EL(_R, 0, 1) * MATD_EL(_objectPoints, i, 1)) +
             static_cast<float>(MATD_EL(_R, 0, 2) * MATD_EL(_objectPoints, i, 2)  + MATD_EL(_t, 0, 0));
        py = static_cast<float>(MATD_EL(_R, 1, 0) * MATD_EL(_objectPoints, i, 0)) +
             static_cast<float>(MATD_EL(_R, 1, 1) * MATD_EL(_objectPoints, i, 1)) +
             static_cast<float>(MATD_EL(_R, 1, 2) * MATD_EL(_objectPoints, i, 2)  + MATD_EL(_t, 1, 0));
        pz = static_cast<float>(MATD_EL(_R, 2, 0) * MATD_EL(_objectPoints, i, 0)) +
             static_cast<float>(MATD_EL(_R, 2, 1) * MATD_EL(_objectPoints, i, 1)) +
             static_cast<float>(MATD_EL(_R, 2, 2) * MATD_EL(_objectPoints, i, 2)  + MATD_EL(_t, 2, 0));
        
        dx = px / pz - MATD_EL(_undistortedPoints, i, 0);
        dy = py / pz - MATD_EL(_undistortedPoints, i, 1);
        
        reprojError += sqrt(dx*dx + dy*dy);
    }
    
    return reprojError;
}

void undistortPoints(matd_t *_imagePoints, matd_t *_undistortedPoints, matd_t *_cameraMatrix, matd_t *_distCoeffs)
{
    
}

} /* namespace ippe */









