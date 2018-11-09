/**
 * Copyright (c) 2018, The Akatsuki(jacob.lsx). All rights reserved.
 * Adaptation of tobycollins/IPPE.
 * https://github.com/tobycollins/IPPE
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

void solvePoseOfMarker(float markerLength, matd_t* &imagePoints, matd_t* &cameraMatrix, matd_t* &distCoeffs,
                       matd_t* &R, matd_t* &t, float &reprojErr)
{
    matd_t *_R = NULL, *_t = NULL;
    float _reprojErr;
    
    ippeSolvePoseOfCentredSquare(markerLength, imagePoints, cameraMatrix, distCoeffs,
                                 R, t, reprojErr, _R, _t, _reprojErr);
    
    matd_destroy(_R);
    matd_destroy(_t);
}

void ippeSolvePoseOfCentredSquare(float squareLength, matd_t* &imagePoints, matd_t* &cameraMatrix, matd_t* &distCoeffs,
                                  matd_t* &_R1, matd_t* &_tvec1, float& reprojErr1, 
                                  matd_t* &_R2, matd_t* &_tvec2, float& reprojErr2)
{
    assert(imagePoints != NULL);
    assert(cameraMatrix != NULL);
    assert(imagePoints->nrows == 4 && imagePoints->ncols == 2);
    assert(cameraMatrix->nrows == 3 && cameraMatrix->ncols == 3);
    assert(distCoeffs == NULL || distCoeffs->nrows * distCoeffs->ncols == 5);
    
    matd_t *undistortedPoints = NULL;    // undistorted version of imagePoints
    matd_t *modelPoints = matd_create(4, 3);
    
    if (_R1 != NULL) matd_destroy(_R1);
    if (_tvec1 != NULL) matd_destroy(_tvec1);
    if (_R2 != NULL) matd_destroy(_R2);
    if (_tvec2 != NULL) matd_destroy(_tvec2);
    
    /* set coordinate system in the middle of the centred square, with Z pointing out */
    float halfLength = squareLength / 2.0f;
    MATD_EL(modelPoints, 0, 0) = -halfLength;   MATD_EL(modelPoints, 0, 1) =  halfLength;   MATD_EL(modelPoints, 0, 2) = 0;
    MATD_EL(modelPoints, 1, 0) =  halfLength;   MATD_EL(modelPoints, 1, 1) =  halfLength;   MATD_EL(modelPoints, 1, 2) = 0;
    MATD_EL(modelPoints, 2, 0) =  halfLength;   MATD_EL(modelPoints, 2, 1) = -halfLength;   MATD_EL(modelPoints, 2, 2) = 0;
    MATD_EL(modelPoints, 3, 0) = -halfLength;   MATD_EL(modelPoints, 3, 1) = -halfLength;   MATD_EL(modelPoints, 3, 2) = 0;
    
    /* (Ra, ta), (Rb, tb) are the two pose solutions from IPPE */
    matd_t *Ra = NULL, *ta = NULL;
    matd_t *Rb = NULL, *tb = NULL;
    matd_t *H = NULL;  // homography matrix pointer
    
    /* undistort the image points (i.e. put them in normalized pixel coordinates) */
    undistortPoints(imagePoints, undistortedPoints, cameraMatrix, distCoeffs);
    
    /* compute the homography mapping the model's four corners to undistortedPoints */
    homographyFromSquarePoints(undistortedPoints, halfLength, H);
    
    /* compute the Jacobian J of the homography at (0,0) */
    float j00, j01, j10, j11, v0, v1;
    
    j00 = MATD_EL(H, 0, 0) - MATD_EL(H, 2, 0) * MATD_EL(H, 0, 2);
    j01 = MATD_EL(H, 0, 1) - MATD_EL(H, 2, 1) * MATD_EL(H, 0, 2);
    j10 = MATD_EL(H, 1, 0) - MATD_EL(H, 2, 0) * MATD_EL(H, 1, 2);
    j11 = MATD_EL(H, 1, 1) - MATD_EL(H, 2, 1) * MATD_EL(H, 1, 2);
    
    /* compute the transformation of (0,0) into the image */
    v0 = MATD_EL(H, 0, 2);
    v1 = MATD_EL(H, 1, 2);
    
    /* compute the two rotation solutions */
    ippeComputeRotations(j00, j01, j10, j11, v0, v1, Ra, Rb);
    
    /* for each rotation solution, compute the corresponding translation solution */
    ippeComputeTranslation(modelPoints, undistortedPoints, Ra, ta);
    ippeComputeTranslation(modelPoints, undistortedPoints, Rb, tb);
    
    /* for each transformation (R, t) solution, compute the reprojection error */
    float reprojErra = ippeEvalReprojError(Ra, ta, modelPoints, undistortedPoints);
    float reprojErrb = ippeEvalReprojError(Rb, tb, modelPoints, undistortedPoints);
    
    if (reprojErra < reprojErrb) {
        _R1 = matd_copy(Ra);
        _tvec1 = matd_copy(ta);
        
        _R2 = matd_copy(Rb);
        _tvec2 = matd_copy(tb);
        
        reprojErr1 = reprojErra;
        reprojErr2 = reprojErrb;
    } else {
        _R1 = matd_copy(Rb);
        _tvec1 = matd_copy(tb);
        
        _R2 = matd_copy(Ra);
        _tvec2 = matd_copy(ta);
        
        reprojErr1 = reprojErrb;
        reprojErr2 = reprojErra;
    }
    
    matd_destroy(modelPoints);
    matd_destroy(undistortedPoints);
    matd_destroy(H);
    matd_destroy(Ra);
    matd_destroy(Rb);
    matd_destroy(ta);
    matd_destroy(tb);
}

void homographyFromSquarePoints(matd_t* &_targetPts, float halfLength, matd_t* &_H)
{
    assert(_targetPts != NULL);
    assert(_targetPts->nrows == 4 && _targetPts->ncols == 2);
    
    if (_H != NULL) matd_destroy(_H);
    _H = matd_create(3, 3);
    
    float p1x = -MATD_EL(_targetPts, 0, 0);
    float p1y = -MATD_EL(_targetPts, 0, 1);
    
    float p2x = -MATD_EL(_targetPts, 1, 0);
    float p2y = -MATD_EL(_targetPts, 1, 1);
    
    float p3x = -MATD_EL(_targetPts, 2, 0);
    float p3y = -MATD_EL(_targetPts, 2, 1);
    
    float p4x = -MATD_EL(_targetPts, 3, 0);
    float p4y = -MATD_EL(_targetPts, 3, 1);
    
    /* analytic solution */
    float detsInv = -1 / (halfLength * (p1x * p2y - p2x * p1y - p1x * p4y + p2x * p3y - p3x * p2y + p4x * p1y
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

void ippeComputeRotations(float j00, float j01, float j10, float j11, float p, float q, matd_t* &_R1, matd_t* &_R2)
{
    /**
     * Note that it is very hard to understand what is going on here from the code, so if you want to have a clear explanation
     * then please refer to the IPPE paper (Algorithm 1 and its description). Or, you can read IPPE matlab code that is certainly 
     * easier to read.
     */
    if (_R1 != NULL) matd_destroy(_R1);
    _R1 = matd_create(3, 3);
    
    if (_R2 != NULL) matd_destroy(_R2);
    _R2 = matd_create(3, 3);
    
    float a00, a01, a10, a11, ata00, ata01, ata11, b00, b01, b10, b11, binv00, binv01, binv10, binv11;
    float rtilde00, rtilde01, rtilde10, rtilde11;
    float rtilde00_2, rtilde01_2, rtilde10_2, rtilde11_2;
    float b0, b1, gamma, dtinv;
    float sp;
    
    /* Finds the rotation Rv that rotates a vector (p, q, 1) to the z axis (0, 0, 1) */
    float s, t, krs0, krs1, krs0_2, krs1_2, costh, sinth;
    
    s = sqrt(p * p + q * q + 1);
    t = sqrt(p * p + q * q);
    costh = 1 / s;
    sinth = sqrt(1 - 1 / (s * s));

    krs0 = p / t;
    krs1 = q / t;
    krs0_2 = krs0 * krs0;
    krs1_2 = krs1 * krs1;
    
    float rv00, rv01, rv02;
    float rv10, rv11, rv12;
    float rv20, rv21, rv22;
    
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

void ippeComputeTranslation(matd_t* &_objectPoints, matd_t* &_imgPoints, matd_t* &_R, matd_t* &_t)
{
    /**
     * This is solved by building the linear system At = b, where t corresponds to the (ALL_DICTS) translation.
     * This is then inverted with the associated normal equations to give t = inv(transpose(A)*A)*transpose(A)*b
     * For efficiency we only store the coefficients of (transpose(A)*A) and (transpose(A)*b)
     */
    assert(_objectPoints != NULL);
    assert(_imgPoints != NULL);
    assert(_R != NULL);
    assert(_R->nrows == 3 && _R->ncols == 3);
    assert(_objectPoints->ncols == 3);
    assert(_imgPoints->ncols == 2);
    assert(_objectPoints->nrows == _imgPoints->nrows);
    
    int numPts = _objectPoints->nrows;
    
    if (_t != NULL) matd_destroy(_t);
    _t = matd_create(3, 1);
    
    /* coefficients of (transpose(A)*A) */
    float ATA00 = numPts;
    float ATA02 = 0;
    float ATA11 = numPts;
    float ATA12 = 0;
    float ATA20 = 0;
    float ATA21 = 0;
    float ATA22 = 0;

    /* coefficients of (transpose(A)*b) */
    float ATb0 = 0;
    float ATb1 = 0;
    float ATb2 = 0;

    /* S  gives inv(transpose(A)*A)/det(A)^2 */
    float S00, S01, S02;
    float S10, S11, S12;
    float S20, S21, S22;

    float rx, ry, rz;
    float a2;
    float b2;
    float bx, by;
    
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
    
    float detAInv = 1.0 / (ATA00 * ATA11 * ATA22 - ATA00 * ATA12 * ATA21 - ATA02 * ATA11 * ATA20);

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

float ippeEvalReprojError(matd_t* &_R, matd_t* &_t, matd_t* &_objectPoints, matd_t* &_undistortedPoints)
{
    assert(_R != NULL);
    assert(_t != NULL);
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

void undistortPoints(matd_t* &_imagePoints, matd_t* &_undistortedPoints, matd_t* &_cameraMatrix, matd_t* &_distCoeffs, int _maxITER)
{
    assert(_imagePoints != NULL);
    assert(_cameraMatrix != NULL);
    assert(_imagePoints->ncols == 2);
    assert(_cameraMatrix->nrows == 3 && _cameraMatrix->ncols == 3);
    assert(_distCoeffs == NULL || _distCoeffs->nrows * _distCoeffs->ncols == 5);
    
    float k1, k2, p1, p2, k3;
    
    if (_distCoeffs == NULL) {
        k1 = k2 = p1 = p2 = k3 = 0;
    } else {
        k1 = _distCoeffs->data[0];
        k2 = _distCoeffs->data[1];
        p1 = _distCoeffs->data[2];
        p2 = _distCoeffs->data[3];
        k3 = _distCoeffs->data[4];
    }
    
    float fx = MATD_EL(_cameraMatrix, 0, 0);
    float fy = MATD_EL(_cameraMatrix, 1, 1);
    float ifx = 1./fx;
    float ify = 1./fy;
    float cx = MATD_EL(_cameraMatrix, 0, 2);
    float cy = MATD_EL(_cameraMatrix, 1, 2);
    
    int numPts = _imagePoints->nrows;
    
    if (_undistortedPoints != NULL) matd_destroy(_undistortedPoints);
    _undistortedPoints = matd_create(numPts, 2);
    
    for (int i = 0; i < numPts; i ++) {
        float x, y, x0 = 0, y0 = 0;
        
        x = MATD_EL(_imagePoints, i, 0);
        y = MATD_EL(_imagePoints, i, 1);
        
        x = (x - cx) * ifx;
        y = (y - cy) * ify;
        
        if (_distCoeffs != NULL) {
            // vector (x, y, 1)
            x0 = x;
            y0 = y;
            
            /* compensate distortion iteratively */
            for (int j = 0; j < _maxITER; j ++) {
                float r2 = x*x + y*y;
                float icdist = 1. / (1 + ((k3*r2 + k2)*r2 + k1)*r2);
                float deltaX = 2*p1*x*y + p2*(r2 + 2*x*x);
                float deltaY = p1*(r2 + 2*y*y) + 2*p2*x*y;
                x = (x0 - deltaX)*icdist;
                y = (y0 - deltaY)*icdist;
            }
        }
        
        MATD_EL(_undistortedPoints, i, 0) = x;
        MATD_EL(_undistortedPoints, i, 1) = y;
    }
}

} /* namespace ippe */









