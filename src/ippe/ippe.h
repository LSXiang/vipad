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

#ifndef __IPPE_H__
#define __IPPE_H__

#include "matd.h"

namespace ippe {

/** 
 * @brief Finds the two possible poses of a square planar object given its four corner correspondences in an image
 * using IPPE. These poses are sorted so that the first one is the one with the lowest reprojection error. The second
 * pose is needed if the problem is ambiguous. The problem is ambiguous when the projection of the model is close to
 * affine, which in practice happens if it is small or viewed from a large distance. In these cases there are two pose
 * solutions that can correctly align the correspondences (up to noise), so it is impossible to select the right one
 * from just the reprojection error. IPPE gives you both the solutions, rather than just a single solution (which in
 * ambiguous cases would be wrong 50% of the time). Geometrically, the two poses roughly correspond to a reflection of
 * the model about a plane whose normal passes through the line-of-sight from the camera centre to the model's centre.
 * For more details about these ambiguities, please refer to the IPPE paper.
 *
 * It is possible to reject the second pose if its reprojection error is significantly worse than the first pose.
 * The correct way to do this is with a likelihood ratio test (https://en.wikipedia.org/wiki/Likelihood-ratio_test).
 *
 * @param squareLength the squares's length (which is also it's width) in object coordinate units (e.g. millimeters,
 *                     meters, etc.)
 * The square is defined in object coordinates on the plane z=0 and centred at the origin. Therefore its four points in
 * object coordinates are given by:
 *   point 0: [-squareLength / 2.0,  squareLength / 2.0, 0]
 *   point 1: [ squareLength / 2.0,  squareLength / 2.0, 0]
 *   point 2: [ squareLength / 2.0, -squareLength / 2.0, 0]
 *   point 3: [-squareLength / 2.0, -squareLength / 2.0, 0]
 * @param imagePoints Array of four corresponding image points, 1x4/4x1 2-channel. Note that the points should be
 *                    ordered to correspond with points 0, 1, 2 and 3.
 * @param cameraMatrix Input camera matrix \f$A = \vecthreethree{fx}{0}{cx}{0}{fy}{cy}{0}{0}{1}\f$ .
 * @param distCoeffs Input vector of distortion coefficients (k_1, k_2, p_1, p_2) of 4 elements. If the vector is
 *                   NULL/empty, the zero distortion coefficients are assumed.
 * @param _rvec1 Output rotation vector (see Rodrigues ) for first pose. That, together with tvec , brings points from
 *               the model coordinate system to the camera coordinate system.
 * @param _tvec1 Output translation vector for first pose.
 * @param reprojErr1 Output reprojection error of first pose
 * @param _rvec2 Output rotation vector (see Rodrigues ) for second pose. That, together with tvec , brings points from
 *               the model coordinate system to the camera coordinate system.
 * @param _tvec2 Output translation vector for second pose.
 * @param reprojErr1 Output reprojection error of second pose
 */
void solvePoseOfCentredSquare(float squareLength, matd_t imagePoints, matd_t *cameraMatrix,
                                matd_t *distCoeffs, matd_t *_rvec1, matd_t *_tvec1,
                                float& reprojErr1, matd_t *_rvec2, matd_t *_tvec2, float& reprojErr2);

/** 
 * @brief Computes the translation solution for a given rotation solution
 * @param _objectPoints Array of corresponding model points, 1xN/Nx1 3-channel where N is the number of points
 * @param _undistortedPoints Array of corresponding image points (undistorted), 1xN/Nx1 2-channel where N is the number of points
 * @param _R1 Rotation solution from IPPE, 3x3 1-channel float
 * @param _t  Translation solution, 3x1 1-channel float
 */
void ippeComputeTranslation(matd_t *_objectPoints, matd_t *_imgPoints, matd_t *_R, matd_t *_t);

/** 
 * @brief Computes the two rotation solutions from the Jacobian of a homography matrix H. For highest accuracy the
 *        Jacobian should be computed at the centroid of the point correspondences (see the IPPE paper for the 
 *        explaination of this). For a point (ux,uy) on the model plane, suppose the homography H maps (ux,uy) to a point (p,q)
 *        in the image (in normalised pixel coordinates). The Jacobian matrix [J00, J01; J10,J11] is the Jacobian
 *        of the mapping evaluated at (ux,uy).
 * @param j00 Jacobian coefficent
 * @param j01 Jacobian coefficent
 * @param j10 Jacobian coefficent
 * @param j11 Jacobian coefficent
 * @param p the x coordinate of point (ux,uy) mapped into the image (undistorted and normalised position)
 * @param q the y coordinate of point (ux,uy) mapped into the image (undistorted and normalised position)
 */
void ippeComputeRotations(double j00, double j01, double j10, double j11, double p, double q, matd_t *_R1, matd_t *_R2);

/**
 * @brief Closed-form solution for the homography mapping with four corner correspondences of a square (it maps
 * source points to target points). The source points are the four corners of a zero-centred squared defined by:
 *  point 0: [-squareLength / 2.0,  squareLength / 2.0]
 *  point 1: [ squareLength / 2.0,  squareLength / 2.0]
 *  point 2: [ squareLength / 2.0, -squareLength / 2.0]
 *  point 3: [-squareLength / 2.0, -squareLength / 2.0]
 *
 * @param _targetPts Array of four corresponding target points, 4x2. Note that the points should be
 *                   ordered to correspond with points 0, 1, 2 and 3.
 * @param halfLength the square's half length (i.e. squareLength/2.0)
 * @param _R1 Rotation solution from IPPE, 3x3 1-channel float
 * @param _H  Homograhy mapping the source points to the target points, 3x3 single channel
 */
void homographyFromSquarePoints(matd_t *_targetPts, float halfLength, matd_t *_H);

/**
 * @brief Finds the rotation _Ra that rotates a vector _a to the z axis (0, 0, 1)
 * @param _a  vector: 3x1 matd_t (double)
 * @param _Ra Rotation: 3x3 matd_t (double)
 */
// void rotateVec2ZAxis(matd_t *_a, matd_t *_Ra);

} /* namespace ippe */

#endif /* __IPPE_H__ */


