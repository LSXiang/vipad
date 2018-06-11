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

#include <iostream>
#include <cmath>
#include <assert.h>

#include "lpe.h"

namespace vipad {

    /**
 * FUNCTION: _quaternion2rotationMatrix
 * 
 * conversion quaternion to rotation matrix
 */
void LocalPositionEstimation::_quaternion2rotationMatrix(const Quaternion &q, apriltag::matd_t* &R)
{
    if (R) apriltag::matd_destroy(R);
    R = apriltag::matd_create(3, 3);
    
    double sqw = q.w * q.w;
    double sqx = q.x * q.x;
    double sqy = q.y * q.y;
    double sqz = q.z * q.z;
    
    /* invs (inverse square length) is only required if quaternion is not already normalised */
    double invs = 1 / (sqx + sqy + sqz + sqw);  /* since sqw + sqx + sqy + sqz = 1/invs*invs */
    MATD_EL(R, 0, 0) = ( sqx - sqy - sqz + sqw) * invs;
    MATD_EL(R, 1, 1) = (-sqx + sqy - sqz + sqw) * invs;
    MATD_EL(R, 2, 2) = (-sqx - sqy + sqz + sqw) * invs;
    
    double tmp1 = q.x * q.y;
    double tmp2 = q.z * q.w;
    MATD_EL(R, 1, 0) = 2.0 * (tmp1 + tmp2) * invs;
    MATD_EL(R, 0, 1) = 2.0 * (tmp1 - tmp2) * invs;
    
    tmp1 = q.x * q.z;
    tmp2 = q.y * q.w;
    MATD_EL(R, 2, 0) = 2.0 * (tmp1 - tmp2) * invs;
    MATD_EL(R, 0, 2) = 2.0 * (tmp1 + tmp2) * invs;
    
    tmp1 = q.y * q.z;
    tmp2 = q.x * q.w;
    MATD_EL(R, 2, 1) = 2.0 * (tmp1 + tmp2) * invs;
    MATD_EL(R, 1, 2) = 2.0 * (tmp1 - tmp2) * invs;
}

/**
 * FUNCTION: _rotationMatrix2quaternion
 * 
 * conversion rotation matrix to quaternion
 */
void LocalPositionEstimation::_rotationMatrix2quaternion(const apriltag::matd_t *R, Quaternion &q)
{
    float tr = MATD_EL(R, 0, 0) + MATD_EL(R, 1, 1) + MATD_EL(R, 2, 2);
    
    if (tr > 0) {
        float S = sqrt(tr + 1.0) * 2;     /* S = 4 * q.w */
        q.w = 0.25 * S;
        q.x = (MATD_EL(R, 2, 1) - MATD_EL(R, 1, 2)) / S;
        q.y = (MATD_EL(R, 0, 2) - MATD_EL(R, 2, 0)) / S;
        q.z = (MATD_EL(R, 1, 0) - MATD_EL(R, 0,1 )) / S;
        
    } else if ((MATD_EL(R, 0, 0) > MATD_EL(R, 1, 1)) & (MATD_EL(R, 0, 0) > MATD_EL(R, 2, 2))) {
        float S = sqrt(1.0 + MATD_EL(R, 0, 0) - MATD_EL(R, 1, 1) - MATD_EL(R, 2, 2)) * 2;    /* S = 4 * q.x */
        q.w = (MATD_EL(R, 2, 1) - MATD_EL(R, 1, 2)) / S;
        q.x = 0.25 * S;
        q.y = (MATD_EL(R, 0, 1) + MATD_EL(R, 1, 0)) / S;
        q.z = (MATD_EL(R, 0, 2) + MATD_EL(R, 2, 0)) / S;
        
    } else if (MATD_EL(R, 1, 1) > MATD_EL(R, 2, 2)) {
        float S = sqrt(1.0 + MATD_EL(R, 1, 1) - MATD_EL(R, 0, 0) - MATD_EL(R, 2, 2)) * 2;    /* S = 4 * q.y */
        q.w = (MATD_EL(R, 0, 2) - MATD_EL(R, 2, 0)) / S;
        q.x = (MATD_EL(R, 0, 1) + MATD_EL(R, 1, 0)) / S;
        q.y = 0.25 * S;
        q.z = (MATD_EL(R, 1, 2) + MATD_EL(R, 2, 1)) / S;
        
    } else {
        float S = sqrt(1.0 + MATD_EL(R, 2, 2) - MATD_EL(R, 0, 0) - MATD_EL(R, 1, 1)) * 2;    /* S = 4 * q.z */
        q.w = (MATD_EL(R, 1, 0) - MATD_EL(R, 0, 1)) / S;
        q.x = (MATD_EL(R, 0, 2) + MATD_EL(R, 2, 0)) / S;
        q.y = (MATD_EL(R, 1, 2) + MATD_EL(R, 2, 1)) / S;
        q.z = 0.25 * S;
    }
}

/**
 * FUNCTION: _get_euler_roll
 * 
 * get euler roll angle from quaternion
 */
float LocalPositionEstimation::_get_euler_roll(const Quaternion &q) const
{
    return (std::atan2<float>(2.0f * (q.w * q.x + q.y * q.z), 1.0f - 2.0 * (q.x * q.x + q.y * q.y)));
}

/**
 * FUNCTION: _get_euler_roll
 * 
 * get euler roll angle from rotation
 */
float LocalPositionEstimation::_get_euler_roll(const apriltag::matd_t *R) const
{
    assert((R != NULL) && (R->ncols == 3) && (R->nrows == 3));
    
    float pitch = std::asin(-MATD_EL(R, 2, 0));
    
    if (std::abs(pitch - M_PI_2) < 1.0e-3f) {
        return 0.0f;
        
    } else if (std::abs(pitch + M_PI_2) < 1.0e-3f) {
        return 0.0f;
        
    } else {
        return std::atan2<float>(MATD_EL(R, 2, 1), MATD_EL(R, 2, 2));
    }
}

/**
 * FUNCTION: _get_euler_pitch
 * 
 * get euler pitch angle from quaternion
 */
float LocalPositionEstimation::_get_euler_pitch(const Quaternion &q) const
{
    return (std::asin(2.0f * (q.w * q.y - q.z * q.x)));
}

/**
 * FUNCTION: _get_euler_pitch
 * 
 * get euler roll angle from rotation
 */
float LocalPositionEstimation::_get_euler_pitch(const apriltag::matd_t *R) const
{
    assert((R != NULL) && (R->ncols == 3) && (R->nrows == 3));
    
    return (std::asin(-MATD_EL(R, 2, 0)));
}

/**
 * FUNCTION: _get_euler_pitch
 * 
 * get euler pitch angle from quaternion
 */
float LocalPositionEstimation::_get_euler_yaw(const Quaternion &q) const
{
    return (std::atan2<float>(2.0f * (q.w * q.z + q.x * q.y), 1.0f - 2.0f * (q.y * q.y + q.z * q.z)));
}

/**
 * FUNCTION: _get_euler_yaw
 * 
 * get euler yaw angle from rotation
 */
float LocalPositionEstimation::_get_euler_yaw(const apriltag::matd_t *R) const
{
    assert((R != NULL) && (R->ncols == 3) && (R->nrows == 3));
    
    float pitch = std::asin(-MATD_EL(R, 2, 0));
    
    if (std::abs(pitch - M_PI_2) < 1.0e-3f) {
        return (std::atan2<float>(MATD_EL(R, 1, 2) - MATD_EL(R, 0, 1), MATD_EL(R, 0, 2) + MATD_EL(R, 1, 1)) + pitch);
        
    } else if (std::abs(pitch + M_PI_2) < 1.0e-3f) {
        return (std::atan2<float>(MATD_EL(R, 1, 2) - MATD_EL(R, 0, 1), MATD_EL(R, 0, 2) + MATD_EL(R, 1, 1)) - pitch);
        
    } else {
        return std::atan2<float>(MATD_EL(R, 1, 0), MATD_EL(R, 0, 0));
    }
}

/**
 * FUNCTION: _eulerAngel2rotationMatrix
 * 
 * conversion euler angel to rotation matrix
 */
void LocalPositionEstimation::_eulerAngel2rotationMatrix(const float roll, const float pitch, const float yaw, apriltag::matd_t* &R)
{
    if (R) apriltag::matd_destroy(R);
    R = apriltag::matd_create(3, 3);
    
    float cp = std::cos(pitch);
    float sp = std::sin(pitch);
    float sr = std::sin(roll);
    float cr = std::cos(roll);
    float sy = std::sin(yaw);
    float cy = std::cos(yaw);
    
    MATD_EL(R, 0, 0) = cp * cy;
    MATD_EL(R, 0, 1) = (sr * sp * cy) - (cr * sy);
    MATD_EL(R, 0, 2) = (cr * sp * cy) + (sr * sy);
    MATD_EL(R, 1, 0) = cp * sy;
    MATD_EL(R, 1, 1) = (sr * sp * sy) + (cr * cy);
    MATD_EL(R, 1, 2) = (cr * sp * sy) - (sr * cy);
    MATD_EL(R, 2, 0) = -sp;
    MATD_EL(R, 2, 1) = sr * cp;
    MATD_EL(R, 2, 2) = cr * cp;
}

} /* namespace vipad */

