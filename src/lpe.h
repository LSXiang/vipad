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

/**
 * This Local Position Estimation (LPE) is a way that used visual marker and IMU data 
 * to accurately compute real-time local position for camera. 
 */

#ifndef __LOCAT_POSITION_ESTIMATION__
#define __LOCAT_POSITION_ESTIMATION__

#include "apriltag.h"
#include "tag_family.h"
#include "ippe.h"

namespace vipad {

struct localPosition
{
    float yaw;  /* unit: rad */
    float x;    /* unit: meter */
    float y;    /* unit: meter */
    float z;    /* unit: meter*/
    int id;     /* using the this ID number to local position estimation now */
    uint32_t tags_num; /* in this frame, have n number be detect */
};

struct Quaternion
{
    float w;    /* Quaternion component 1, w (1 in null-rotation) */
    float x;    /* Quaternion component 2, x (0 in null-rotation) */
    float y;    /* Quaternion component 3, y (0 in null-rotation) */
    float z;    /* Quaternion component 4, z (0 in null-rotation) */
};

/**
 * The angle between the camera and the flight controller head
 * (clockwise rotation of the camera to the flight control head coincident angle)
 */
enum cam2droneRotation
{
	ClockwiseAngle_0,
	ClockwiseAngle_90,
	ClockwiseAngle_180,
	ClockwiseAngle_270,
	ClockwiseAngle360
};

struct lpe_params
{
	uint8_t *input;
	uint32_t width;
	uint32_t height;

	double fx, fy, cx, cy;
	double k1, k2, p1, p2, k3;

	cam2droneRotation angle;
	Quaternion *q;

	localPosition *locat;

	enum apriltag::tag_family_type tag_type;
	float marker_length;
};

class LocalPositionEstimation
{
public:
    /* Constructor */
    LocalPositionEstimation(void *param);
    /* Destructor */
    ~LocalPositionEstimation(void);

    /* estimate local position */
    void estimateLocalPosition(void);
    
private:
    /**
     * @brief conversion quaternion to rotation matrix 
     * 
     * @param q quaternion (in)
     * @param R rotation matrix (out)
     */
    void _quaternion2rotationMatrix(const Quaternion &q, apriltag::matd_t* &R);
    
    /**
     * @brief conversion rotation matrix to quaternion
     * 
     * @param R rotation matrix (in)
     * @param q quaternion (out)
     */
    void _rotationMatrix2quaternion(const apriltag::matd_t *R, Quaternion &q);
    
    /**
     * @brief get euler roll angle
     * 
     * @param q quaternion
     * @return euler roll angle
     */
    float _get_euler_roll(const Quaternion &q) const;
    /**
     * @brief get euler roll angle
     * 
     * @param R rotation matrix
     * @return euler roll angle
     */
    float _get_euler_roll(const apriltag::matd_t *R) const;
    
    /**
     * @brief get euler pitch angle
     * 
     * @param q quaternion
     * @return euler pitch angle
     */
    float _get_euler_pitch(const Quaternion &q) const;
    /**
     * @brief get euler pitch angle
     * 
     * @param R rotation matrix
     * @return euler pitch angle
     */
    float _get_euler_pitch(const apriltag::matd_t *R) const;
    
    
    /**
     * @brief get euler yaw angle
     * 
     * @param q quaternion
     * @return euler yaw angle
     */
    float _get_euler_yaw(const Quaternion &q) const;
    /**
     * @brief get euler yaw angle
     * 
     * @param R rotation matrix
     * @return euler yaw angle
     */
    float _get_euler_yaw(const apriltag::matd_t *R) const;
    
    /**
     * @brief conversion euler angel to rotation matrix
     * 
     * @param roll euler roll angle
     * @param pitch euler pitch angle
     * @param yaw euler yaw angle
     * @param R rotation matrix
     */
    void _eulerAngel2rotationMatrix(const float roll, const float pitch, const float yaw, apriltag::matd_t* &R);
    
    apriltag::apriltag_family_t *_marker_family;
    apriltag::apriltag_detector_t *_marker_detector;

    apriltag::matd_t *_camera_matrix, *_dist_param;

    struct lpe_params *_param;
};


} /* namespace vipad */

#endif /* __LOCAT_POSITION_ESTIMATION__ */

