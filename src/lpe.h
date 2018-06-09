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
 * This Local Position Estimatiom (LPE) is a way that used visual marker and IMU data 
 * to accurately compute real-time local position for camera. 
 */

#ifndef __LOCAT_POSITION_ESTIMATION__
#define __LOCAT_POSITION_ESTIMATION__

#include "apriltag.h"
#include "ippe.h"

namespace vipad {

struct localPosition
{
    float yaw;  /* unit: rad */
    float x;    /* unit: meter */
    float y;    /* unit: meter */
    float z;    /* unit: meter*/
};

struct Quaternion
{
    float w;    /* Quaternion component 1, w (1 in null-rotation) */
    float x;    /* Quaternion component 2, x (0 in null-rotation) */
    float y;    /* Quaternion component 3, y (0 in null-rotation) */
    float z;    /* Quaternion component 4, z (0 in null-rotation) */
};

class LocalPositionEstimation
{
public:
    typedef LocalPositionEstimation* ptr;
    
    /* Constructor */
    LocalPositionEstimation();
    /* Destructor */
    ~LocalPositionEstimation();
    
    /* initialize */
    init();
    
    /* estimate local position */
    estimateLocalPosition();
    
private:
    
    apriltag::apriltag_family_t *_marker_family;
    apriltag::apriltag_detector_t *_marker_detector;
};


} /* namespace vipad */

#endif /* __LOCAT_POSITION_ESTIMATION__ */

