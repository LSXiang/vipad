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

#include "ippe.h"

namespace ippe {



void homographyFromSquarePoints(matd_t *_targetPts, float halfLength, matd_t *_H)
{
    assert(_targetPts->nrows == 4 && _targetPts->ncols == 2);
    
    if (_H != nullptr)
        matd_destroy(_H);
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

} /* namespace ippe */

