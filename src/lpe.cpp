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
#include <limits.h>

#include "lpe.h"

using namespace apriltag;

namespace vipad {
    

/* Constructor */
LocalPositionEstimation::LocalPositionEstimation(void *param) :
	_marker_family(NULL),
	_marker_detector(NULL),
	_camera_matrix(NULL),
	_dist_param(NULL),
	_param(NULL)
{
	_param = (struct lpe_params *)param;

	_marker_family = tags_create(_param->tag_type);
	assert(_marker_family);

	_marker_detector = apriltag_detector_create();
	assert(_marker_detector);
	apriltag_detector_add_family(_marker_detector, _marker_family);

// 	_marker_detector->quad_decimate = 1.0f;
// 	_marker_detector->quad_sigma = 0.8f;
	_marker_detector->refine_edges = 1;

	_camera_matrix = matd_identity(3);
	MATD_EL(_camera_matrix, 0, 0) = _param->fx;
	MATD_EL(_camera_matrix, 0, 2) = _param->cx;
	MATD_EL(_camera_matrix, 1, 1) = _param->fy;
	MATD_EL(_camera_matrix, 1, 2) = _param->cy;

	double dist[] = {_param->k1, _param->k2, _param->p1, _param->p2, _param->k3};
	_dist_param = matd_create_data(1, 5, dist);
}

/* Destructor */
LocalPositionEstimation::~LocalPositionEstimation(void)
{
	matd_destroy(_camera_matrix);
	matd_destroy(_dist_param);
	apriltag_detector_destroy(_marker_detector);
	tags_destroy(_marker_family);
}

/* estimate local position */
void LocalPositionEstimation::estimateLocalPosition(void)
{
	image_u8_t image = { static_cast<int32_t>(_param->width),
						 static_cast<int32_t>(_param->height),
						 static_cast<int32_t>(_param->width),
						 _param->input
					   };

	zarray_t *detections = apriltag_detector_detect(_marker_detector, &image);

	int tags_number = zarray_size(detections);
	_param->locat->tags_num = tags_number;

	printf("%d tags detected \r\n", tags_number);

	/* if at least one marker detected */
	if (tags_number > 0) {

		bool update_usedId = true;
		int index = 0;

		apriltag_detection_t *det;

		if (tags_number > 1) {
			int width_boundary = image.width / 8;
			int height_boundary = image.height / 8;

			/* whether last used_id in the current image */
			for (int i = 0; i < tags_number; i++) {
				zarray_get(detections, i, &det);

				if ((det->id == _param->locat->id) && (det->id < 900)) {
					if ((det->c[0] > width_boundary) &&
						(det->c[0] < image.width - width_boundary) &&
						(det->c[1] > height_boundary) &&
						(det->c[1] < image.height - height_boundary)) {

						update_usedId = false;
						index = i;
					}
					break;
				}
			}

		} else {
			zarray_get(detections, index, &det);
			update_usedId = false;
			_param->locat->id = det->id;
		}

		/* whether update trick marker id */
		if (update_usedId) {
			/* look for the marker at the very center of the camera's field of view */
			unsigned int distance_square = UINT_MAX, mini_distance_square = UINT_MAX;
			for (int i = 0; i < tags_number; i++) {
				zarray_get(detections, i, &det);
//				if (det->id < 900) {
					/* compute distance for center of the marker to center of the  */
					distance_square = std::pow(det->c[0] - image.width/2, 2) + std::pow(det->c[1] - image.height/2, 2);
					if (distance_square < mini_distance_square) {
						mini_distance_square = distance_square;
						index = i;
					}
//				}
			}

			zarray_get(detections, index, &det);
			_param->locat->id = det->id;
		}

//		float x_offset = .0f;
//		float y_offset = .0f;
//		int marker_length_scale = 1;
//
//		if (_param->locat->id >= 900) {
//			/* small marker location offset. Reserve code !!!!!!!!!!!!!! */
//			x_offset = .0f;
//			y_offset = .0f;
//			marker_length_scale = 1;
//		}

		/* estimate location */
		matd_t *image_points;
		image_points = matd_create_data(4, 2, det->p[0]);
        
//         /* ----------------test---------------- */
//         printf("image_points : \r\n %f, %f, %f, %f, %f, %f, %f, %f \r\n", det->p[0][0], det->p[0][1], det->p[1][0], det->p[1][1], det->p[2][0], det->p[2][1], det->p[3][0], det->p[3][1]);
//         std::cout << std::endl;

		matd_t *R = NULL, *t = NULL;
		float reprojErr;

		ippe::solvePoseOfMarker(_param->marker_length, image_points, _camera_matrix, _dist_param, R, t, reprojErr);

		matd_t *R_t = matd_transpose(R);
		_param->locat->yaw = _get_euler_yaw(R_t);
		
		printf("camera yaw angle = %.4f rad (%.2f deg) \r\n", _param->locat->yaw, _param->locat->yaw / M_PI * 180);

		if (_param->q == NULL) {
			matd_t *t_inv = matd_op("-M*M", R_t, t);
			_param->locat->x = MATD_EL(t_inv, 0, 0);
			_param->locat->y = MATD_EL(t_inv, 1, 0);
			_param->locat->z = MATD_EL(t_inv, 2, 0);
			matd_destroy(t_inv);

			printf("camera local position X Y Z = %.4f, %.4f, %.4f \r\n", _param->locat->x, _param->locat->y, _param->locat->z);

		} else {
			/* Construct a new rotation matrix using the pose of the IMU data */
			switch (_param->angle) {
			case ClockwiseAngle_0:
				break;

			case ClockwiseAngle_90:
				break;

			case ClockwiseAngle_180:
				break;

			case ClockwiseAngle_270: {
				float roll  = _get_euler_roll(*_param->q);
				float pitch = _get_euler_pitch(*_param->q);

				matd_t *R_inv;
				_eulerAngel2rotationMatrix(roll < 0 ? roll + M_PI : roll - M_PI, -pitch, _param->locat->yaw, R_inv); // maybe need change!!!!!!!

				matd_t *t_inv = matd_op("-M*M", R_inv, t);
				_param->locat->x = MATD_EL(t_inv, 0, 0);
				_param->locat->y = MATD_EL(t_inv, 1, 0);
				_param->locat->z = MATD_EL(t_inv, 2, 0);

				matd_destroy(R_inv);
				matd_destroy(t_inv);
				}
				break;

			case ClockwiseAngle360:
				break;

			default: {
				printf("Warning: Does not support this rotation angle/ \r\n");
				matd_t *t_inv = matd_op("-M*M", R_t, t);
				_param->locat->x = MATD_EL(t_inv, 0, 0);
				_param->locat->y = MATD_EL(t_inv, 1, 0);
				_param->locat->z = MATD_EL(t_inv, 2, 0);
				matd_destroy(t_inv);
				}
				break;
			}

			printf("camera local position X Y Z = %.4f, %.4f, %.4f \r\n", _param->locat->x, _param->locat->y, _param->locat->z);
		}

		matd_destroy(image_points);
		matd_destroy(R);
        matd_destroy(R_t);
		matd_destroy(t);
	}
	
    zarray_destroy(detections);
    
    std::cout << std::endl;
}

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

