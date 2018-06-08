

#include <iostream>
#include <chrono>

#include <opencv2/opencv.hpp>

#include "apriltag.h"

#include "tag_family.h"

#include "ippe/ippe.h"

using namespace std;
using namespace cv;
using namespace apriltag;

void getRelativeTransform(double tag_size, double fx, double fy, double px, double py, double p[4][2]);
void getRelativeTransform2(double tag_size, double fx, double fy, double cx, double cy, double p[4][2]);

int main(int argc, char *argv[])
{
    // Initialize camera
    VideoCapture cap(0);
    if (!cap.isOpened()) {
        cerr << "Couldn't open video capture device" << endl;
        return -1;
    }

    // Initialize tag detector with options
    apriltag_family_t *tf = NULL;
    
    tf = tags_create(tag36h11);

    apriltag_detector_t *td = apriltag_detector_create();
    apriltag_detector_add_family(td, tf);
    
    td->quad_decimate = 1.0f;
    td->quad_sigma = 0.8f;
    td->refine_edges = 1;
    
    

    Mat frame, gray;
    while (true) {
        cap >> frame;
        cvtColor(frame, gray, COLOR_BGR2GRAY);

        // Make an image_u8_t header for the Mat data
        image_u8_t im = { .width = gray.cols,
            .height = gray.rows,
            .stride = gray.cols,
            .buf = gray.data
        };

        zarray_t *detections = apriltag_detector_detect(td, &im);
        cout << zarray_size(detections) << " tags detected" << endl;

        // Draw detection outlines
        for (int i = 0; i < zarray_size(detections); i++) {
            apriltag_detection_t *det;
            zarray_get(detections, i, &det);
            line(frame, Point(det->p[0][0], det->p[0][1]),
                     Point(det->p[1][0], det->p[1][1]),
                     Scalar(0, 0xff, 0), 2);
            line(frame, Point(det->p[0][0], det->p[0][1]),
                     Point(det->p[3][0], det->p[3][1]),
                     Scalar(0, 0, 0xff), 2);
            line(frame, Point(det->p[1][0], det->p[1][1]),
                     Point(det->p[2][0], det->p[2][1]),
                     Scalar(0xff, 0, 0), 2);
            line(frame, Point(det->p[2][0], det->p[2][1]),
                     Point(det->p[3][0], det->p[3][1]),
                     Scalar(0xff, 0, 0), 2);

            stringstream ss;
            ss << det->id;
            String text = ss.str();
            int fontface = FONT_HERSHEY_SCRIPT_SIMPLEX;
            double fontscale = 1.0;
            int baseline;
            Size textsize = getTextSize(text, fontface, fontscale, 2,
                                            &baseline);
            putText(frame, text, Point(det->c[0]-textsize.width/2,
                                       det->c[1]+textsize.height/2),
                    fontface, fontscale, Scalar(0xff, 0x99, 0), 2);
        }
        
        /* -------------------------------------------------------------------- */
        if (zarray_size(detections) > 0) {
            apriltag_detection_t *det;
            zarray_get(detections, 0, &det);
            
//             matd_t *pose = homography_to_pose(det->H, 2.4908279215754123e+03, 2.4935314568583112e+03, 3.4745731382095448e+02, 2.4094331871742105e+02);
//             
// //             cout << "x y z : " << MATD_EL(pose, 0, 3) * 0.05 <<  ",\t" << MATD_EL(pose, 1, 3) * 0.05 << ",\t" << MATD_EL(pose, 2, 3) * 0.05 << endl;
// //             cout << "R : " << MATD_EL(pose, 0, 0) <<  ",\t" << MATD_EL(pose, 0, 1) << ",\t" << MATD_EL(pose, 0, 2) << "\r\n"
// //                            << MATD_EL(pose, 1, 0) <<  ",\t" << MATD_EL(pose, 1, 1) << ",\t" << MATD_EL(pose, 1, 2) << "\r\n"
// //                            << MATD_EL(pose, 2, 0) <<  ",\t" << MATD_EL(pose, 2, 1) << ",\t" << MATD_EL(pose, 2, 2) << endl;
// //             cout << "yaw :" << std::atan2<float>(MATD_EL(pose, 1, 0), MATD_EL(pose, 0, 0)) * 57.3 << endl;
//             
//             matd_t *R = matd_select(pose, 0, 2, 0, 2);
// //             cout << "R **** : " << MATD_EL(R, 0, 0) <<  ",\t" << MATD_EL(R, 0, 1) << ",\t" << MATD_EL(R, 0, 2) << "\r\n"
// //                         << MATD_EL(R, 1, 0) <<  ",\t" << MATD_EL(R, 1, 1) << ",\t" << MATD_EL(R, 1, 2) << "\r\n"
// //                         << MATD_EL(R, 2, 0) <<  ",\t" << MATD_EL(R, 2, 1) << ",\t" << MATD_EL(R, 2, 2) << endl;
//             
//             matd_t *t_ = matd_select(pose, 0, 2, 3, 3);
// //             cout << "x y z $$$$$$$ : " << MATD_EL(t_, 0, 0) * 0.05 <<  ",\t" << MATD_EL(t_, 1, 0) * 0.05 << ",\t" << MATD_EL(t_, 2, 0) * 0.05 << endl;
//             matd_t *R_t = matd_transpose(R);
// //             cout << "R ^^^^ : " << MATD_EL(R_t, 0, 0) <<  ",\t" << MATD_EL(R_t, 0, 1) << ",\t" << MATD_EL(R_t, 0, 2) << "\r\n"
// //                                 << MATD_EL(R_t, 1, 0) <<  ",\t" << MATD_EL(R_t, 1, 1) << ",\t" << MATD_EL(R_t, 1, 2) << "\r\n"
// //                                 << MATD_EL(R_t, 2, 0) <<  ",\t" << MATD_EL(R_t, 2, 1) << ",\t" << MATD_EL(R_t, 2, 2) << endl;
//             matd_t *t = matd_multiply(R_t, t_);
// //             cout << "size :::::::::::: " << t->nrows << ", " << t->ncols << endl;;
//             
//             cout << "x y z **** : " << MATD_EL(t, 0, 0) * 0.05 <<  ",\t" << MATD_EL(t, 1, 0) * 0.05 << ",\t" << MATD_EL(t, 2, 0) * 0.05 << endl;
//             
// //             cout << "goodness : " << det->goodness << endl;
//             
//             matd_destroy(R);
//             matd_destroy(t_);
//             matd_destroy(t);
//             matd_destroy(R_t);
//             matd_destroy(pose);
            
            chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
            getRelativeTransform(0.1, 2.4908279215754123e+03, 2.4935314568583112e+03, 3.4745731382095448e+02, 2.4094331871742105e+02, det->p);
            chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
            chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>> (t2 - t1);
            std::cout << "opencv pnp cost time: " << time_used.count() << std::endl;
            
            t1 = chrono::steady_clock::now();
            getRelativeTransform2(0.1, 2.4908279215754123e+03, 2.4935314568583112e+03, 3.4745731382095448e+02, 2.4094331871742105e+02, det->p);
            t2 = chrono::steady_clock::now();
            chrono::duration<double> time_used2 = chrono::duration_cast<chrono::duration<double>> (t2 - t1);
            std::cout << "ippe cost time: " << time_used2.count() << std::endl;
            
            std::cout << "ippe speed / opencv speed :" << time_used.count() / time_used2.count() << endl;
            
            apriltag_detection_destroy(det);
        }
        
        zarray_destroy(detections);

        imshow("Tag Detections", frame);
        if (waitKey(1) == 27)
            break;
    }

    apriltag_detector_destroy(td);
    
    tags_destroy(tf);

    return 0;
}

void getRelativeTransform(double tag_size, double fx, double fy, double px, double py, double p[4][2])
{
    std::vector<cv::Point3f> objPts;
    std::vector<cv::Point2f> imgPts;
    double s = tag_size/2.;
    objPts.push_back(cv::Point3f(-s,-s, 0));
    objPts.push_back(cv::Point3f( s,-s, 0));
    objPts.push_back(cv::Point3f( s, s, 0));
    objPts.push_back(cv::Point3f(-s, s, 0));

    imgPts.push_back(cv::Point2f(p[0][0], p[0][1]));
    imgPts.push_back(cv::Point2f(p[1][0], p[1][1]));
    imgPts.push_back(cv::Point2f(p[2][0], p[2][1]));
    imgPts.push_back(cv::Point2f(p[3][0], p[3][1]));

//     cout << imgPts << endl;
    
    cv::Mat rvec, tvec;
    cv::Matx33d cameraMatrix(
                            fx, 0, px,
                            0, fy, py,
                            0,  0,  1);
//     cv::Vec4d distParam(0,0,0,0);
    double dist[] = {-5.0968287369808518e-02, -8.0252844113471298e+01, -1.5334326534795978e-03, -1.8098396142340031e-02,  -1.0045140113684745e+00};
    cv::Mat distParam = cv::Mat(1, 5, CV_64FC1, dist);
    cv::solvePnP(objPts, imgPts, cameraMatrix, distParam, rvec, tvec);
    cv::Matx33d r;
    cv::Rodrigues(rvec, r);
    
//   Eigen::Matrix3d wRo;
//   wRo << r(0,0), r(0,1), r(0,2), r(1,0), r(1,1), r(1,2), r(2,0), r(2,1), r(2,2);
// 
//   Eigen::Matrix4d T; 
//   T.topLeftCorner(3,3) = wRo;
//   T.col(3).head(3) << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//   T.row(3) << 0,0,0,1;
// 
//   return T;
    
    cv::Vec3d t(tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2));
    
    /* T.inv -> t = -r.inv * t */
//     cv::Matx44d T(r(0, 0), r(0, 1), r(0, 2), tvec.at<double>(0),
//                   r(1, 0), r(1, 1), r(1, 2), tvec.at<double>(1),
//                   r(2, 0), r(2, 1), r(2, 2), tvec.at<double>(2),
//                         0,       0,       0,                 1);
//     cv::Matx44d Tinv = T.inv();
//     cout << "x y z @@@@@@@ :" << Tinv(0, 3) <<  ",\t" << Tinv(1, 3) << ",\t" << Tinv(2, 3) << endl;
    cv::Vec3d ttt = - r.t() * t;
    
    cout << "x y z #### : " << ttt[0] <<  ",\t" << ttt[1] << ",\t" << ttt[2] << endl;
  
//     cout << "x y z ------ : " << tvec.at<double>(0) <<  ",\t" << tvec.at<double>(1) << ",\t" << tvec.at<double>(2) << endl;
//     cout << "R ------ : " << r << endl;
//     cout << "yaw ------ :" << std::atan2<float>(r(1, 0), r(0, 0)) * 57.3 << endl;    
}

void getRelativeTransform2(double tag_size, double fx, double fy, double cx, double cy, double p[4][2])
{
    matd_t *imagePoints;
    matd_t *cameraMatrix, *distParam;
    
    imagePoints = matd_create_data(4, 2, p[0]);
    
    cameraMatrix = matd_identity(3);
    MATD_EL(cameraMatrix, 0, 0) = fx;
    MATD_EL(cameraMatrix, 0, 2) = cx;
    MATD_EL(cameraMatrix, 1, 1) = fy;
    MATD_EL(cameraMatrix, 1, 2) = cy;
    
    double dist[] = {-5.0968287369808518e-02, -8.0252844113471298e+01, -1.5334326534795978e-03, -1.8098396142340031e-02,  -1.0045140113684745e+00};
    distParam = matd_create_data(1, 5, dist);
    
    matd_t *R = NULL, *t = NULL;
    float reprojErr;
    
    ippe::solvePoseOfMarker(tag_size, imagePoints, cameraMatrix, distParam, R, t, reprojErr);
    
    matd_t *tinv = matd_op("-M'*M", R, t);
    cout << "x y z ^^^^ : " << MATD_EL(tinv, 0, 0) <<  ",\t" << MATD_EL(tinv, 1, 0) << ",\t" << MATD_EL(tinv, 2, 0) << endl;
    
    matd_destroy(imagePoints);
    matd_destroy(cameraMatrix);
    matd_destroy(distParam);
    matd_destroy(tinv);
}

























