
#include <iostream>
#include <chrono>

#include <opencv2/opencv.hpp>

#include "apriltag.h"

#include "tag_family.h"

#include "ippe/ippe.h"

#include "lpe.h"

using namespace std;
using namespace cv;
using namespace apriltag;
using namespace vipad;

int main(int argc, char **argv)
{
    CvMat *g_mat = (CvMat *)cvLoad("marker.yaml");
    Mat gray = cvarrToMat(g_mat);
    
//     Mat src_gray = imread("pyramid_level1_320_240_2.png");
    
    /* -------------------------------------- */
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    pyrDown(gray, gray, Size(gray.cols / 2, gray.rows / 2));
    cout << gray.size() << endl;
    
    imshow("marker_gray", gray);
    
    vipad::lpe_params param;
    localPosition locat;
    
    param.angle = ClockwiseAngle_270;
    
//     param.fx = 2.4908279215754123e+03;
//     param.fy = 2.4935314568583112e+03;
//     param.cx = 3.4745731382095448e+02;
//     param.cy = 2.4094331871742105e+02;
//     
//     param.k1 = -5.0968287369808518e-02;
//     param.k2 = -8.0252844113471298e+01;
//     param.p1 = -1.5334326534795978e-03;
//     param.p2 = -1.8098396142340031e-02;
//     param.k3 = -1.0045140113684745e+00;
    
    param.fx = 4.4208943041154714e+02 / 2;
    param.fy = 4.4241613130483989e+02 / 2;
    param.cx = 3.2103104558965657e+02 / 2;
    param.cy = 2.4049214702014103e+02 / 2;
    
    param.k1 = 2.2974810891907303e-02;
    param.k2 = -1.6576845327552969e-01;
    param.p1 = -4.0782409518393070e-05;
    param.p2 = -5.7000919704701235e-03;
    param.k3 = 2.2850168365880927e-01;
    
    param.tag_type = tag36h11;
    param.marker_length = 0.1f;
    
    param.q = NULL;
    param.locat = &locat;

    param.width = gray.cols;
    param.height = gray.rows;
    param.input = gray.data;
        
    LocalPositionEstimation est(&param);
    
    est.estimateLocalPosition();
    
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<float> time_used = chrono::duration_cast<chrono::duration<float>> (t2 - t1);
    std::cout << "cost time: " << time_used.count() << std::endl;

    waitKey(0);
    
    CvMat mat = gray;
    cvSave("marker_pyramid.yaml", &mat);
    
//     cvReleaseMat(&g_mat);
    
    return 0;
}





















