
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
    
    imshow("marker_gray", gray);
    
    vipad::lpe_params param;
    localPosition locat;
    
    param.angle = ClockwiseAngle_270;
    
    param.fx = 2.4908279215754123e+03;
    param.fy = 2.4935314568583112e+03;
    param.cx = 3.4745731382095448e+02;
    param.cy = 2.4094331871742105e+02;
    
    param.k1 = -5.0968287369808518e-02;
    param.k2 = -8.0252844113471298e+01;
    param.p1 = -1.5334326534795978e-03;
    param.p2 = -1.8098396142340031e-02;
    param.k3 = -1.0045140113684745e+00;
    
    param.tag_type = tag36h11;
    param.marker_length = 0.1f;
    
    param.q = NULL;
    param.locat = &locat;

    param.width = gray.cols;
    param.height = gray.rows;
    param.input = gray.data;
        
    LocalPositionEstimation est(&param);
    
    est.estimateLocalPosition();

    waitKey(0);
    
    cvReleaseMat(&g_mat);
    
    return 0;
}





















