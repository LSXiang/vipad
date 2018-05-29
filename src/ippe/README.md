# IPPE

Infinitesimal Plane-based Pose Estimation (IPPE) is a very fast and accurate way to compute a camera's pose
from a single image of a planar object using point correspondences. This has uses in several applications, 
including augmented reality, 3D tracking and pose estimation with planar markers, and 3D scene understanding.

This package is free and covered by the BSD licence without any warranty. We hope you find this code useful and if so
please cite our paper in your work:

    @article{ year={2014}, issn={0920-5691}, journal={International Journal of Computer Vision}, volume={109}, number={3},
              doi={10.1007/s11263-014-0725-5}, title={Infinitesimal Plane-Based Pose Estimation},
              url={http://dx.doi.org/10.1007/s11263-014-0725-5}, publisher={Springer US}, 
              keywords={Plane; Pose; SfM; PnP; Homography}, 
              author={Collins, Toby and Bartoli, Adrien}, 
              pages={252-286}, language={English} }

Please contact Toby (toby.collins@gmail.com) if you have any questions about the code, paper and IPPE.

You can go [here](https://github.com/tobycollins/IPPE) for a more detailed.

## Brief Introduction

The sources is implemented in C++ with no external dependencies, you only include the common directory.
