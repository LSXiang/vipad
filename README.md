# vipad
Visual-Inertial Positioning with AprilTag for Drones

---

# Overview
The main purpose of the project is to move away from opencv and migrate the MARKER-based visual positioning to CMU or DSP.

Programs in the src folder can be ported to the platform you need without relying on any other libraries. For different 
platforms, you may need to modify the dynamic memory request, specific modifications src/apriltag/inc/apriltag_allocator.h

Maybe you also need to add or modify the contents of TODO in lpe.cpp according to the actual situation.

Give applications example in app folder, maybe can help you.

I modified some programs and ported this function to the DSP core of RV1108 to successfully implement the marker visual 
positioning function. I hope the project can help you.

These codes are under MIT license. You don't need permission to use it or change it. If you have any questions about 
the code, please add an issue so I can see it.
