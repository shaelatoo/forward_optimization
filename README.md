# forward_optimization

The current version of this software is very, very untested, but the goal is to adapt some code I have been working on 
that optimizes PFSS coronal magnetic field models using constraints derived from images of the corona.  In the existing version 
(not on Github) the image-based constraints are derived from image analysis that searches for quasi-radial intensity gradients.
This optimization software takes those image-based constraints and compares the azimuth angle of the intensity gradient to 
the local magnetic field direction projected onto the image plane.  (See publication in ApJ April 2016 by Jones, Davila, and 
Uritsky.)

In this new version of the optimization I am using the FORWARD library from SolarSoft (IDL) to create a forward-modeled image
based on a PFSS/hydrostatic model.  The model images are then compared directly to the actual images, avoiding the image analysis
and the hairy geometry required to sample the model magnetic field at a particular point in a particular plane, project it 
onto the plane, and calculate the azimuth angle.  This will also eventually allow us to use images other than pB white light
images, as FORWARD can model pB, B and EUV emission.

