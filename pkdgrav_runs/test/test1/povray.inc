/*
** povray.inc -- DCR 94-10-19
** ==========================
** Definitions for use with POV-Ray scene files.
*/

#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "stones.inc"

/* definitions */

#declare ParticleFinish = finish { phong 0.2 diffuse 1.0 ambient 0.02 }

#declare ParticleSolid = texture { finish { ParticleFinish } }
#declare ParticleFancy = texture { T_Stone5 finish { ambient 0.5 } }
#declare ParticleBright = texture { finish { ambient 1.0 } }

#declare StyleParticle = texture { ParticleSolid }

// Some patterns for spheres to track orientations...
// (use rgbt or rgtf for transparency, allowing ball color to show through)

// 1. Beach Ball pattern...
#declare StyleOrient = texture { pigment { radial frequency 6 color_map { [0.0 color rgb <1,1,1>] [0.5 color rgb <1,1,1>] [0.5 color rgbt <0,0,0,1>] [1.0 color rgbt <0,0,0,1>] } } finish { ParticleFinish } }

// 2. Repeated Circular pattern...
// #declare StyleOrient = texture { pigment { boxed frequency 15 turbulence 0 color_map{[0.0 color rgb <0,0,0>] [0.7 color rgbf <1,1,1,1>] [1.0 color rgbf <1,1,1,1>] } } finish { ParticleFinish } }

// 3. Other texture ideas: http://www.f-lohmueller.de/pov_tut/tex/tex_750e.htm

#declare WallPlain = texture {pigment {transmit 1.0} finish {ambient 0.02}}
#declare WallAgate = texture {pigment {agate}}
#declare WallBrick = texture {pigment {brick brick_size 0.03 scale <1.5,1,1> mortar 0.003}}

/* defaults */

#declare StylePlane		= texture {WallPlain}
#declare StyleTriangle		= texture {WallPlain}
#declare StyleRectangle		= texture {WallPlain}
#declare StyleDisk		= texture {WallPlain}
#declare StyleCylinderInfinite	= texture {WallPlain}
#declare StyleCylinderFinite	= texture {WallPlain}
#declare StyleShell		= texture {WallPlain}

/* default background color */

background {color Black}

/* povray.inc */
