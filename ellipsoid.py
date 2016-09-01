# vim: filetype=python et sw=8 ts=8
"""Pymol functions for creating superquadric objects, including superellipsoids
and supertoroids.

Ellipsoids, toroids, and hyperboloids are 3D shapes with quadratic
cross-sections.  For instance, an ellipsoid consists of an ellipse rotated
around an orthogonal ellipse. A superellipsoid generallizes this to
non-elliptical cross sections, such as rounded rectangles and 4-pointed stars.

A good reference for these shapes is Chapter 2 of the book

        Segmentation and Recovery of Superquadrics
        by Ales Jaklic, Ales Leonardis and Franc Solina
        (Computational imaging and vision, Vol. 20, Kluwer, Dordrecth, 2000)
        Available at http://lrv.fri.uni-lj.si/~franc/SRSbook/SRS.html

History
=======

Original C++ code by Jonathan Metzgar
http://www.gamedev.net/page/resources/_/technical/opengl/superquadric-ellipsoids-and-toroids-opengl-lig-r1172

Python port by user Gilleain at
http://www.pymolwiki.org/index.php/Ellipsoid,

Further adaptations by Spencer Bliven
"""

from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
from pymol import cmd
import math

def signOfFloat(f):
        if f < 0: return -1
        if f > 0: return 1
        return 0

def sqC(v, n):
        """Signed power function of cosine. Cosine-like function, but with an
        additional polynomial power applied.
        """
        return signOfFloat(math.cos(v)) *  math.pow(math.fabs(math.cos(v)), n)

def sqCT(v, n, alpha):
        "Offset signed power function of cosine"
        return alpha + sqC(v, n)

def sqS(v, n):
        """Signed power function of sine. Sine-like function, but with an
        additional polynomial power applied.
        """
        return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)

def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e):
        """Convert a polar coordinate on a superellipsoid to a cartesian one

        x,y,z:      center position
        a1,a2,a3:   axes length
        u:          longitude to evaluate at (east-west)
        v:          latitude to evaluate at (north-south)
        n,e:        parameters of the superquadric ellipsoid

        returns a 6-tuple:
                x,y,z     position corresponding to (u,v)
                nx,ny,nz  normal vector at that point
        """
        x = a1 * sqC(u, n) * sqC(v, e) + x
        y = a2 * sqC(u, n) * sqS(v, e) + y
        z = a3 * sqS(u, n) + z
        nx = sqC(u, 2 - n) * sqC(v, 2 - e) / a1
        ny = sqC(u, 2 - n) * sqS(v, 2 - e) / a2
        nz = sqS(u, 2 - n) / a3
        return x, y, z, nx, ny, nz


def sqToroid(x, y, z, a1, a2, a3, u, v, n, e, alpha):
        """Convert a polar coordinate on a supertoroid to a cartesian one

        x,y,z:      center position
        a1,a2,a3:   axes length
        u:          longitude to evaluate at (east-west)
        v:          latitude to evaluate at (north-south)
        n,e:        parameters of the superquadric ellipsoid

        returns a 6-tuple:
                x,y,z     position corresponding to (u,v)
                nx,ny,nz  normal vector at that point
        """
        a1prime = 1.0 / (a1 + alpha)
        a2prime = 1.0 / (a2 + alpha)
        a3prime = 1.0 / (a3 + alpha)
        x = a1prime * sqCT(u, e, alpha) * sqC(v, n)
        y = a2prime * sqCT(u, e, alpha) * sqS(v, n)
        z = a3prime * sqS(u, e)
        nx = sqC(u, 2 - e) * sqC(v, 2 - n) / a1prime
        ny = sqC(u, 2 - e) * sqS(v, 2 - n) / a2prime
        nz = sqS(u, 2 - e) / a3prime
        return x, y, z, nx, ny, nz

def transformPoint(x,y,z,s,transformation):
        """transforms a point (x,y,z) according to the 16-element homogeneous
        transformation matrix.
        """
        if len(transformation) != 16:
                raise ValueError("Improperly formed transformation matrix")
        t = map(float,transformation)
        # Basic matrix multiplication
        f = t[12]*x + t[13]*y + t[14]*z + t[15] #should be 1 normally
        return (
                ( t[0]*x + t[1]*y + t[2]*z + t[3]*s )/f,
                ( t[4]*x + t[5]*y + t[6]*z + t[7]*s )/f,
                ( t[8]*x + t[9]*y +t[10]*z +t[11]*s )/f,
                )

def parseColor(color):
        """Parse various ways of specifying colors. Returns an array of 3-tuples

        Examples:
        "red" -> [(1,0,0)]
        "red green" -> [(1,0,0),(0,1,0)]
        (1,0,0) -> (1,0,0)
        [(1,0,0),(0,1,0)] -> [(1,0,0),(0,1,0)]

        """
        if type(color) == str:
                colors = [cmd.get_color_tuple(c) for c in color.split()]
                return colors

        colors = None
        if len(color) == 3:
                # Try treating as a 3-tuple
                try:
                        colors = [map(float,color)]
                except TypeError:
                        pass # not a 3-tuple
        if colors is None:
                # Assume list of 3-tuples
                if any([len(c) != 3 for c in color]):
                        raise ValueError("Unrecognized color format")
                try:
                        colors = [map(float,c) for c in color]
                except TypeError:
                        raise ValueError("Unrecognized color format")

        # validate color
        for col in colors:
                if any([c < 0 or 1 < c for c in col]):
                        raise ValueError("Color value out of bounds: %s"%col)
        return colors

def makeSuperQuadricEllipsoid(x, y, z,
                a1, a2, a3,
                n, e,
                u1=-math.pi/2, u2=math.pi/2,
                v1=-math.pi, v2=math.pi,
                u_segs=20, v_segs=40,
                color="gray",
                transformation=None,
                _self=cmd):
        """Create a superellipsoid, a generalization of the ellipsoid.

        This satisfies the equation:

                ( (x/a1)^(2/e) + (y/a2)^(2/e) )^(e/n) + (z/a3)^(2/n) = 1

        where e and n are the east-west and north-south exponents.
        This is an ellipsoid with n=1, e=1.

        As n approaches 0, the equatorial profile becomes more rectangular. At
        n=1, the equator is an ellipse. At the maximum of n=2, The equator
        forms a diamond. Higher values of n are forbidden as the would lead to
        concave shapes. The e exponent behaves equivalently for latitudinal
        cross-sections.

        Arguments
        =========
        x,y,z:          center coordinates
        a1,a2,a3:       length of axes
        n:              north-south exponent
        e:              east-west exponent
        u1,u2:          start and end angles in the east-west direction
        v1,v2:          start and end angles in the north-south direction
        u_segs:         Number of segments in the east-west direction
        v_segs:         Number of segments in the north-south direction
        color:          Either one (solid color) or six (-x,-y,-z,+x,+y,+z)
                        colors. Colors may be either a space-separated list of
                        pymol color names, or a list of tuples with RGB color
                        for the ellipsoid (0-1.0).
        transformation: a 16-element list giving the 4x4 transformation matrix,
                        in row-major order, as described in get_object_matrix()
                        {default: identity matrix}
        """
        # Cast arguments
        x = float(x)
        y = float(y)
        z = float(z)
        a1 = float(a1)
        a2 = float(a2)
        a3 = float(a3)
        n = float(n)
        e = float(e)
        u1 = float(u1)
        u2 = float(u2)
        v1 = float(v1)
        v2 = float(v2)
        u_segs = int(u_segs)
        v_segs = int(v_segs)

        # Parse color parameter
        colors = parseColor(color)
        if len(colors) == 1:
                colors = colors*6
        elif len(colors) != 6:
                raise ValueError("Unexpected number of colors")
        r,g,b = colors[0]


        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs

        o = [ BEGIN, TRIANGLES ]

        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # Compute points and normals for the square north-east of U,V
                        x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e)
                        x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e)
                        x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e)
                        x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e)

                        # Compute color
                        def blendcolor(nx,ny,nz):
                                #Normalize normal
                                l = math.sqrt(nx*nx+ny*ny+nz*nz)
                                fx = nx/l
                                fy = ny/l
                                fz = nz/l
                                #fx = signOfFloat(fx) * abs(fx)**1/2.
                                #fy = signOfFloat(fy) * abs(fy)**1/2.
                                #fz = signOfFloat(fz) * abs(fz)**1/2.
                                def pos(x):
                                        return x if x>=0 else 0

                                return [pos(colors[1][d]*fx) + pos(colors[0][d]*(-fx)) +
                                        pos(colors[3][d]*fy) + pos(colors[2][d]*(-fy)) +
                                        pos(colors[5][d]*fz) + pos(colors[4][d]*(-fz))
                                        for d in xrange(3)]

                        r1,g1,b1 = blendcolor(n1x,n1y,n1z)
                        r2,g2,b2 = blendcolor(n2x,n2y,n2z)
                        r3,g3,b3 = blendcolor(n3x,n3y,n3z)
                        r4,g4,b4 = blendcolor(n4x,n4y,n4z)

                        #Apply transformation, if any
                        if transformation != None:
                                x1,y1,z1 = transformPoint(x1,y1,z1,1,transformation)
                                x2,y2,z2 = transformPoint(x2,y2,z2,1,transformation)
                                x3,y3,z3 = transformPoint(x3,y3,z3,1,transformation)
                                x4,y4,z4 = transformPoint(x4,y4,z4,1,transformation)

                                n1x,n1y,n1z = transformPoint(n1x,n1y,n1z,0,transformation)
                                n2x,n2y,n2z = transformPoint(n2x,n2y,n2z,0,transformation)
                                n3x,n3y,n3z = transformPoint(n3x,n3y,n3z,0,transformation)
                                n4x,n4y,n4z = transformPoint(n4x,n4y,n4z,0,transformation)


                        # Right triangle with corner UV and legs going north and east
                        o.extend([COLOR, r1, g1, b1, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([COLOR, r2, g2, b2, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r4, g4, b4, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        # The other half of the rectangle, sharing a diagonal the first triangle
                        o.extend([COLOR, r2, g2, b2, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r3, g3, b3, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([COLOR, r4, g4, b4, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])

                        #o.extend([COLOR, r1, g1, b1, NORMAL, x1, y1, z1, VERTEX, n1x, n1y, n1z])
                        #o.extend([COLOR, r2, g2, b2, NORMAL, x2, y2, z2, VERTEX, n2x, n2y, n2z])
                        #o.extend([COLOR, r4, g4, b4, NORMAL, x4, y4, z4, VERTEX, n4x, n4y, n4z])
                        #o.extend([COLOR, r2, g2, b2, NORMAL, x2, y2, z2, VERTEX, n2x, n2y, n2z])
                        #o.extend([COLOR, r3, g3, b3, NORMAL, x3, y3, z3, VERTEX, n3x, n3y, n3z])
                        #o.extend([COLOR, r4, g4, b4, NORMAL, x4, y4, z4, VERTEX, n4x, n4y, n4z])

                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o

def makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, n, e,
                u1, u2, v1, v2, u_segs, v_segs, color="gray",
                transformation=None, _self=cmd):
        """Create a supertoroid, a generalization of the toroid.

        This satisfies the equation:

                ( [ (x/a1)^(2/e) + (y/a2)^(2/e) ]^(e/2) - alpha )^(2/n) + (z/a3)^(2/n) = 1

        where e and n are the east-west and north-south exponents.
        This is an ellipsoid with n=1, e=1.

        As n approaches 0, the equatorial profile becomes more rectangular. At
        n=1, the equator is an ellipse. At the maximum of n=2, The equator
        forms a diamond. Higher values of n are forbidden as the would lead to
        concave shapes. The e exponent behaves equivalently for latitudinal
        cross-sections.

        The radius of the supertoroid is given by
        
                R = a4 * sqrt( a1^2 + a2^2 )


        Arguments
        =========
        x,y,z:          center coordinates
        a1,a2,a3:       length of axes
        alpha:          radius parameter
        n:              north-south exponent
        e:              east-west exponent
        u1,u2:          start and end angles in the east-west direction
        v1,v2:          start and end angles in the north-south direction
        u_segs:         Number of segments in the east-west direction
        v_segs:         Number of segments in the north-south direction
        color:          Color, either as a pymol color names or a tuple with
                        RGB color for the toroid (0-1.0).
        transformation: a 16-element list giving the 4x4 transformation matrix,
                        in row-major order, as described in get_object_matrix()
                        {default: identity matrix}
        """
        # Cast arguments
        x = float(x)
        y = float(y)
        z = float(z)
        a1 = float(a1)
        a2 = float(a2)
        a3 = float(a3)
        alpha = float(alpha)
        n = float(n)
        e = float(e)
        u1 = float(u1)
        u2 = float(u2)
        v1 = float(v1)
        v2 = float(v2)
        u_segs = int(u_segs)
        v_segs = int(v_segs)

        # Parse color parameter
        r, g, b = parseColor(color)[0]

        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs

        o = [ BEGIN, TRIANGLES ]

        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = sqToroid(x, y, z, a1, a2, a3, U, V, n, e, alpha)
                        x2, y2, z2, n2x, n2y, n2z = sqToroid(x, y, z, a1, a2, a3, U + dU, V, n, e, alpha)
                        x3, y3, z3, n3x, n3y, n3z = sqToroid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e, alpha)
                        x4, y4, z4, n4x, n4y, n4z = sqToroid(x, y, z, a1, a2, a3, U, V + dV, n, e, alpha)

                        #Apply transformation, if any
                        if transformation != None:
                                x1,y1,z1 = transformPoint(x1,y1,z1,1,transformation)
                                x2,y2,z2 = transformPoint(x2,y2,z2,1,transformation)
                                x3,y3,z3 = transformPoint(x3,y3,z3,1,transformation)
                                x4,y4,z4 = transformPoint(x4,y4,z4,1,transformation)

                                n1x,n1y,n1z = transformPoint(n1x,n1y,n1z,0,transformation)
                                n2x,n2y,n2z = transformPoint(n2x,n2y,n2z,0,transformation)
                                n3x,n3y,n3z = transformPoint(n3x,n3y,n3z,0,transformation)
                                n4x,n4y,n4z = transformPoint(n4x,n4y,n4z,0,transformation)

                        o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])

                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o

def makeEllipsoid(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)

def ellipsoid(name,x,y,z, rx,ry,rz,segs=40,color="red",transformation=None,_self=cmd):
        """Create an ellipsoid

        Arguments
        =========
        x,y,z:          center coordinates
        rx,ry,rz:       length of axes
        n:              north-south exponent
        e:              east-west exponent
        u1,u2:          start and end angles in the east-west direction
        v1,v2:          start and end angles in the north-south direction
        u_segs:         Number of segments in the east-west direction
        v_segs:         Number of segments in the north-south direction
        color:          Either one (solid color) or six (-x,-y,-z,+x,+y,+z)
                        colors. Colors may be either a space-separated list of
                        pymol color names, or a list of tuples with RGB color
                        for the ellipsoid (0-1.0).
        transformation: a 16-element list giving the 4x4 transformation matrix,
                        in row-major order, as described in get_object_matrix()
                        {default: identity matrix}

        """
        cgo = makeSuperQuadricEllipsoid(x, y, z, rx,ry,rz, 1.0, 1.0,
                        u_segs=int(segs)/2,v_segs=segs,color=color,transformation=transformation)
        cmd.load_cgo(cgo,name)
cmd.extend("ellipsoid",ellipsoid)

def makeCylinder(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)

def makeSpindle(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)

def makeDoublePyramid(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 2.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)

def makePillow(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 0.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)

def makeRoundCube(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.2, 0.2, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)

def makeToroid(x, y, z, a1, a2, a3, alpha):
                return makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, 1.0, 1.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)

#x, y, z, rx, ry, rz = 1, 1, 1, 1, 2, 3
#cmd.load_cgo(makeEllipsoid(x, y, z, rx, ry, rz), 'ellipsoid-cgo')
#x, y, z, rx, ry, rz = 1, 1, 1, 8, 2, 2
#cmd.load_cgo(makeToroid(x, y, z, rx, ry, rz, 3), 'toroid-cgo')

def ellipsoid(name,x,y,z, rx,ry,rz,segs=40,color="red",transformation=None,_self=cmd):
        cgo = makeSuperQuadricEllipsoid(x, y, z, rx,ry,rz, 1.0, 1.0,
                        u_segs=int(segs)/2,v_segs=segs,color=color,transformation=transformation)
        cmd.load_cgo(cgo,name)
cmd.extend("ellipsoid",ellipsoid)


