#! /usr/bin/env python

# Author: Pierre Poulain
# Contributors: Justine Guegan, Edithe Selwa, Spencer Bliven
 
# This Python script computes principal axes from a PDB file

#==========================================================================
# import required modules
#==========================================================================
import sys
import os.path
import numpy
from pymol import cmd, stored
from pymol.cgo import *

#--------------------------------------------------------------------------
# compute principal axes
#--------------------------------------------------------------------------
def computeprincipalaxes(coord):
    """
    xyz: numpy nx3 array giving the input points
    returns: axis1, axis2, axis3, and center point
    """

    # compute geometric center
    center = numpy.mean(coord, 0)
    #print "geometric center coordinates:\n", center

    # center with geometric center
    coord = coord - center

    # compute principal axis matrix
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)
    # warning eigen values are not necessary ordered!
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html

    #--------------------------------------------------------------------------
    # order eigen values (and eigen vectors)
    #
    # axis1 is the principal axis with the biggest eigen value (eval1)
    # axis2 is the principal axis with the second biggest eigen value (eval2)
    # axis3 is the principal axis with the smallest eigen value (eval3)
    #--------------------------------------------------------------------------
    for i in xrange(len(e_values)):
            # find biggest eigen value
            if e_values[i] == max(e_values):
                    eval1 = e_values[i]
                    axis1 = e_vectors[:,i]
            # find smallest eigen value
            elif e_values[i] == min(e_values):
                    eval3 = e_values[i]
                    axis3 = e_vectors[:,i]
            # middle eigen value
            else:
                    eval2 = e_values[i]
                    axis2 = e_vectors[:,i]

    # Normalize axes to the bounds of the coordinates
    proj1 = numpy.dot(coord,axis1)
    range1 = numpy.max(numpy.abs(proj1))
    proj2 = numpy.dot(coord,axis2)
    range2 = numpy.max(numpy.abs(proj2))
    proj3 = numpy.dot(coord,axis2)
    range3 = numpy.max(numpy.abs(proj3))
    return axis1*range1,axis2*range2,axis3*range3,center

def principalaxes(selection,name="axes",scale_factor=1,state=1,radius=.5,hlength=-1,hradius=-1):
    """Display the principal axes of the selection as arrows
    """

    xyz = []
    cmd.iterate_state(state,selection,"xyz.append( (x,y,z) )", space={'xyz':xyz} )

    #create coordinates array
    coord = numpy.array(xyz, float)

    axis1,axis2,axis3,center = computeprincipalaxes(coord)

    if hlength < 0:
        hlength = radius * 3
    if hradius < 0:
        hradius = hlength * .6

    start1 = -scale_factor * axis1 + center
    end1 = scale_factor * axis1 + center
    mid1 = -axis1/numpy.linalg.norm(axis1)*hlength + end1

    start2 = -scale_factor * axis2 + center
    end2 = scale_factor * axis2 + center
    mid2 = -axis2/numpy.linalg.norm(axis2)*hlength + end2

    start3 = -scale_factor * axis3 + center
    end3 = scale_factor * axis3 + center
    mid3 = -axis3/numpy.linalg.norm(axis3)*hlength + end3

    axis1 =  [
            CYLINDER, start1[0],start1[1],start1[2],
            mid1[0],mid1[1],mid1[2],
            radius,
            1.0,0.0,0.0,
            1.0,0.0,0.0,
            CONE, mid1[0],mid1[1],mid1[2],
            end1[0],end1[1],end1[2],
            hradius, 0.0,
            1.0,0.0,0.0,
            1.0,0.0,0.0,
            1.0,0.0
            ]
    axis2 =  [
            CYLINDER, start2[0],start2[1],start2[2],
            mid2[0],mid2[1],mid2[2],
            radius,
            0.0,1.0,0.0,
            0.0,1.0,0.0,
            CONE, mid2[0],mid2[1],mid2[2],
            end2[0],end2[1],end2[2],
            hradius, 0.0,
            0.0,1.0,0.0,
            0.0,1.0,0.0,
            1.0,0.0
            ]
    axis3 =  [
            CYLINDER, start3[0],start3[1],start3[2],
            mid3[0],mid3[1],mid3[2],
            radius,
            0.0,0.0,1.0,
            0.0,0.0,1.0,
            CONE, mid3[0],mid3[1],mid3[2],
            end3[0],end3[1],end3[2],
            hradius, 0.0,
            0.0,0.0,1.0,
            0.0,0.0,1.0,
            1.0,0.0
            ]
    cmd.load_cgo(axis1+axis2+axis3, name)

cmd.extend("principalaxes",principalaxes)
cmd.auto_arg[0]["principalaxes"] = [cmd.object_sc,"selection",", "]


