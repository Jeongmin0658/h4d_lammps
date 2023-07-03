#############################
# Basic modules
#############################
#from __future__ import print_function
from ctypes import *        
import sys,random,math
import numpy as np
import os
from numpy import linalg as LA

#############################
# Import modules
#############################
from define_lmp_misc import LMP_misc
mymisc=LMP_misc()

class LMP_flying_init():
    def __init__(self):
        self.intro="Initial orientation/momentum for flying water"

    def rot_water(self,x):
        # ref: http://www.cognitive-antics.net/uniform-random-orientation/
        # ref: https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
        # draw three random numbers
        u=random.random()
        v=random.random()
        w=random.random()
        theta=2.*np.pi*u        # first - z
        phi=np.arcsin(v)        # second - y
        roll=np.pi*(2.*w-1.)    # third - x
        # rotation matrix
        Euler=self.euler_to_rotMat(roll, phi, theta)
        # rotate x vector
        x_rot=[]
        for i in np.arange(3): # three atoms
            xi=x[3*i:3*i+3]
            x_rot.extend(np.dot(Euler,xi))
        return x_rot

    def diagonalize_mom_iner(self,mass=[16,1,1],x_res=None):
        # x_res: residual position relavtive to com
        # mom_iner: matrix before diagonalization
        #from numpy import linalg as LA
        # calculate moment of inertia
        # 3 x 3 matrix
        mom_iner=[[0. for l in np.arange(3)] for m in np.arange(3)]
        # for each atom
        for i,k in enumerate(mass):
            xi_res=np.array(x_res[3*i:3*i+3])    # residual vector of i site
            xi_res_sq=np.dot(np.transpose(xi_res),xi_res)
            for l,rl in enumerate(xi_res):
                for m,rm in enumerate(xi_res):
                    mom_iner[l][m]-=k*rl*rm
                    if l==m:
                        mom_iner[l][m]+=k*xi_res_sq
        # diagonalize the moment of inertia
        mom_iner_diag, qmat = LA.eig(np.array(mom_iner))
        #def isRotationMatrix(R) :
        #    Rt = np.transpose(R)
        #    shouldBeIdentity = np.dot(Rt, R)
        #    I = np.identity(3 , dtype = R.dtype)
        #    n = np.linalg.norm(I - shouldBeIdentity)
        #    #print ("Norm:",n)
        #    return n < 1e-6
        #if not isRotationMatrix(qmat):
        #    print ("Wrong diagonalization of moment of inertia")
        #    print (format(qmat))
        #    print (format(mom_iner_diag))
        #    print (format(np.subtract(np.dot(qmat,np.dot(np.diag(mom_iner_diag),np.transpose(qmat))),mom_iner)))
        #    #exit()
        #else:
        #    print("Diagonal matrix : {} ".format(mom_iner_diag))
        return np.transpose(qmat),mom_iner_diag

    def get_residual_pos_water(self,mass=[16,1,1],x=[1.,1.,1.]):
        xcom=[[] for i in np.arange(3)] # center of mass
        totalmass=np.sum(mass)
        #print ("Total mass:",totalmass)
        for i,r in enumerate(x):
            axis=i%3        # axis 
            m=mass[i//3]    # mass for each atom
            xcom[axis].append(m*r)
        xcom=[np.sum(xcom[j])/totalmass for j in np.arange(3)]
        xc=xcom[:]
        xcom.extend(xc)
        xcom.extend(xc)
        #print ("COM position:",xcom)
        x_res=np.subtract(x,xcom)
        return x_res

    def random_atom_vel_rigid(self,mass=[16,1,1],x=[1.,1.,1.],trans=1.,rot=1.):
        # x: position of each atom in space frame
        # x_res: residual position of each atom relative to com (body frame)
        if len(x) != 9:
            print ("Should be a single water molecule")
            exit()
        x_res=self.get_residual_pos_water(mass,x)
        # First, diagonalize moment of inertia
        # qmat: rotational matrix
        # mom_iner_diag: diagonal matrix of moment of inertia
        qmat,mom_iner_diag=self.diagonalize_mom_iner(mass=mass,x_res=x_res)
        # Second, draw random velocities
        # trans: translational temperature
        # rot: rotational temperature
        # linear - three dofs
        vcom=[]
        totalmass=np.sum(mass)
        for j in np.arange(3):
            u=random.gauss(0.,1.)/math.sqrt(totalmass/trans)
            vcom.append(u)
        # angular - three dofs
        vang=[]
        for j in np.arange(3):
            v=random.gauss(0.,1.)/math.sqrt(mom_iner_diag[j]/rot)
            vang.append(v)
        # Finally, back to space frame for each atom
        # angular velocities in space frame
        vang_space=np.dot(qmat,vang)
        # velocity for each atom
        v_atom=(3*3*c_double)()
        for i in np.arange(3): # three atoms in a water molecule
            xi_res=x_res[3*i:3*i+3] # residual vector of i site
            angular=np.cross(vang_space,xi_res)
            v=np.add(vcom,angular)
            for j in np.arange(3):  # each axis
                v_atom[j+3*i]=v[j]
        if len(x) != len(v_atom):
            print ("Wrong velocity length")
            exit()
        return v_atom

    def euler_to_rotMat(self,yaw,pitch,roll):
        Rz_yaw = np.array([
            [np.cos(yaw), -np.sin(yaw), 0],
            [np.sin(yaw),  np.cos(yaw), 0],
            [          0,            0, 1]])
        Ry_pitch = np.array([
            [ np.cos(pitch), 0, -np.sin(pitch)],
            [             0, 1,             0],
            [ np.sin(pitch), 0, np.cos(pitch)]])
        Rx_roll = np.array([
            [1,            0,             0],
            [0, np.cos(roll), -np.sin(roll)],
            [0, np.sin(roll),  np.cos(roll)]])
        # R = RzRyRx
        rotMat = np.dot(Rz_yaw, np.dot(Ry_pitch, Rx_roll))
        return rotMat

