# Problem statement

The task is to reconstruct the 3D geometry of points by a set of point correspondences between two images. 

The goal of the triangulation method function is to compute the Fundamental F and Essential E matrices and subsequently decompose E into a rotation matrix R and a skew-symmetric matrix [t]x, that describe the relative pose of the cameras, using a set of known correspondences between two images. Once we obtain these, we proceed with the 3D reconstruction using the linear method for triangulation. 


# The team
Author: Ioanna Panagiotidou, Irina Gheorghiu, Cynthia Cai

Supervisor: Liangliang Nan

Date: May 27, 2022
