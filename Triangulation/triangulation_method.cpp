/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    /*/// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission */

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    int i;
    double smx0=0, smy0=0, smx1=0, smy1=0;
    double p1x, p1y, p2x, p2y;
    for(i=0;i<points_0.size();i++){
        smx0=smx0+points_0[i].x();
        smy0=smy0+points_0[i].y();
    }

    p1x=smx0/points_0.size();
    p1y=smy0/points_0.size();

    for(i=0;i<points_1.size();i++){
        smx1=smx1+points_1[i].x();
        smy1=smy1+points_1[i].y();
    }
    p2x=smx1/points_1.size();
    p2y=smy1/points_1.size();

    std::cout<<"new center of the 1st image "<<p1x<<" "<<p1y<<"\n";
    std::cout<<"new center of the 2nd image "<<p2x<<" "<<p2y<<"\n";

    double tx0=-p1x, ty0=-p1y, tx1=-p2x, ty1=-p2y;
    //Matrix33 T0(0, 0, tx0, 0, 0, ty0, 0, 0, 1);
    int j;
    double dist0=0.0, dist1=0.0;
    for( i=0;i<points_0.size();i++)
    { dist0=sqrt((points_0[i].x()-p1x)*(points_0[i].x()-p1x)+(points_0[i].y()-p1y)*(points_0[i].y()-p1y))+dist0;
    }

    double s0=sqrt(2)/(dist0/points_0.size());

    for( i=0;i<points_1.size();i++)
    {dist1=sqrt((points_1[i].x()-p2x )*(points_1[i].x()-p2x)+(points_1[i].y()-p2y)*(points_1[i].y()-p2y))+dist1;}

    double s1=sqrt(2)/(dist1/points_1.size());

    std::cout<<"scale of the 1st image "<<s0<<"\n";
    std::cout<<"scale of the 2nd image "<<s1<<"\n";

    Matrix33 T0(s0, 0, tx0, 0, s0, ty0, 0, 0, 1);
    Matrix33 T1(s1, 0, tx1, 0, s1, ty1, 0, 0, 1);

    std::cout<<"T0: \n"<<T0<<" \n";
    std::cout<<"T1: \n"<<T1;
    int d0=points_0.size();
    int d1=points_1.size();
    Matrix q00(3,1,0.0);
    Matrix q01(3, 1, 0.0);
    Matrix W0 (points_0.size(), 9, 0.0);
    Matrix W1 (points_1.size(), 9, 0.0);
    for(i=0;i<points_0.size();i++){
        //int rd=rand()%points_0.size();
        Matrix p0(3,1,0.0);
        p0[0][0]=points_0[i].x();
        p0[1][0]=points_0[i].y();
        p0[2][0]=1;

        q00=T0 * p0;
        //int k=0;

        W0[i][0]=q00[0][0]*p0[0][0];
        W0[i][1]=q00[0][0]*p0[1][0];
        W0[i][2]=q00[0][0]*p0[2][0];
        W0[i][3]=q00[1][0]*p0[0][0];
        W0[i][4]=q00[1][0]*p0[1][0];
        W0[i][5]=q00[1][0]*p0[2][0];
        W0[i][6]=q00[2][0]*p0[0][0];
        W0[i][7]=q00[2][0]*p0[1][0];
        W0[i][8]=q00[2][0]*p0[2][0];
        //k++;
        }

    for(i=0;i<points_1.size();i++){
        //int rd=rand()%points_1.size();
        Matrix p1(3,1,0.0);
        p1[0][0]=points_1[i].x();
        p1[1][0]=points_1[i].y();
        p1[2][0]=1;

        q01=T1 * p1;
        //int k=0;

        W1[i][0]=q01[0][0]*p1[0][0];
        W1[i][1]=q01[0][0]*p1[1][0];
        W1[i][2]=q01[0][0]*p1[2][0];
        W1[i][3]=q01[1][0]*p1[0][0];
        W1[i][4]=q01[1][0]*p1[1][0];
        W1[i][5]=q01[1][0]*p1[2][0];
        W1[i][6]=q01[2][0]*p1[0][0];
        W1[i][7]=q01[2][0]*p1[1][0];
        W1[i][8]=q01[2][0]*p1[2][0];
        //k++;
    }

    //std::cout<<"W0: \n"<<W0<<"\n";
    //std::cout<<"W1: \n"<<W1<<"\n";

    //SVD of W0 - USV
    Matrix U0(d0,d0,0.0);
    Matrix S0(d0,9,0.0);
    Matrix V0(9,9,0.0);
    svd_decompose(W0,U0,S0,V0);


    // take the last column of V as the solution for f
    Vector fs = V0.get_column(V0.cols()-1);
    //std::cout << "fs: \n" << fs << std::endl;

    Matrix33 Fq ;
    Fq.set_row(0, {fs[0], fs[1], fs[2]});
    Fq.set_row(1, {fs[3], fs[4], fs[5]});
    Fq.set_row(2, {fs[6], fs[7], fs[8]});

    //std::cout << "Fq: \n" << Fq << std::endl;

    // Enforce Rank 2 constraint by again using SVD decomposition
    Matrix Uq(3, 3, 0.0);
    Matrix Sq(3, 3, 0.0);
    Matrix Vq(3, 3, 0.0);
    svd_decompose(Fq, Uq, Sq, Vq);

    //std::cout << "Uq: \n" << Uq << std::endl;
    //std::cout << "Sq: \n" << Sq << std::endl;
    //std::cout << "Vq: \n" << Vq << std::endl;

    // enforce rank 2 --> manipulate the diagonal S
    Sq[2][2] = 0;
    // recompute Fq after enforcement
    Fq = Uq * Sq * Vq.transpose();
    //std::cout << "Fq: \n" << Fq << std::endl;

//    // TODO: ????????? Intermediate step --- Find the closest rank-2 matrix --- (ΤΙ ΕΙΝΑΙ ΑΥΤΟ;;;;)
//    // ΝΑ ΔΩ ΤΙ ΕΙΝΑΙ . ΑΝ ΤΟ ΕΧΩ ΚΑΝΕΙ Ή ΑΝ ΧΡΕΙΑΖΕΤΑΙ ΚΑΤΙ ΠΑΡΑΠΑΝΩ .
//
    // Last step: DENORMALIZATION Fq to be F. F = T′TFqT
    Matrix33 F = T1.transpose() * Fq * T0;
    std::cout << "F: \n" << F << std::endl;
//
//    // TODO: Intermediate step - The recovered F is up to scale. Please scale F such that F(2, 2) = 1.0 after denormalization.
//    // so probably take the last element and divide everything with that
//    F = F / F[2][2];
//    std::cout << "F: \n" << F << std::endl;
//    //std::cout << "Norm_F: " <<  norm(F) << std::endl;
//    //std::cout << "Norm_F: " <<  F[0][0] + F[0][1] + F[0][2] + F[1][0] + F[1][1] + F[1][2] + F[2][0] + F[2][1] + F[2][2] << std::endl;
//
//
//    // to be deleted
//    // TODO step 2: compute the essential matrix E; E = KTFK
//
    // consruct the K matrix --> same for both cameras
    Matrix33 K(fx, 0, cx, 0, fy, cy, 0, 0, 1);
    //std::cout << "K: \n" << K << std::endl;

    Matrix33 E = K.transpose() * F * K;  // 5 degrees of freedom --> it encodes R, and t (extrinsics).
    std::cout << "E: \n" << E << std::endl;

    // Decomposition of E into R and t.
    // we define two matrices that we will use in the decomposition of E --> W and Z
    Matrix33 We(0,-1,0,1,0,0,0,0,1);
    Matrix33 Ze(0,1,0,-1,0,0,0,0,0);

    Matrix Ue(3, 3, 0.0);
    Matrix Se(3, 3, 0.0);
    Matrix Ve(3, 3, 0.0);
    svd_decompose(E, Ue, Se, Ve);
    Matrix Tx(3, 3,0.0);
    Tx=Ue*Ze* transpose(Ue);
    Matrix Re1(3,3,0.0);
    Matrix Re2(3,3,0.0);
    //Re1=Ue*We*transpose(Ve);
    //Re2=Ue*transpose(We)*transpose(Ve);
    Re1=determinant(Ue*We*transpose(Ve))*Ue*We*transpose(Ve);
    Re2=determinant(Ue*transpose(We)*transpose(Ve))*Ue*transpose(We)*transpose(Ve);
    Matrix t1(3,1,0.0);
    t1[0][0]=-Ue[0][2];
    t1[1][0]=-Ue[1][2];
    t1[2][0]=-Ue[2][2];

    Matrix t2(3,1,0.0);
    t2[0][0]=Ue[0][2];
    t2[1][0]=Ue[1][2];
    t2[2][0]=Ue[2][2];
//UZUTt
    Matrix te1(3,1,0.0);
    te1=Tx*t1;
    Matrix te2(3,1,0.0);
    te2=Tx*t2;

//solutions of E

    Matrix A1(3,1,0.0);
    Matrix A2(3,1,0.0);
    Matrix A3(3,1,0.0);
    Matrix A4(3,1,0.0);

    A1=t1*Re1;
    A2=t1*Re2;
    A3=t2*Re1;
    A4=t2*Re2;




//
//    std::cout << "Ue: \n" << Ue << std::endl;
//    std::cout << "Se: \n" << Se << std::endl;
//    std::cout << "Ve: \n" << Ve << std::endl;




    //std::cout<<"S0:" <<S0;

    /*std::cout<<"new center of the 1st image "<<p1x<<" "<<p1y<<"\n";
    std::cout<<"new center of the 2nd image "<<p2x<<" "<<p2y<<"\n";*/


    return points_3d.size() > 0;



}
