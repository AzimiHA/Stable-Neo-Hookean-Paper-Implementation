#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the energy density
//Output:
//  dphi - the 4x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    //Alternative way to calculate
    Eigen::Vector3d X02,X12, X22, X32;
    X02 = V.row(element(0));
    X12 = V.row(element(1));
    X22 = V.row(element(2));
    X32 = V.row(element(3));

    Eigen::Vector3d dX,dY,dZ;
    dX = X12-X02;
    dY = X22-X02;
    dZ = X32-X02;

    Eigen::Matrix3d T2;
    T2.col(0) = dX;
    T2.col(1) = dY;
    T2.col(2) = dZ;


    Eigen::Vector3d ones;
    ones<<-1,-1,-1;
    Eigen::Matrix3d T2inv;
    T2inv = T2.inverse();

    Eigen::Matrix43d dp;
    dp.row(0) = ones.transpose()*T2inv;
    dp.row(1) = T2inv.row(0);
    dp.row(2) = T2inv.row(1);
    dp.row(3) = T2inv.row(2);

    dphi = dp;


}