#include <phi_linear_tetrahedron.h>
#include <iostream>
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the energy density
//Output:
//  phi - the 4x1 values the basis functions
void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    Eigen::Vector3d X0,X1, X2, X3;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    Eigen::Vector3d dX1, dX2, dX3;
    dX1 = X1-X0;
    dX2 = X2-X0;
    dX3 = X3-X0;

    Eigen::Matrix3d T;
    T.col(0) = dX1;
    T.col(1) = dX2;
    T.col(2) = dX3;

    Eigen::Vector3d phi3;
    phi3 = T.inverse()*(x-X0);
    double phi0 = 1-phi3(0)-phi3(1)-phi3(2);
    Eigen::Vector4d phi4;
    phi4(0) = phi0;
    phi4(1) = phi3(0);
    phi4(2) = phi3(1);
    phi4(3) = phi3(2);
    phi = phi4;
}