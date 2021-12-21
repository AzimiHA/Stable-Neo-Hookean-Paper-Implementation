#include <V_linear_tetrahedron.h>
#include <iostream>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron
void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::Matrix43d dphi;
        Eigen::Vector3d temp;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        Eigen::Matrix3d F;
        Eigen::Vector3d X0,X1, X2, X3;
        X0 << q(3*element(0)),q(3*element(0)+1),q(3*element(0)+2);
        X1 << q(3*element(1)),q(3*element(1)+1),q(3*element(1)+2);
        X2 << q(3*element(2)),q(3*element(2)+1),q(3*element(2)+2);
        X3 << q(3*element(3)),q(3*element(3)+1),q(3*element(3)+2);

        Eigen::Matrix34d tetraIndices;
        tetraIndices.col(0) = X0;
        tetraIndices.col(1) = X1;
        tetraIndices.col(2) = X2;
        tetraIndices.col(3) = X3;
        F = tetraIndices * dphi;
        double ps;
        psi_neo_hookean(ps, F, C, D);

        e = ps;
    };


    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  

}