#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    double J = F.determinant();
    Eigen::Matrix3d I3;
    I3 = F.transpose()*F;
    double a = 1.0+(2.0*C)/(2.0*D)-(2.0*C)/(8.0*D);
    double temp_psi = C*(I3.trace()-3.0)+D*pow((J-a),2)-C* log(I3.trace()+1);
    psi = temp_psi;
}