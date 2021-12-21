#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    double kenergy;
    Eigen::Matrix1212d M0;
    mass_matrix_linear_tetrahedron(M0,qdot,element,density,volume);
    Eigen::Vector12d qdotTetra;
    qdotTetra << qdot(3*element(0)) ,qdot(3*element(0)+1), qdot(3*element(0)+2),
            qdot(3*element(1)),qdot(3*element(1)+1),qdot(3*element(1)+2),
            qdot(3*element(2)),qdot(3*element(2)+1), qdot(3*element(2)+2),
            qdot(3*element(3)),qdot(3*element(3)+1),qdot(3*element(3)+2);
    kenergy = 1.0/2.0 * qdotTetra.transpose() * M0 *qdotTetra;
    T = kenergy;
}