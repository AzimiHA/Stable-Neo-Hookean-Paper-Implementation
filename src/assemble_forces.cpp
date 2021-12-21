#include <assemble_forces.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 


        int n = qdot.rows()/3;
        Eigen::VectorXd F(3*n,1);
        F.setZero();
        //Iterate over tetrahedrons
        Eigen::RowVectorXi tetraIndicies(1,4);
        for(int a = 0;a<T.rows(); a++){
            tetraIndicies = T.row(a);
            Eigen::Vector12d fj;
            dV_linear_tetrahedron_dq(fj,q,V,tetraIndicies,v0(a), C,D);
            //Iterate over Fj elements
            for(int i =0;i<4;i++){
                F(3*tetraIndicies(i)) += fj(3*i);
                F(3*tetraIndicies(i)+1) += fj(3*i+1);
                F(3*tetraIndicies(i)+2) += fj(3*i+2);
            }
        }
        f = -F;
    };