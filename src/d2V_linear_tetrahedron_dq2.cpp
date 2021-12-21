#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include<iostream>
void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
      //Code to compute non-integrated hessian matrix goes here
       Eigen::Matrix43d DM;
       dphi_linear_tetrahedron_dX(DM,V,element,X);
       Eigen::Vector9d c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;
       c1 << DM(0,0), DM(0,1), DM(0,2),0,0,0,0,0,0;
       c2 << 0,0,0,DM(0,0), DM(0,1), DM(0,2),0,0,0;
       c3 << 0,0,0,0,0,0,DM(0,0), DM(0,1), DM(0,2);
       c4 << DM(1,0), DM(1,1), DM(1,2),0,0,0,0,0,0;
       c5 << 0,0,0,DM(1,0), DM(1,1), DM(1,2),0,0,0;
       c6 << 0,0,0,0,0,0,DM(1,0), DM(1,1), DM(1,2);
       c7 << DM(2,0), DM(2,1), DM(2,2),0,0,0,0,0,0;
       c8 << 0,0,0,DM(2,0), DM(2,1), DM(2,2),0,0,0;
       c9 << 0,0,0,0,0,0,DM(2,0), DM(2,1), DM(2,2);
       c10 << DM(3,0), DM(3,1), DM(3,2),0,0,0,0,0,0;
       c11 << 0,0,0,DM(3,0), DM(3,1), DM(3,2),0,0,0;
       c12 << 0,0,0,0,0,0,DM(3,0), DM(3,1), DM(3,2);
       Eigen::MatrixXd B(9,12);
       B.col(0) = c1;
       B.col(1) = c2;
       B.col(2) = c3;
       B.col(3) = c4;
       B.col(4) = c5;
       B.col(5) = c6;
       B.col(6) = c7;
       B.col(7) = c8;
       B.col(8) = c9;
       B.col(9) = c10;
       B.col(10) = c11;
       B.col(11) = c12;

       Eigen::Matrix3d F;
       Eigen::Vector3d X0,X1, X2, X3;
       X0 << q(3*element(0)),q(3*element(0)+1),q(3*element(0)+2);
       X1 << q(3*element(1)),q(3*element(1)+1),q(3*element(1)+2);
       X2 << q(3*element(2)),q(3*element(2)+1),q(3*element(2)+2);
       X3 << q(3*element(3)),q(3*element(3)+1),q(3*element(3)+2);


       Eigen::Matrix34d tetraIndice;
       tetraIndice.col(0) = X0;
       tetraIndice.col(1) = X1;
       tetraIndice.col(2) = X2;
       tetraIndice.col(3) = X3;
       F = tetraIndice*DM;
       Eigen::Matrix99d nineD;
       d2psi_neo_hookean_dF2(nineD,F,C,D);

       Eigen::Matrix1212d deriv;
       deriv = B.transpose()*nineD*B;
       dV = deriv;
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
