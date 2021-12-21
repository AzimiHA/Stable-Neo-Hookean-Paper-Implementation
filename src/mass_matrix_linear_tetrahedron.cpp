 
 #include <mass_matrix_linear_tetrahedron.h>

#include <iostream>
 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    Eigen::Matrix1212d Mtemp;

     Eigen::Matrix3d identity;
     identity.setIdentity();

     Mtemp << identity,identity,identity,identity,identity,identity,identity,identity,identity,identity,identity,
             identity,identity,identity,identity,identity;
     for(int i =0;i<12;i++){
         Mtemp(i,i)+=1;
     }
     M = Mtemp*density*volume/20.0;

 }