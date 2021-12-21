#include <assemble_stiffness.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) {

    int n = qdot.rows()/3;
    typedef Eigen::Triplet<double> Tr;
    Eigen::SparseMatrixd Stiff(3*n,3*n);
    //Thus we iterate over T (i.e m times as there are m tetras)
    std::vector<Tr> tripletList;

    for(int a =0; a<T.rows();a++){
        Eigen::RowVectorXi tetraIndicies(1,4);
        tetraIndicies = T.row(a);
        Eigen::Matrix1212d Hj;
        d2V_linear_tetrahedron_dq2(Hj,q,V,tetraIndicies, v0(a),C,D);
        //Goes over each entry in the Hj matrix and based on what block of the Hj matrix the entry is,
        //Adds it to the global stiffness matrix
        for(int i =0;i<12;i++){
            for(int j =0;j<12;j++){
                if(i<3){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(0)+j,Hj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(1)-1*(3-j),Hj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(2)-1*(6-j),Hj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(3)-1*(9-j),Hj(i,j));
                    }
                }
                else if(2<i && i<6){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(0)+j,Hj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(1)-1*(3-j),Hj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(2)-1*(6-j),Hj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(3)-1*(9-j),Hj(i,j));
                    }
                }
                else if(5<i && i<9){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(0)+j,Hj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(1)-1*(3-j),Hj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(2)-1*(6-j),Hj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(3)-1*(9-j),Hj(i,j));
                    }
                }
                else if(8<i && i<12){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(0)+j,Hj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(1)-1*(3-j),Hj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(2)-1*(6-j),Hj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(3)-1*(9-j),Hj(i,j));
                    }
                }

            }
        }

    }
    //Finish iterating over T.
    Stiff.setFromTriplets(tripletList.begin(), tripletList.end());
    K=-Stiff;



    };
