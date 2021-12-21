#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    int n = qdot.rows()/3;
    //Want to construct selection matrix for each tetra's mass matrix
    typedef Eigen::Triplet<double> Tr;
    Eigen::SparseMatrixd Mass(3*n,3*n);
    //Thus we iterate over T (i.e m times as there are m tetras)
    std::vector<Tr> tripletList;

    for(int a =0; a<T.rows();a++){
        Eigen::RowVectorXi tetraIndicies(1,4);
        tetraIndicies = T.row(a);

        Eigen::Matrix1212d Mj;
        mass_matrix_linear_tetrahedron(Mj,qdot,tetraIndicies,density, v0(a));
        for(int i =0;i<12;i++){
            for(int j =0;j<12;j++){
                if(i<3){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(0)+j,Mj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(1)-1*(3-j),Mj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(2)-1*(6-j),Mj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(0)+i,3*tetraIndicies(3)-1*(9-j),Mj(i,j));
                    }
                }
                else if(2<i && i<6){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(0)+j,Mj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(1)-1*(3-j),Mj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(2)-1*(6-j),Mj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(1)-1*(3-i),3*tetraIndicies(3)-1*(9-j),Mj(i,j));
                    }
                }
                else if(5<i && i<9){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(0)+j,Mj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(1)-1*(3-j),Mj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(2)-1*(6-j),Mj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(2)-1*(6-i),3*tetraIndicies(3)-1*(9-j),Mj(i,j));
                    }
                }
                else if(8<i && i<12){
                    if(j<3){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(0)+j,Mj(i,j));
                    }
                    else if(2<j && j<6){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(1)-1*(3-j),Mj(i,j));
                    }
                    else if(5<j && j<9){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(2)-1*(6-j),Mj(i,j));
                    }
                    else if(8<j && j<12){
                        tripletList.emplace_back(3*tetraIndicies(3)-1*(9-i),3*tetraIndicies(3)-1*(9-j),Mj(i,j));
                    }
                }

            }
        }

    }

    //Finish iterating over T.
    Mass.setFromTriplets(tripletList.begin(), tripletList.end());
    M=Mass;


}