#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix
void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {

    int n = V.rows();
    int m = T.rows();
    int l = V_skin.rows();

    typedef Eigen::Triplet<double> Tr;
    Eigen::SparseMatrixd Skin(l,n);
    std::vector<Tr> tripletList;
    //First we iterate over the points/rows in V_skin
    for(int i=0;i<l;i++){
        //Then we iterate over the tetrahedrons
        for(int t=0;t<m;t++){
            //Check if the point i is in tetrahedron j
            //To do this we multiply we check the bounds using phi
            //Note psi(X) \in [0,1] means the point is in the tetrahedron
            Eigen::Vector4d inTetra;
            //Element is 4 vertex indices for tetrahedron
            Eigen::RowVectorXi element;
            element = T.row(t);
            //X is the point we want to check
            Eigen::Vector3d x;
            x = V_skin.row(i);
            phi_linear_tetrahedron(inTetra,V,element,x);
            if((0<=inTetra(0) && inTetra(0)<=1) && (0<=inTetra(1) && inTetra(1)<=1) && (0<=inTetra(2) && inTetra(2)<=1)
            && (0<=inTetra(3) && inTetra(3)<=1)){
                //Then point is in tetrahedron
                //So we insert the 4 points into sparse matrix
                tripletList.emplace_back(i, element(0), inTetra(0));
                tripletList.emplace_back(i, element(1), inTetra(1));
                tripletList.emplace_back(i, element(2), inTetra(2));
                tripletList.emplace_back(i, element(3), inTetra(3));
                break;
            }
        }
    }
    Skin.setFromTriplets(tripletList.begin(), tripletList.end());
    N = Skin;
}