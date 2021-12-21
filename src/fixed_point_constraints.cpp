#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int n = q_size/3;
    int m = indices.size();
    Eigen::SparseMatrixd temp(3*(n-m),3*n);
    int numFixed = 0;

    for(int i =0;i<n;i++){
        if(!(std::find(indices.begin(), indices.end(),i)!=indices.end())){
            tripletList.emplace_back(3*i-3*numFixed,3*i,1.0);
            tripletList.emplace_back(3*i-3*numFixed+1,3*i+1,1.0);
            tripletList.emplace_back(3*i-3*numFixed+2,3*i+2,1.0);
        }
        else{
            numFixed+=1;
        }
    }
    temp.setFromTriplets(tripletList.begin(), tripletList.end());
    P = temp;

}