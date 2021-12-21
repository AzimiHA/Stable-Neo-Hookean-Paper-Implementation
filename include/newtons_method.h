#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
        //Perform newton's method and return q*, so that qdot(t+1) = qdot(t+1) + q*
        //std::cout<<"Newton pls"<<"\n";
        Eigen::VectorXd negative(tmp_g.size());

        double alpha = 1.0;
        double p= 0.5;
        double c = 1e-8;
        Eigen::VectorXd d(tmp_g.size());
        int i = 0;
        g(tmp_g,x0);
        while  (((tmp_g).norm()>c) && (i<maxSteps)){
            H(tmp_H,x0);//Gets the H, assuming this is calculated using the passed in function in implicit
            g(tmp_g,x0); //Gets the g, assuming this is calculated using the passed in function in implicit
            negative=-1*tmp_g;
            Eigen::SimplicialLDLT<Eigen::SparseMatrixd > solverb(tmp_H);
            if(solverb.info()!=Eigen::Success) {
                std::cout<<"Couldn't compute A";
            }
            d = solverb.solve(negative);
            if(solverb.info()!=Eigen::Success) {
                std::cout<<"Couldn't solve";
            }
            alpha=1.0;
            while( !((f(x0+alpha*d) <= f(x0) + c*d.transpose()*tmp_g) || (alpha<1e-3))){
                alpha = p*alpha;
            }
            x0 = x0+alpha*d;
            i+=1;
        }
        return 0.0;
        //Finished newton's step
}
