#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    //First we make the q matrix by doing {q0,q1) into a 6x1 matrix
    Eigen::VectorXd q(6,1);
    q(0,0) = q0(0,0);
    q(1,0) = q0(1,0);
    q(2,0) = q0(2,0);
    q(3,0) = q1(0,0);
    q(4,0) = q1(1,0);
    q(5,0) = q1(2,0);

    //Then we make the B matrix, which is a 3x6 matrix
    Eigen::MatrixXd B(3,6);
    B(0,0) = -1.0;
    B(0,1) = 0.0;
    B(0,2) = 0.0;
    B(0,3) = 1.0;
    B(0,4) = 0.0;
    B(0,5) = 0.0;
    B(1,0) = 0.0;
    B(1,1) = -1.0;
    B(1,2) = 0.0;
    B(1,3) = 0.0;
    B(1,4) = 1.0;
    B(1,5) = 0.0;
    B(2,0) = 0.0;
    B(2,1) = 0.0;
    B(2,2) = -1.0;
    B(2,3) = 0.0;
    B(2,4) = 0.0;
    B(2,5) = 1.0;

    Eigen::Vector6d result;
    result = Eigen::Vector6d::Zero(6,1);
    double a =q.transpose()*B.transpose()*B*q;
    result = stiffness*(sqrt(a)-l0)*(1/sqrt(a))*B.transpose()*B*q;
    f = result;
}