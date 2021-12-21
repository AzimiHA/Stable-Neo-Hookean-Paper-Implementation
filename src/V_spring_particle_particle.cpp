#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    V = 0.0;
    Eigen::VectorXd q(6,1);
    q = Eigen::VectorXd::Zero(6,1);

    q(0,0) = q0(0,0);
    q(1,0) = q0(1,0);
    q(2,0) = q0(2,0);
    q(3,0) = q1(0,0);
    q(4,0) = q1(1,0);
    q(5,0) = q1(2,0);

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
    double l = sqrt(q.transpose()*B.transpose()*B*q);
    V = (1.0/2.0)*stiffness*(l-l0)*(l-l0);
    
}