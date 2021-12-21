#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

//Assignment 3 code 
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <init_state.h>
#include <find_min_vertices.h>
#include <fixed_point_constraints.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <T_linear_tetrahedron.h>
#include <V_linear_tetrahedron.h>
#include <V_spring_particle_particle.h>
#include <dV_linear_tetrahedron_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <mass_matrix_mesh.h>
#include <assemble_forces.h>
#include <assemble_stiffness.h>
#include <linearly_implicit_euler.h>
#include <implicit_euler.h>
#include <build_skinning_matrix.h>

#include <chrono>
#include <thread>


//Variable for geometry
Eigen::MatrixXd V; //vertices of simulation mesh 
Eigen::MatrixXi T; //faces of simulation mesh
Eigen::MatrixXi F; //faces of simulation mesh

//variables for skinning
Eigen::MatrixXd V_skin;
Eigen::MatrixXi F_skin;
Eigen::SparseMatrixd N; 

//material parameters
double density = 0.1;
double YM = 6e5; //young's modulus
double mu = 0.40; //poissons ratio, was 0.4
double D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
double C = 0.5*YM/(2.0*(1.0+mu));

//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::SparseMatrixd P;
Eigen::VectorXd x0; 

//mass matrix
Eigen::SparseMatrixd M;
Eigen::VectorXd v0;

//scratch memory for assembly
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::SparseMatrixd tmp_stiffness;

std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points;

bool skinning_on = true;
bool fully_implicit = true;
bool bunny = true;
bool cube=false;
bool plane=false;
bool point = false;
bool line = false;
bool cube2=false;

//selection spring
double k_selected = 1e5;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {
    if(plane){
        for(int i =0; i<q.size()/3;i++){
            q(3*i+1)=0;
        }
        plane=!plane;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }else if(point){
        for(int i =0; i<q.size()/3;i++){
            q(3*i)=0;
            q(3*i+1)=0;
            q(3*i+2)=0;
        }
        point=!point;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }else if(line){
        for(int i =0; i<q.size()/3;i++){
            q(3*i)=0;
            q(3*i+1)=0;
        }
        line=!line;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }




    double V_ele, T_ele, KE,PE;

    spring_points.clear();

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);

    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
    }

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> qdot_1)->double {
        double E = 0;
        Eigen::VectorXd newq = P.transpose()*(q+dt*qdot_1)+x0;

        for(unsigned int ei=0; ei<T.rows(); ++ei) {
            
            V_linear_tetrahedron(V_ele,newq , V, T.row(ei), v0(ei), C, D);
            E += V_ele;
        }

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {   
            V_spring_particle_particle(V_ele, spring_points[pickedi].first, newq.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            E += V_ele;
        }

        E += 0.5*(qdot_1 - qdot).transpose()*M*(qdot_1 - qdot);

        return E;
    };

    auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) { 
            //std::cout<<"Before assemble forces";
            assemble_forces(f, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D);

            //std::cout<<"after assemble forces";
            for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
                dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose()*q2+x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
            }
            //std::cout<<"picked vertices";
            f = P*f;
        };

        //assemble stiffness matrix,
        auto stiffness = [&](Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) {
            //std::cout<<"before assembly stiff";
            assemble_stiffness(K, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D);
            //std::cout<<"after assemble stiff";
            K = P*K*P.transpose();
        };

        if(fully_implicit)
            implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_qdot, tmp_force, tmp_stiffness);
        else
            linearly_implicit_euler(q, qdot, dt, M, force, stiffness, tmp_force, tmp_stiffness);
        //std::cout<<"after implicit forces";


    KE = 0;
    PE = 0;
    //std::cout<<"before energy calculations";
    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D);
        PE += V_ele;
    }
    //std::cout<<"after energy calculation";
    Visualize::add_energy(t, KE, PE);
    //std::cout<<"after visualize";

}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    Visualize::update_vertex_positions(0, P.transpose()*q + x0);

}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    //if(key =='N') {
   //     std::cout<<"toggle integrators \n";
    //    fully_implicit = !fully_implicit;
    //} else
    if(key == 'S') {
        skinning_on = !skinning_on;
        Visualize::toggle_skinning(skinning_on);
    } else if(key=='J' && cube){
        plane=!plane;
    }else if(key=='P' && cube){
        point=!point;
    }else if(key=='L' && cube){
        line=!line;
    }
    return false;
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 
    igl::readMESH("../data/coarser_bunny.mesh",V,T, F);
    igl::readOBJ("../data/bunny_skin.obj", V_skin, F_skin);

    if(argc > 1) {
        if(strcmp(argv[1], "arma") == 0) {
            read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");
            igl::readOBJ("../data/armadillo.obj", V_skin, F_skin);
            bunny = false;
            fully_implicit = true;
        }
        if(strcmp(argv[1], "cube") == 0) {
            igl::readMESH("../data/cube86.mesh",V,T, F);
            bunny = false;
            fully_implicit = true;
            V_skin = V;
            F_skin = F;
            cube=true;
        }
        if(strcmp(argv[1], "cube2") == 0) {
            igl::readMESH("../data/cube86.mesh",V,T, F);
            bunny = false;
            fully_implicit = true;
            V_skin = V;
            F_skin = F;
            cube2=true;
        }
    }
    
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();
    build_skinning_matrix(N, V, T, V_skin);
    //setup simulation 
    init_state(q,qdot,V);

    //add geometry to scene
    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::toggle_skinning(false);
    //bunny
    if(bunny)
        Visualize::set_picking_tolerance(1.);
    else
        Visualize::set_picking_tolerance(0.01);

    //volumes of all elements
    igl::volume(V,T, v0);

    //Mass Matrix
    mass_matrix_mesh(M, qdot, T, density, v0);
    
    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }
    
    //setup constraint matrix
    if(bunny)
        find_min_vertices(fixed_point_indices, V, 3);
    else
        find_min_vertices(fixed_point_indices, V, 0.1);
    if(cube){
        fixed_point_indices.clear();
    }
    if(cube2){
        fixed_point_indices.clear();
        fixed_point_indices.push_back(0);
        fixed_point_indices.push_back(1);
        fixed_point_indices.push_back(21);
        fixed_point_indices.push_back(23);
    }

    //material properties
    //bunny
    if(bunny) {
        YM = 6e6; //young's modulus
        mu = 0.49; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e8;
    } else {
        //arma
        YM = 6e5; //young's modulus
        mu = 0.49; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e5;
    }
    if(cube||cube2){
        YM = 6e2; //young's modulus
        mu = 0.49; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e3;
    }
    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixed_point_constraints(P, q.rows(), fixed_point_indices);
    
    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else    
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();

    //igl additional menu setup
    // Add content to the default menu window
    Visualize::viewer_menu().callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * Visualize::viewer_menu().menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
            "Energy Plot", nullptr,
            ImGuiWindowFlags_NoSavedSettings

        );

        ImVec2 min = ImGui::GetWindowContentRegionMin();
        ImVec2 max = ImGui::GetWindowContentRegionMax();

        max.x = ( max.x - min.x ) / 2;
        max.y -= min.y + ImGui::GetItemsLineHeightWithSpacing() * 3;

        Visualize::plot_energy("T", 1, ImVec2(-15,10), ImVec2(0,2e8), ImGui::GetColorU32(ImGuiCol_PlotLines));
        Visualize::plot_energy("V", 2, ImVec2(-15,10), ImVec2(0,2e7), ImGui::GetColorU32(ImGuiCol_HeaderActive));
        Visualize::plot_energy("T+V", 3, ImVec2(-15,10), ImVec2(0,4e8), ImGui::GetColorU32(ImGuiCol_ColumnActive));

        ImGui::End();
    };

    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

