function [A,b]=treat_Dirichlet_boundary_time_theta_triangle(Dirichlet_boundary_function_name,current_time,dt,theta,A,b,boundary_nodes,M_basis)
%Xiaoming He, 07/17/2009.
%Deal with Dirichlet boundary nodes.
%We will use "FE" to replace "finite element" in the comments.
%boundary_nodes(1,k): specifiy the type of the kth boundary node.
%boundary_nodes(1,k)=-1: Dirichlet boundary node;
%boundary_nodes(1,k)=-2: Neumann boundary node;
%boundary_nodes(1,k)=-3: Robin boundary node. 
%boundary_nodes(2,k): global index of the kth boundary node among all nodes of FE. 
%                     That is, the index of FE is used here.
%M_basis: store the coordinates of all the nodes for the FE,not the partition. 
%See "generate_M_T_triangular.m" for more explanation about M_basis.
%Dirichelet_boundary_fucntion_name: the name of the Dirichelet boundary function.
%current_time: the time when we want to evaluate the Dirichlet boundary function.
%dt is the step size of the time partition.
%theta: decide the time discretization algorithm.
%0<=theta<=1.
%theta=0: forward Euler.
%theta=1: backward Euler.
%theta=0.5: Crank-Nicolson.
%More explanation is in my "Notes for tool box of standard triangular FE" section 3-2-1.

%nbn: the total number of all the boundary nodes of FE.

nbn=size(boundary_nodes,2);

%Check all boundary nodes of FE.
for k=1:nbn

%If the ith FE node X_i is a Dirichlet boundary node,then we reset the ith equation in the linear sysytem to be 1*u_i=u(X_i).
%Here u_i is the unknown associated with the node X_i and u(X_i)=Dirichelet_boundry_function_name(X_i).
%Also i=boundary_nodes(2,k) if boundary_nodes(1,k)==-1.

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=theta*feval(Dirichlet_boundary_function_name,current_time,M_basis(1,i),M_basis(2,i))+(1-theta)*feval(Dirichlet_boundary_function_name,current_time-dt,M_basis(1,i),M_basis(2,i));
    end

end