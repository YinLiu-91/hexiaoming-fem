function b=treat_Robin_boundary_for_vector_time_triangle(Neumann_boundary_function_name,time_n_plus_1,time_n,theta,b,boundary_edges,M_partition,T_partition,T_basis_test,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,test_basis_type,test_derivative_degree_x,test_derivative_degree_y)
%Xiaoming He, 07/23/2009.
%Deal with Robin boundary edges for the vector, which depends on time.
%We will use "FE" to replace "finite element" in the comments.
%Neumann_boundary_fucntion_name: the name of the Neumann boundary function q(x,y) in my notes "Notes for tool box of standard triangular FE" section 3-1-6.
%time_n_plus_1,time_n: two neighboring time levels in the time iteration.
%theta: decide the time discretization algorithm.
%0<=theta<=1.
%theta=0: forward Euler.
%theta=1: backward Euler.
%theta=0.5: Crank-Nicolson.%b: the vector affected by the Robin boundary condition.
%boundary_edges(1,k): specifiy the type of the kth boundary edge.
%boundary_edges(1,k)=-1: Dirichlet boundary edge;
%boundary_edges(1,k)=-2: Neumann boundary edge;
%boundary_edges(1,k)=-3: Robin boundary edge. 
%boundary_edges(2,k): index of the element which contains the kth boundary edge.
%boundary_edges(3:4,k): indices of the two end points of the kth boundary edge among all grid points, not the nodes of FE.
%                       That is, the index of partition is used here.
%M_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%%T_basis_test: T_basis for the test basis function.
%The explanation for M_partition,T_partition,T_basis is in generate_M_T_triangular.m.
%%number_of_test_local_basis: the number of local FE basis functions for the test function in a local element.
%Gauss_coefficient_reference_1D,Gauss_point_reference_1D:the Gauss coefficients and Gauss points on the reference interval [-1,1].
%%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=2:2D Lagrange quadratic FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%More explanation is in my "Notes for tool box of standard triangular FE" section 3-2-3.

%nbe: the total number of all the boundary edges.

nbe=size(boundary_edges,2);

%Check all boundary edges.
for k=1:nbe

%If the kth boundary edge is a Robin boundary edge,then we add the corresponding line integral to A, b.

    if boundary_edges(1,k)==-3 
        
        element_index=boundary_edges(2,k);
        vertices=M_partition(:,T_partition(:,element_index));
        end_point_1=M_partition(:,boundary_edges(3,k));
        end_point_2=M_partition(:,boundary_edges(4,k));
        for beta=1:number_of_test_local_basis
            temp1=Gauss_quadrature_for_line_integral_test_time_triangle(Neumann_boundary_function_name,time_n_plus_1,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,vertices,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
            temp2=Gauss_quadrature_for_line_integral_test_time_triangle(Neumann_boundary_function_name,time_n,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,vertices,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
            b(T_basis_test(beta,element_index),1)=b(T_basis_test(beta,element_index),1)+theta*temp1+(1-theta)*temp2;
        end        
    end

end