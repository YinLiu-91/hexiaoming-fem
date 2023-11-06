function result=poisson_solver_1D(left,right,h_partition,basis_type,Gauss_point_number)
%Xiaoming He, 10/08/2011.
%Finite element solver for 1D Poisson equation.
%We will use "FE" to replace "finite element" in the comments.
%The trial FE functions and test FE functions need to be the same.
%So far the code can only handle 1D linear FE and 1D Lagrange quadratic FE.
%For other types of FE, we need to add the corresponding information to code.
%The problem domain is [left,right].
%h_partition is the step size of the partition.
%basis_type:the type of the FE basis function.
%basis_type=101:1D linear FE.
%basis_type=102:1D quadratic FE.

%N_basis:The N for the FE basis functions,not the partition.
%N_partition:The N for the partition,not the FE basis functions.
%N is the number of the sub-intervals.
%M_partition,T_partition, M_basis,T_basis: see the note in "generate_M_T_1D.m".
%function_a: the coefficient function on the left side of the poisson equation.
%funciton_f: the right hand side function of the poisson equation.
%function_g: the Dirichelet boundary function in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_q_tilde: the name of the Neumann boundary function q(x,y) when p(x,y)=0 in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_q: the name of the Neumann boundary function q(x,y) when p(x,y) is nonzero in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_p: the name of the Robin coefficient function p(x,y) in my notes "Notes for tool box of standard triangular FE" section 1-1.

%h_basis is the step size for the finite element nodes, not the partition.

% 单元数
N_partition=(right-left)/h_partition;

% 根据单元类型，得到基函数个数
if basis_type==102
    N_basis=N_partition*2;
elseif basis_type==101
    N_basis=N_partition;
end

%Mesh information for partition and finite element basis functions.
[M_partition,T_partition]=generate_M_T_1D(left,right,h_partition,101);

if basis_type==102
    [M_basis,T_basis]=generate_M_T_1D(left,right,h_partition,102);
elseif basis_type==101
    M_basis=M_partition;
    T_basis=T_partition;
end 


%Guass quadrature's points and weights on the refenrece interval [-1,1].
% 产生1D 单元的高斯点坐标与系数
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number);

%Assemble the stiffness matrix.
number_of_elements=N_partition; % 注意单元数与N_basis 的区别
matrix_size=[N_basis+1 N_basis+1]; % N_basis是有限元节点的个数
vector_size=N_basis+1;
if basis_type==102
    number_of_trial_local_basis=3; % 2阶的话每个单元三个有限元节点，所以 trial和test 都是3
    number_of_test_local_basis=3;
elseif basis_type==101
    number_of_trial_local_basis=2; % 同上，但是1阶只有两个有限元节点，trial 和test都是2
    number_of_test_local_basis=2;
end

%Assemble the stiffness matrix 组装刚度矩阵，因为是椭圆方程，所以trial和test都是一阶导数，二次元的一阶导数
A=assemble_matrix_from_1D_integral('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,1,basis_type,1);


%Assemble the load vector.

b=assemble_vector_from_1D_integral('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0);

%Get the information matrix for boundary nodes.
boundary_nodes=generate_boundary_nodes_1D(N_basis);

%Deal with Dirichlet boundary condition. 注意，dirchlet边界要在最后处理，即robin，neumann先于它处理
[A,b]=treat_Dirichlet_boundary_1D('function_g',A,b,boundary_nodes,M_basis);

%Compute the numerical solution
result=A\b;

%compute the maximum error at all nodes.
if basis_type==102
    h_basis=h_partition/2;
elseif basis_type==101
    h_basis=h_partition;
end
maxerror=get_maximum_error_1D(result,N_basis,left,h_basis);
maximum_error_at_all_nodes_of_FE=maxerror







