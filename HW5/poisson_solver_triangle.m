function r=poisson_solver_triangle(left,right,bottom,top,h_partition,basis_type)
%Xiaoming He, 07/04/2009.
%Finite element solver for Poisson equation on a triangular mesh.
%We will use "FE" to replace "finite element" in the comments.
%The trial FE functions and test FE functions need to be the same.
%So far the code can only handle 2D linear FE and 2D Lagrange quadratic FE.
%For other types of FE, we need to add the corresponding information to code.
%The problem domain is [left,right]*[bottom,top].
%h_partition is the step size of the partition.
%basis_type:the type of the FE basis function.
%basis_type=1:2D linear FE.  
%basis_type=2:2D Lagrange quadratic FE.

%N1_basis,N2_basis:The N1 and N2 for the FE basis functions,not the partition.
%N1_partition,N2_partition:The N1 and N2 for the partition,not the FE basis functions.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%M_partition,T_partition, M_basis,T_basis: see the note in "generate_M_T_triangular.m".
%function_a: the coefficient function on the left side of the poisson equation.
%funciton_f: the right hand side function of the poisson equation.
%function_g: the Dirichelet boundary function in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_q_tilde: the name of the Neumann boundary function q(x,y) when p(x,y)=0 in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_q: the name of the Neumann boundary function q(x,y) when p(x,y) is nonzero in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_p: the name of the Robin coefficient function p(x,y) in my notes "Notes for tool box of standard triangular FE" section 1-1.

%h_basis is the step size for the finite element nodes, not the partition.

N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);

if basis_type==2
    N1_basis=N1_partition*2;
    N2_basis=N2_partition*2;
elseif basis_type==1
    N1_basis=N1_partition;
    N2_basis=N2_partition;
end

%Mesh information for partition and finite element basis functions.
[M_partition,T_partition]=generate_M_T_triangle(left,right,bottom,top,h_partition,1);

if basis_type==2 %对于二阶型函数，需要产生有限元新的节点信息，一阶的话则与网格节点重合
    [M_basis,T_basis]=generate_M_T_triangle(left,right,bottom,top,h_partition,2);
elseif basis_type==1
    M_basis=M_partition;
    T_basis=T_partition;
end 

% 产生标准单元的积分点坐标与权重信息
%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);


%Assemble the stiffness matrix.
number_of_elements=2*N1_partition*N2_partition;
matrix_size=[(N1_basis+1)*(N2_basis+1) (N1_basis+1)*(N2_basis+1)];    % 注意，这里的矩阵大小是有限元节点的个数
if basis_type==2    % 根据阶次不同，确定每个单元不同的trial/test函数个数
    number_of_trial_local_basis=6;
    number_of_test_local_basis=6;
elseif basis_type==1
    number_of_trial_local_basis=3;
    number_of_test_local_basis=3;
end   % 下面这里由于有了分别对x，y坐标的基函数的积分，所以生成了两个矩阵，然后进行矩阵的对位相加即可得到体积分的矩阵结果
A1=assemble_matrix_from_volume_integral_triangle('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,1,0,basis_type,1,0); % 这里表示trial，和test都是对x求导
A2=assemble_matrix_from_volume_integral_triangle('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,1,basis_type,0,1); % 这里表示trial，test都是对y求导
A=A1+A2;

%Assemble the load vector. 注意load 是加在每个单元上的，所以还是要对单元进行面积分，只不过只有test 型函数了，没有trial了，加在方程的组的左端，为向量
vector_size=(N1_basis+1)*(N2_basis+1);
b=assemble_vector_from_volume_integral_triangle('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,0);

%Get the information matrices for boundary nodes and boundary edges. 
% boundary_nodes 记录了每个边界点是属于什么边界类型，
% boundary_edges记录了1. 第一行为所有边的边界条件类型，第二行为单元编号，第三行为节点a编号，第四行为节点b编号
[boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1_basis,N2_basis,N1_partition,N2_partition);

%Deal with Dirichlet boundary condition.
[A,b]=treat_Dirichlet_boundary_triangle('function_g',A,b,boundary_nodes,M_basis);

%Compute the numerical solution
r=A\b;

%Transfer the 1D solution into 2D solution and compute the maximum error at all nodes.
if basis_type==2
    h_basis=h_partition/2;
elseif basis_type==1
    h_basis=h_partition;
end
[solution_2D,maxerror]=get_2D_solution_and_maximum_error(r,N1_basis,N2_basis,left,bottom,h_basis);
maximum_error_at_all_nodes_of_FE=maxerror





