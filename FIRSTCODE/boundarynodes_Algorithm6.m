function [A,B]=boundarynodes_Algorithm6(N,A,B)
%ʹ�ú���ʦ����֮����㷨6������Լ��
%boundarynodes��Լ���ڵ��б�
%1:��ʾ�������ױ߽�����
%2:Ŧ���߽�����
%3���ޱ��߽�����
%�ż���
boundarynodes=[1,1;1,N+1;0,cos(1)];
nbn=size(boundarynodes,2);
for k=1:nbn
if boundarynodes(1,k)==1
    i=boundarynodes(2,k);
    A(i,:)=0;
    A(i,i)=1;
    B(i)=boundarynodes(3,k);
end
end
end