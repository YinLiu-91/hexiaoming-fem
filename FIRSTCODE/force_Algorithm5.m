function [B]=force_Algorithm5(Pb,Tb,r,s,element_order,N)
%ʹ�ú���ʦ����֮����㷨5���γ��Ҳ��غɾ���
%Pb:����Ԫ�ڵ�������Ϣ����for test function
%Tb:����Ԫ��Ԫ��Ϣ����for test function
%r:��̽�������󵼽״�
%s:���Ժ������󵼽״�
%element_order:��Ԫ�״�
%Pb��Tb��ǰ���information_input����
%KΪ�����ȫ�ָնȾ���
%�ż���
%syms x;
%v=int(-exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x))*x,0,1)
%vpa(v)
%ֱ�Ӳ�������symsx;G=int(-exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x))*(1-x))
%H=int(-exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x))*(x))
%���G =-(exp(x)*(2*cos(x) + sin(x) + 2*x^2*sin(x) - x*cos(x) - 3*x*sin(x)))/2
%���H=(exp(x)*(sin(x) + 2*x^2*sin(x) - x*cos(x) - x*sin(x)))/2
if r==1 && s==1    %��ʱû�������ôд
if element_order==1
    Nb_force=size(Pb,2);%����δ֪�ڵ���߻���������
    B=zeros(Nb_force,1);
    for i=1:N
        T_current=Tb(:,i);
        eta1=-0.577350269;           %��˹���ֵ�1
        eta2=0.577350269;            %��˹���ֵ�2
        eta=[eta1,eta2];
        L=size(eta,2);
        w1=1;                     %��˹���ֵ�Ȩ��
        w2=1;                     %��˹���ֵ�Ȩ��
        w=[w1,w2];
        Nb_test=size(Tb,1); %ÿ����Ԫ�ڵ�����Ļ�����������������Tb��������ȡ
            for n=1: Nb_test
                r=0;%�����ڸ�˹���õ�����ͣ�����Ҫ��r���г�ʼ�����Ӷ��������ڵ���             
             for g=1:L
        x=(Pb(T_current(2))+Pb(T_current(1)))/2+...
            ((Pb(T_current(2))-Pb(T_current(1)))*eta(g))/2; %etaΪ��Ȼ�����
        f=-exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x))*...
            ((Pb(T_current(2))-Pb(T_current(1)))/2)*w(g);%���������Ϊ[-1,1]���Ӷ�����ʹ�ø�˹����
        N_shape1=(Pb(T_current(2))-x)/(Pb(T_current(2))-Pb(T_current(1)));
        N_shape2=(x-Pb(T_current(1)))/(Pb(T_current(2))-Pb(T_current(1)));
        N=[N_shape1,N_shape2];
        r=r+f*N(n);%ע���κ���Ӧ��Ҫ��N_test��Ӧ�����ſ���
              end
            B(Tb(n,i),1)=r+B(Tb(n,i),1); 
            end
    end
end
end 
end