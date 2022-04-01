    subroutine ComputeExpH2(h,ExpH,h_N,tau,div_2_N,taylor_N)
    implicit none
    !hΪָ�����󣬷���
    !h_NΪָ���������״
    !tauΪʱ�䲽��
    !div_2_NΪ2��ָ��
    !taylor_NΪ̩��չ��ʽ�Ľ�ȡ����
    !���ݽ����ı���������
    integer(kind=4)h_N,taylor_N,div_2_N
    real(kind=8)tau
    real(kind=8)h(h_N,h_N),ExpH(h_N,h_N)
    !��������������������
    integer(kind=4)nfactorial,i,j,k,m
    real(kind=8)delta_t
    real(kind=8)hdt(h_N,h_N),taylor(h_N,h_N,taylor_N),taylor_temp(h_N,h_N,taylor_N)

    ExpH=0.0D0
    nfactorial=1
    delta_t=tau
    do i=1,div_2_N,1
        delta_t=delta_t/2.0d0
    end do

    !��ʼ��hdt����,�þ���洢 (H*delta_t)**taylor_N
    hdt=0.0d0
    forall(i=1:h_N:1)
        hdt(i,i)=1.0d0
    end forall
    Taylor(:,:,1)=hdt(:,:)
    !Taylorչ��ʽ�ļ��㣬����ӵڶ��ʼ����,Ϊ������ܵı�������
    do i=2,taylor_N,1
        hdt=matmul(hdt,h*delta_t)
        nfactorial=nfactorial*(i-1)
        Taylor(:,:,i)=hdt/dble(nfactorial)
    end do
    !���ݽ״������ֵģ�ͬ��С������һ�����������ڵ���taylor_N�״ε�С��ͳһ��Ϊ���һ������
    !���������һ���������ʾ�Ƚ�ֱ�ۣ��Ը��Խ���Ϊ�ֽ���
    do m=1,div_2_N,1
        taylor_temp=0.0d0
        do i=taylor_N,1,-1
            do j=taylor_N,1,-1
                k=i+j-1
                if (k>=taylor_N) then
                    taylor_temp(:,:,taylor_N)=taylor_temp(:,:,taylor_N)+matmul(taylor(:,:,i),taylor(:,:,j))
                else
                    taylor_temp(:,:,k)=taylor_temp(:,:,k)+matmul(taylor(:,:,i),taylor(:,:,j))
                end if
            end do
        end do
        taylor=taylor_temp
    end do
    !���������
    do i=1,taylor_N,1
        ExpH=ExpH+taylor(:,:,i)
    end do
    end subroutine ComputeExpH2

