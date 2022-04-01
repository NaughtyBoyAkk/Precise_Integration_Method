    subroutine ComputeExpH(h,ExpH,h_N,tau,div_2_N,taylor_N)
    implicit none
    !hΪָ�����󣬷���
    !h_NΪָ���������״
    !tauΪʱ�䲽��
    !div_2_NΪ2��ָ��
    !taylor_NΪ̩��չ��ʽ�Ľ�ȡ����
    !���ݽ����ı���������
    integer(kind=4)::h_N,taylor_N,div_2_N
    real(kind=8)::tau
    real(kind=8)::h(h_N,h_N),ExpH(h_N,h_N)
    !��������������������
    integer(kind=4)::nfactorial,i
    real(kind=8)::delta_t
    real(kind=8)::hdt(h_N,h_N),Taylor(h_N,h_N,2:taylor_N)
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
    !Taylorչ��ʽ�ļ��㣬����ӵڶ��ʼ����
    do i=2,taylor_N,1
        hdt=matmul(hdt,h*delta_t)
        nfactorial=nfactorial*(i-1)
        ExpH=ExpH+hdt/dble(nfactorial)
    end do
    !��ʼָ����˵Ĺ���
    do i=1,div_2_N,1
        ExpH=2.0d0*ExpH+matmul(ExpH,ExpH)
    end do
    forall(i=1:h_N:1)
        ExpH(i,i)=ExpH(i,i)+1.0d0
    end forall
    end subroutine ComputeExpH
