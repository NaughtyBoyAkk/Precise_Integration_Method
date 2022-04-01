    subroutine ComputeExpH2(h,ExpH,h_N,tau,div_2_N,taylor_N)
    implicit none
    !h为指数矩阵，方阵
    !h_N为指数矩阵的形状
    !tau为时间步长
    !div_2_N为2的指数
    !taylor_N为泰勒展开式的截取项数
    !传递进来的变量的声明
    integer(kind=4)h_N,taylor_N,div_2_N
    real(kind=8)tau
    real(kind=8)h(h_N,h_N),ExpH(h_N,h_N)
    !本程序其他变量的声明
    integer(kind=4)nfactorial,i,j,k,m
    real(kind=8)delta_t
    real(kind=8)hdt(h_N,h_N),taylor(h_N,h_N,taylor_N),taylor_temp(h_N,h_N,taylor_N)

    ExpH=0.0D0
    nfactorial=1
    delta_t=tau
    do i=1,div_2_N,1
        delta_t=delta_t/2.0d0
    end do

    !初始化hdt矩阵,该矩阵存储 (H*delta_t)**taylor_N
    hdt=0.0d0
    forall(i=1:h_N:1)
        hdt(i,i)=1.0d0
    end forall
    Taylor(:,:,1)=hdt(:,:)
    !Taylor展开式的计算，这里从第二项开始计算,为极大可能的保留精度
    do i=2,taylor_N,1
        hdt=matmul(hdt,h*delta_t)
        nfactorial=nfactorial*(i-1)
        Taylor(:,:,i)=hdt/dble(nfactorial)
    end do
    !根据阶次来划分的，同阶小量放在一个变量，大于等于taylor_N阶次的小量统一作为最后一个变量
    !这个操作用一个表格来表示比较直观，以副对角线为分界线
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
    !级数的求和
    do i=1,taylor_N,1
        ExpH=ExpH+taylor(:,:,i)
    end do
    end subroutine ComputeExpH2

