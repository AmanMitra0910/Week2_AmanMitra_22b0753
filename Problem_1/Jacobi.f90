program P1
    implicit none
    real :: d(3 , 3) , lu(3 , 3) , x(1 , 3) , b (1 , 3) , luxx(1 , 3) , s(1 , 3) , di(3 , 3) , y(1 , 3) , a(3 , 3)
    real :: t
    integer :: i , n_max

    d = reshape([4 ,0 ,0 ,0 ,6 ,0 ,0 ,0 ,7], shape(d))
    lu = reshape([0 ,-1 ,-1 ,-2 ,0 ,1 , -1, 1, 0], shape(lu))
    b = reshape([3 ,9 ,-6], shape(b))
    x = reshape([0, 0, 0], shape(x))
    t = 1.0e-5
    n_max = 100
    a = d + lu

    do i = 1 , n_max
        
        luxx = mlt13(lu , x)
        s = b - luxx
        di = inv33(d)

        x = mlt13(di , s)
        y = mlt13(a , x) - b

        if( (abs(y(1,1)) < t) .AND. (abs(y(1,2)) < t) .AND. (abs(y(1,3)) < t)) exit
        
    end do

    print *, "The solution for this set of linear equations using Jacobi Method is,"
    print *, "x =", x(1 , 1)
    print *, "y =", x(1, 2)
    print *, "z =", x(1, 3)

contains

    function mlt13(m1 , m2) result(mlt)
        implicit none
        real, intent(in) :: m1(3 , 3) , m2(1 , 3)
        real :: mlt(1 , 3)
        integer :: i

        do i = 1, 3
            mlt(1, i) = sum(m1(:, i) * m2(1, :))
        end do
    end function mlt13

    function inv33(m) result(inv)
        implicit none
        real, intent(in) :: m(3 , 3)
        real :: inv(3 , 3)
        real :: adj(3 , 3) , detm

        detm = det(m)
        
        adj(1,1) =  m(2,2)*m(3,3) - m(2,3)*m(3,2)
        adj(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
        adj(1,3) =  m(1,2)*m(2,3) - m(1,3)*m(2,2)
        adj(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
        adj(2,2) =  m(1,1)*m(3,3) - m(1,3)*m(3,1)
        adj(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
        adj(3,1) =  m(2,1)*m(3,2) - m(2,2)*m(3,1)
        adj(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
        adj(3,3) =  m(1,1)*m(2,2) - m(1,2)*m(2,1)

        inv = transpose(adj) / detm

    end function inv33

    function det(m) result(deter)
        implicit none
        real, intent(in) :: m(3 , 3)
        real :: deter

        deter = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2))&
               - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1))&
               + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))

    end function det

end program P1