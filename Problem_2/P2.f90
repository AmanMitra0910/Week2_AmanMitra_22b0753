program poisson_solver
    implicit none

    integer, parameter :: n = 50
    integer :: i, j, k
    real(8), dimension(n, n) :: cd_grid, u
    real(8), dimension(n*n) :: rho_flat, u_flat
    real(8), dimension(n*n, n*n) :: opr
    real(8), parameter :: h = 1.0 / (n - 1)
    real(8), parameter :: tol = 1.0e-10
    integer, parameter :: n_max = 1000

    cd_grid = 0.0
    u = 0.0

    cd_grid(n/2 + 10, n/2) = 1.0
    cd_grid(n/2 - 10, n/2) = -1.0

    u(1, :) = 1.0
    u(n, :) = -1.0

    call grid2flat(cd_grid, rho_flat)
    call grid2flat(u, u_flat)

    call make_laplacian_operator(opr)

    do i = 1, n
        do j = 1, n
            k = (i - 1) * n + j
            if (i == 1 .or. i == n .or. j == 1 .or. j == n) then
                rho_flat(k) = u_flat(k)
                opr(k, :) = 0.0
                opr(k, k) = 1.0
            else
                rho_flat(k) = -rho_flat(k) * h**2
            end if
        end do
    end do

    call solve_system_lineq(opr, rho_flat, u_flat, n_max, tol)

    call flat2grid(u_flat, u)

    open(unit=10, file='potential_distribution.txt')
    do i = 1, n
        write(10, *) u(i, :)
    end do
    close(10)

contains

    subroutine grid2flat(grid, flat)
        implicit none
        integer, parameter :: n = 50
        real(8), dimension(n, n), intent(in) :: grid
        real(8), dimension(n*n), intent(out) :: flat
        integer :: i, j, k

        k = 0
        do i = 1, n
            do j = 1, n
                k = k + 1
                flat(k) = grid(i, j)
            end do
        end do
    end subroutine grid2flat

    subroutine flat2grid(flat, grid)
        implicit none
        integer, parameter :: n = 50
        real(8), dimension(n*n), intent(in) :: flat
        real(8), dimension(n, n), intent(out) :: grid
        integer :: i, j, k

        k = 0
        do i = 1, n
            do j = 1, n
                k = k + 1
                grid(i, j) = flat(k)
            end do
        end do
    end subroutine flat2grid

    subroutine make_laplacian_operator(operator)
        implicit none
        integer, parameter :: n = 50
        real(8), dimension(n*n, n*n), intent(out) :: operator
        integer :: i, j, k, l

        operator = 0.0

        do i = 1, n
            do j = 1, n
                k = (i - 1) * n + j
                if (i > 1) then
                    l = (i - 2) * n + j
                    operator(k, l) = 1.0
                end if
                if (i < n) then
                    l = i * n + j
                    operator(k, l) = 1.0
                end if
                if (j > 1) then
                    l = (i - 1) * n + (j - 1)
                    operator(k, l) = 1.0
                end if
                if (j < n) then
                    l = (i - 1) * n + (j + 1)
                    operator(k, l) = 1.0
                end if
                operator(k, k) = -4.0
            end do
        end do
    end subroutine make_laplacian_operator

    subroutine solve_system_lineq(opr, rhs, sol, n_max, tol)
        implicit none
        integer, parameter :: n = 50
        integer, parameter :: n2 = n*n
        integer, intent(in) :: n_max
        real(8), intent(in) :: tol
        real(8), dimension(n2, n2), intent(in) :: opr
        real(8), dimension(n2), intent(in) :: rhs
        real(8), dimension(n2), intent(out) :: sol
        real(8) :: sum
        integer :: iter, i, j

        real(8), dimension(n2) :: old_sol

        sol = 0.0

        do iter = 1, n_max
            old_sol = sol
            do i = 1, n2
                sum = 0.0
                do j = 1, n2
                    if (i /= j) then
                        sum = sum + opr(i, j) * old_sol(j)
                    end if
                end do
                sol(i) = (rhs(i) - sum) / opr(i, i)
            end do
            if (maxval(abs(sol - old_sol)) < tol) exit
        end do
    end subroutine solve_system_lineq

end program poisson_solver

