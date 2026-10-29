module plane42rectskewd

    !! This module contains subroutines specific to the plane42 element 
    !!        
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o 
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!             
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc  
    !! `E1` = element face 1, `E2` = element face 2, etc  
    
    use types
    implicit none
    save
    
    private
    public :: plane42rect_ke, plane42rect_re, plane42rect_ss

contains

    subroutine plane42rect_shape(ex, ey, xi, eta, n, B, J_det, J)
        implicit none

        real(wp), intent(in) :: ex(4), ey(4), xi, eta
        real(wp), intent(out) :: B(3, 8), J_det, J(2, 2)

        real(wp) :: N_thilde(4, 8), L(3, 4), Gamma(2, 2), Gamma_thilde(4, 4), n(2, 8)
        real(wp) :: xi_i(4) = [-1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp]
        real(wp) :: eta_i(4) = [-1.0_wp, -1.0_wp, 1.0_wp, 1.0_wp]
        real(wp) :: dNdxi(4), dNdeta(4)
        integer :: i

        L = 0.0_wp
        L(1, 1) = 1
        L(2, 4) = 1
        L(3,2) = 1
        L(3,3) = 1

        n = 0.0_wp
        N_thilde = 0.0_wp
        do i = 1, 4
            dNdxi(i) = 0.25_wp * xi_i(i) * (1 + eta_i(i) * eta)
            dNdeta(i) = 0.25_wp * eta_i(i) * (1 + xi_i(i) * xi)
            N_thilde(1, 2*i-1) = dNdxi(i)
            N_thilde(2, 2*i-1) = dNdeta(i)
            N_thilde(3, 2*i) = dNdxi(i)
            N_thilde(4, 2*i) = dNdeta(i)
            n(1, 2*i-1) = 0.25_wp * (1+xi_i(i)*xi)*(1+eta_i(i)*eta)
            n(2, 2*i) = 0.25_wp * (1+xi_i(i)*xi)*(1+eta_i(i)*eta)
        end do
        
        J(1, 1) = sum(dNdxi * ex)
        J(2, 1) = sum(dNdeta * ex)
        J(1, 2) = sum(dNdxi * ey)
        J(2, 2) = sum(dNdeta * ey)

        J_det = J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)
        Gamma(1, 1) = J(2, 2) / J_det
        Gamma(2, 1) = -J(2, 1) / J_det
        Gamma(1, 2) = -J(1, 2) / J_det
        Gamma(2, 2) = J(1, 1) / J_det

        Gamma_thilde = 0.0_wp
        Gamma_thilde(1:2, 1:2) = Gamma
        Gamma_thilde(3:4, 3:4) = Gamma

        B = matmul(L, matmul(Gamma_thilde, N_thilde))

    end subroutine plane42rect_shape
!
!--------------------------------------------------------------------------------------------------
!

    subroutine plane42rect_ke(xe, young, nu, thk, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.
        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
        real(wp) :: cmat(3,3), fact
        real(wp), parameter :: gauss_points(2) = [-1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp)]
        real(wp), parameter :: gauss_weight(2) = [1.0_wp, 1.0_wp]
        real(wp) :: ex(4), ey(4)
        real(wp) :: B(3, 8), J_det, J(2, 2), n(2, 8)
        integer :: i, y, n_g

        ! Resetting element stiffness matrix
        ke = 0.0_wp

        ! build constitutive matrix (plane stress)
        cmat = 0.0_wp
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        ! Gaus-points, 2*n_g - 1 = p, p = 2 => n_g >= 1.5 => n_g = 2
        n_g = 2
        ex = xe(1::2)
        ey = xe(2::2)
        
        do i = 1, n_g
            do y = 1, n_g
                call plane42rect_shape(ex, ey, gauss_points(i), gauss_points(y), n, B, J_det, J)
                ke = ke + gauss_weight(i)*gauss_weight(y)*thk*matmul(transpose(B), matmul(cmat, B))*J_det
            end do
        end do

    end subroutine plane42rect_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42rect_re(xe, eface, fe, thk, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp), parameter :: gauss_weight(2) = [1.0_wp, 1.0_wp]
        real(wp) :: ex(4), ey(4)
        real(wp) :: B(3, 8), J_det, J(2, 2), J_cal(2), n(2, 8)
        real(wp), parameter :: gauss_points(2) = [-1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp)]
        integer :: i
        integer, parameter :: n_g = 2
        
        ex = xe(1::2)
        ey = xe(2::2)
        J_cal = 0.0_wp
        re = 0.0_wp

        if (eface == 1) then
            do i = 1, n_g
                call plane42rect_shape(ex, ey, gauss_points(i), -1.0_wp, n, B, J_det, J)
                J_cal(1) = -J(1, 2)
                J_cal(2) = J(1, 1)
                re = re + gauss_weight(i) * thk * fe * matmul(transpose(n), J_cal)
            end do
        elseif (eface == 2) then
            do i = 1, n_g
                call plane42rect_shape(ex, ey, 1.0_wp, gauss_points(i), n, B, J_det, J)
                J_cal(1) = -J(2, 2)
                J_cal(2) = J(2, 1)
                re = re + gauss_weight(i) * thk * fe * matmul(transpose(n), J_cal)
            end do
        elseif (eface == 3) then
            do i = 1, n_g
                call plane42rect_shape(ex, ey, gauss_points(i), 1.0_wp, n, B, J_det, J)
                J_cal(1) = J(1, 2)
                J_cal(2) = -J(1, 1)
                re = re + gauss_weight(i) * thk * fe * matmul(transpose(n), J_cal)
            end do
        elseif (eface == 4) then
            do i = 1, n_g
                call plane42rect_shape(ex, ey, -1.0_wp, gauss_points(i), n, B, J_det, J)
                J_cal(1) = J(2, 2)
                J_cal(2) = -J(2, 1)
                re = re + gauss_weight(i) * thk * fe * matmul(transpose(n), J_cal)
            end do
        endif
    end subroutine plane42rect_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42rect_ss(xe, de, young, nu, estress, estrain)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu 
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)

        ! New variables
        real(wp) :: bmat(3, 8), cmat(3, 3) 
        real(wp) :: aa, bb, dndx(4), dndy(4), fact
        integer :: i

        ! Build strain-displacement matrix

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2
        dndx = [-1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp]/(4*aa)
        dndy = [-1.0_wp, -1.0_wp, 1.0_wp, 1.0_wp]/(4*bb)

        bmat = 0
        do i = 1, 4
            bmat(1,2*i-1) = dndx(i)
            bmat(2,2*i) = dndy(i)
            bmat(3,2*i-1) = dndy(i)
            bmat(3,2*i) = dndx(i)
        end do

        ! Compute element strain
        estrain = matmul(bmat, de)

        ! Build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1.0_wp-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1.0_wp-nu)/2

        ! Compute element stress
        estress = matmul(cmat, estrain)

        ! Compute principal stress and direction
        ! Does this in fea.f90 to not change the subroutine

    end subroutine plane42rect_ss

end module plane42rectskewd
