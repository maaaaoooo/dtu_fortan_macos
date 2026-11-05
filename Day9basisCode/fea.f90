module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, compliance

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial

        !! This subroutine is mainly used to allocate vectors and matrices

        use fedata
        use link1
        use plane42rectskewd

! Hint for continuum elements:
        !integer, parameter :: mdim = 8
        integer, dimension(mdim) :: edof
        integer :: e, i, nen

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            ! Making sure bw = 0
            bw = 0
            ! Looping through the different elements
            do e = 1, ne
            nen = element(e)%numnode
                ! Picking out the different edofs
                do i = 1, nen
                    edof(2*i-1) = 2 * element(e)%ix(i) - 1
                    edof(2*i)   = 2 * element(e)%ix(i)
                end do
                ! Updating bw if the bandwith of the element is bigger than the previous ones
                bw = max(bw, (maxval(edof(1:2*nen)) - minval(edof(1:2*nen)) + 1))
            end do
            allocate (kmat(bw, neqn))
        end if
        allocate (p(neqn), d(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
        allocate (stress_eq(ne, 4))

        ! Initial stress and strain
        strain = 0
        stress = 0
        stress_eq = 0
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e
        real(wp), dimension(:), allocatable :: plotval

        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
            ! Factor banded stiffness matrix
            call bfactor(kmat)
            ! Solve for displacement vector
            call bsolve(kmat, p)
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)

        ! Recover stress
        call recover

        ! Output results
        call output

        ! Plot deformed structure
        call plotmatlabdef('Deformed')

        ! Plot element values
        allocate (plotval(ne))
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                plotval(e) = stress(e,3)
            end if
        end do
        call plotmatlabeval('Shear stress',plotval)
        call plotmatlabeval('von Mises', stress_eq(:, 1))
        call plotmatlabevec('Principal stresses', stress_eq(:, 2), stress_eq(:, 3), stress_eq(:, 4))

    end subroutine displ

!
!--------------------------------------------------------------------------------------------------
!
    subroutine compliance

        !! This subroutine calculates compliance

        use fedata

        real(wp) :: c

        c = dot_product(p, d)
        write (*, '("Compliance C = p^T d = ", es15.8)') c        
        open (20, file='compliance.dat', position='append')
        write (20, '(i8, 1x, i8, 1x, es20.10)') ne, neqn, c
        close (20)

    end subroutine compliance


!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42rectskewd
        
! Hint for continuum elements:

        integer :: i, j, e, nen, idof, eface
        real(wp) :: fe, thk
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re

        ! Build load vector
        p(1:neqn) = 0
        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
                idof = int(2*(loads(i,2)-1) + loads(i,3))
                p(idof) = p(idof) + loads(i, 4)
            case( 2 )
            	! Build uniformly distributed surface (pressure) load contribution
                e = int(loads(i, 2))
                eface = int(loads(i,3))
                fe    = loads(i,4)
                thk   = mprop(element(e)%mat)%thk
                nen   = element(e)%numnode

                do j = 1, nen
                    xe(2*j-1) = x(element(e)%ix(j),1)
                    xe(2*j  ) = x(element(e)%ix(j),2)
                    edof(2*j-1) = 2 * element(e)%ix(j) - 1
                    edof(2*j)   = 2 * element(e)%ix(j)
                end do
                call plane42rect_re(xe, eface, fe, thk, re)
                p(edof(1: 2*nen)) = p(edof(1: 2*nen)) + re(1:2*nen)
            case default
                print *, 'ERROR in fea/buildload 2'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rectskewd

        integer :: e, i, j
        integer :: nen
! Hint for system matrix in band form:
        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
! Hint for modal analysis:
        real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, area
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        kmat = 0
        !if (.not. banded) then
        !    kmat = 0
        !else
        !    kmat = 0
        !end if

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 nu = mprop(element(e)%mat)%nu
                 thk = mprop(element(e)%mat)%thk
                 call plane42rect_ke(xe, young, nu, thk, ke)
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do
! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
            else
                do i = 1, 2*nen
                    do j = 1, 2*nen
                    if (edof(j) >= edof(i)) then
                        kmat(edof(j)-edof(i)+1, edof(i)) = kmat(edof(j)-edof(i)+1, edof(i)) + ke(i, j)
                    end if
                    end do
                end do
            end if  
        end do
    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof, k, r
        real(wp) :: penal

        ! Correct for supports
        if (.not. banded) then
                if (.not. penalty) then
                    do i = 1, nb
                        idof = int(2*(bound(i,1)-1) + bound(i,2))
                        p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                        p(idof) = bound(i, 3)
                        kmat(1:neqn, idof) = 0
                        kmat(idof, 1:neqn) = 0
                        kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))

                    do k = (idof-min(idof-1, bw-1)), (idof-1)
                        P(k) = P(k) - kmat(idof+1-k, k) * bound(i, 3)
                        kmat(idof+1-k, k) = 0
                    end do
                    
                    do r = 2, min(bw, neqn - idof + 1)
                        P(idof+r-1) = P(idof+r-1) - kmat(r, idof) * bound(i, 3)
                        kmat(r, idof) = 0
                    end do

                    p(idof) = bound(i, 3)
                    kmat(1, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(1, idof) = kmat(1, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        end if
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain,
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42rectskewd

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area
! Hint for continuum elements:
        real(wp):: nu, dens, thk, sigma_vm, sigma_1, sigma_2, psi
        real(wp), dimension(3) :: estrain, estress

        ! Reset force vector
        p = 0

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
                case( 2 )
                    young = mprop(element(e)%mat)%young
                    nu = mprop(element(e)%mat)%nu
                    thk = mprop(element(e)%mat)%thk
                    call plane42rect_ke(xe, young, nu, thk, ke)
                    p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                    call plane42rect_ss(xe, de, young, nu, estress, estrain)
                    stress(e, 1:3) = estress
                    strain(e, 1:3) = estrain

                    sigma_vm = sqrt((estress(1)**2) + (estress(2)**2) - (estress(1)*estress(2)) + 3*(estress(3)**2))
                    sigma_1 = 0.5 * (estress(1) + estress(2))+sqrt(((estress(1)-estress(2))/2)**2 + (estress(3)**2))
                    sigma_2 = 0.5 * (estress(1) + estress(2))-sqrt(((estress(1)-estress(2))/2)**2 + (estress(3)**2))
                    psi = atan2(((-2*estress(3))), ((estress(1)-estress(2))))/2.0_wp

                    stress_eq(e, 1) = sigma_vm
                    stress_eq(e, 2) = sigma_1
                    stress_eq(e, 3) = sigma_2
                    stress_eq(e, 4) = psi
            end select
        end do
    end subroutine recover

end module fea
