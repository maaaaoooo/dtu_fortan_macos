program main

    !! The "main" routine of the program
    !!
    !! This is where execution starts when running the program
    
    use fedata
    use processor
    use fea

    implicit none

    ! Read model data
    call input

    ! Initialize problem
    call initial

    ! Split in dynamic and static
    select case (trim(antype))
    ! If static -> displacement and compliance calculated
    case ('STATIC', 'static')
        call displ
        call compliance
    case ('MODAL', 'modal')
    ! If dynamic -> the eigenvalues are calculated
        call eigen
    end select

end program main
