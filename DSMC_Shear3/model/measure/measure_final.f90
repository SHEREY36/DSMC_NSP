  subroutine measure_final()	
    use output
    use run_param
    use particle
    use constants
    use geometry
    use DSMC
	
    implicit none

    integer :: I,J
    double precision :: vsq(3), Tg
    double precision :: VMEAN, chi, CF0
    double precision :: kIso, NMean

    ! Close files to be replaced in case a re-run
    ! Should seperate file opening from intialization to avoid repeated file
    ! opening and closing
    close(1001); close(1002); close(1003); close(1005); close(1006); close(1007); close(1008);
    close(1011); close(1012); close(1013)

    RETURN
  end subroutine measure_final