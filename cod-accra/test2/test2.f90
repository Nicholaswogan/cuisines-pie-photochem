program main
  use clima, only: AdiabatClimate, dp, s_str_len, Radtran
  use clima_const, only: k_boltz
  implicit none

  call test()

contains

  subroutine test()
    character(:), allocatable :: err
    type(AdiabatClimate) :: c
    real(dp), parameter :: surface_emissivity = 0.9_dp
    real(dp) :: T_surf
    real(dp), allocatable :: P_i_surf(:)
    logical :: converged

    c = AdiabatClimate( &
      "../test2/species.yaml", &
      "../test2/settings.yaml", &
      "../test1/Sun_now.txt", &
      "/Users/nicholas/Documents/Research_local/PhotochemPy/photochem_clima_data", &
      err=err &
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    c%rad%surface_emissivity = surface_emissivity
    c%P_top = 1.0e2_dp

    block
      real(dp) :: stellar_radiation
      integer :: i
      stellar_radiation = 0.0_dp
      do i = 1,c%rad%sol%nw
        stellar_radiation = stellar_radiation + c%rad%photons_sol(i)*(c%rad%sol%freq(i) - c%rad%sol%freq(i+1))
      enddo
      stellar_radiation = stellar_radiation/1.0e3_dp

      print*,1361.0_dp/stellar_radiation

    endblock

    if (size(c%species_names) /= 3) then
      stop 1
    endif
    if (any(c%species_names /= ['H2O','CO2','N2 '])) then
      stop 1
    endif

    P_i_surf = [260.0_dp, 0.393201E-03_dp, 0.983002E+00_dp]
    P_i_surf = P_i_surf*1.0e6_dp

    T_surf = c%surface_temperature(P_i_surf, 280.0_dp, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    converged = c%RCE(P_i_surf, c%T_surf, c%T, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    call write_output_file(c)

  end subroutine

  subroutine write_output_file(c)
    type(AdiabatClimate), target, intent(in) :: c
    type(Radtran), pointer :: rad
    real(dp) :: OLR, ASR
    character(1000) :: tmp
    character(s_str_len), allocatable :: labels(:)
    character(:), allocatable :: fmt_label, fmt_number, fmt_number1
    integer :: io, i
    real(dp) :: mmw

    rad => c%rad

    fmt_label = "(a14)"
    fmt_number = "(es14.4e3)"
    fmt_number1 = "(f14.4)"

    OLR = rad%wrk_ir%fup_n(rad%nz+1)/1.0e3_dp
    ASR = rad%wrk_sol%fdn_n(rad%nz+1)/1.0e3_dp - rad%wrk_sol%fup_n(rad%nz+1)/1.0e3_dp

    open(unit=2,file='../test2/output_ModernEarth_codaccra_photochemclima_test2.txt',status='replace',iostat=io)
    write(2,'(a)') "# Model name = photochem.clima"
    write(2,'(a)') "# Diameter: 12742                # Diameter of the planet [km]"
    write(2,'(a)') "# Gravity: 9.81                  # Gravity at the surface of the atmosphere [m/s2]"
    write(2,'(a)') "# Star-distance: 1               # Planet-star distance [ua]"
    write(2,'(a)') "# Star-type: Sun                 # Host-star type"
    write(2,'(a)') "# Star-temperature: 5772         # Host-star temperature [K]"
    write(2,'(a)') "# Wavelength-min: N/A            # Minimum wavelength [um]"
    write(2,'(a)') "# Wavelength-max: N/A            # Maximum wavelength [um]"
    write(2,'(a)') "# Radiance-unit: W/m^2           # Radiation output unit"
    write(tmp,fmt_number1) c%T_surf
    write(2,'(a)') "# Surface-temperature: "//trim(adjustl(tmp))//"  # Surface temperature [K]"
    write(2,'(a)') "# Surface-albedo: 0.32           # Surface albedo"
    write(2,'(a)') "# Surface-emissivity: 0.9        # Surface emissivity"
    write(2,'(a)') "# Opacity: corrk + continuum     # corrk + continuum, other ? "
    write(2,'(a)') "# Line opacities = H2O, CO2"
    write(2,'(a)') "# CIA opacities = H2O-H2O, H2O-N2, CO2-CO2, N2-N2"
    write(2,'(a)') "# Rayleigh opacities = N2, CO2, H2O"
    write(tmp,'(f12.3)') OLR
    write(2,'(a)') "# Outgoing longwave radiation (W/m^2) = "//trim(adjustl(tmp))
    write(tmp,'(f12.3)') ASR
    write(2,'(a)') "# Absorbed stellar radiation (W/m^2) = "//trim(adjustl(tmp))
    write(2,'(a)') "# Atmosphere-description: Modern Earth"
    write(2,'(a)') "# Atmosphere-columns: P T MMW Alt N2 CO2 H2O OTR ASR"
    write(2,'(a)') "# Atmosphere-fit-T: N/A"
    write(2,'(a)') "# Atmosphere-fit-H2: N/A"
    write(2,'(a)') "# Atmosphere-fit-He: N/A"
    write(2,'(a)') "# Atmosphere-fit-H2O: N/A"
    write(2,'(a)') "# Atmosphere-fit-CH4: N/A"
    write(tmp,'(i12)') c%nz
    write(2,'(a)') "# Atmosphere-layers: "//trim(adjustl(tmp))

    ! labels
    labels = [ &
      'P (Pa)      ', &
      'T (K)       ', &
      'MMW (g/mol) ', &
      'Alt (km)    ', &
      'N2          ', &
      'CO2         ', &
      'H2O         ', &
      'OTR         ', &
      'ASR         ' &
    ]
    do i = 1,size(labels)
      write(2,fmt_label,advance='no') labels(i)
    enddo
    write(2,*)

    ! surface

    write(tmp,fmt_number) c%P_surf/10.0_dp
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number1) c%T_surf
    write(2,fmt_label,advance='no') adjustl(tmp)

    mmw = sum(c%f_i(1,:) * [18.0_dp, 44.0_dp, 28.0_dp])
    write(tmp,fmt_number1) mmw
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number1) 0.0_dp
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) c%f_i(1,3)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) c%f_i(1,2)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) c%f_i(1,1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) -(rad%wrk_ir%fdn_n(1) - rad%wrk_ir%fup_n(1))/1.0e3_dp
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) (rad%wrk_sol%fdn_n(1) - rad%wrk_sol%fup_n(1))/1.0e3_dp
    write(2,fmt_label,advance='no') adjustl(tmp)
    write(2,*)

    ! atmospheric layers

    do i = 1,c%nz
      write(tmp,fmt_number) c%P(i)/10.0_dp
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number1) c%T(i)
      write(2,fmt_label,advance='no') adjustl(tmp)

      mmw = sum(c%f_i(i,:) * [18.0_dp, 44.0_dp, 28.0_dp])
      write(tmp,fmt_number1) mmw
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number1) c%z(i)/1.0e5_dp
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) c%f_i(i,3)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) c%f_i(i,2)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) c%f_i(i,1)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) -(rad%wrk_ir%fdn_n(2*i+1) - rad%wrk_ir%fup_n(2*i+1))/1.0e3_dp
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) (rad%wrk_sol%fdn_n(2*i+1) - rad%wrk_sol%fup_n(2*i+1))/1.0e3_dp
      write(2,fmt_label,advance='no') adjustl(tmp)
      write(2,*)
    enddo

    close(2)

  end subroutine

end program
