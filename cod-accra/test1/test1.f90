program main
  use clima, only: Radtran, dp, s_str_len
  use clima_const, only: k_boltz
  implicit none

  type :: InputFile
    real(dp), allocatable :: alt_low(:), alt_high(:), alt(:)
    real(dp), allocatable :: P(:)
    real(dp), allocatable :: T(:)
    real(dp), allocatable :: N2(:)
    real(dp), allocatable :: CO2(:)
    real(dp), allocatable :: H2O(:)

    real(dp) :: T_surface
    real(dp), allocatable :: T_layer(:)
    real(dp), allocatable :: P_layer(:)
    real(dp), allocatable :: densities(:,:)
    real(dp), allocatable :: dz(:)
  end type

  call test()

contains

  function create_InputFile() result(input)
    type(InputFile) :: input

    integer :: io, i
    character(1000) :: tmp
    real(dp) :: alt_low, alt_high, alt, P, T, N2, CO2, H2O
    real(dp) :: density

    open(unit=2,file='../test1/input_Modern_Earth.txt',status='old',iostat=io)
    if (io /= 0) then
      print*,io
      stop 1
    endif

    do
      read(2,'(a)',iostat=io) tmp
      if (tmp(1:1) /= '#') exit
    enddo
    backspace(2)

    allocate(input%alt_low(0),input%alt_high(0))
    allocate(input%alt(0),input%P(0),input%T(0),input%N2(0),input%CO2(0),input%H2O(0))
    do
      read(2,*,iostat=io) alt_low, alt_high, alt, P, T, N2, CO2, H2O
      if (io /= 0) exit
      input%alt_low = [input%alt_low, alt_low]
      input%alt_high = [input%alt_high, alt_high]
      input%alt = [input%alt, alt]
      input%P = [input%P, P]
      input%T = [input%T, T]
      input%N2 = [input%N2, N2]
      input%CO2 = [input%CO2, CO2]
      input%H2O = [input%H2O, H2O]
    enddo

    close(2)

    input%T_surface = input%T(1)
    allocate(input%dz(size(input%alt)-1))
    allocate(input%T_layer(size(input%alt)-1))
    allocate(input%P_layer(size(input%alt)-1))
    allocate(input%densities(size(input%alt)-1,3))
    do i = 1,size(input%dz)
      input%dz(i) = (input%alt_high(i+1) - input%alt_low(i+1))*1.0e5_dp
      input%T_layer(i) = input%T(i+1)
      input%P_layer(i) = input%P(i+1)/1.0e5_dp ! bars

      density = input%P_layer(i)*1.0e6_dp/(k_boltz*input%T_layer(i))
      input%densities(i,1) = input%N2(i+1)*density
      input%densities(i,2) = input%CO2(i+1)*density
      input%densities(i,3) = input%H2O(i+1)*density
    enddo

  end function

  subroutine test()
    character(:), allocatable :: err
    type(Radtran) :: rad
    type(InputFile) :: input
    integer, parameter :: num_zenith_angles = 4
    real(dp), parameter :: surface_albedo = 0.32_dp
    real(dp), parameter :: surface_emissivity = 0.9_dp

    input = create_InputFile()

    rad = Radtran( &
      "../test1/settings_codaccra.yaml", &
      "../test1/Sun_now.txt", &
      num_zenith_angles, &
      surface_albedo, &
      size(input%dz), &
      "/Users/nicholas/Documents/Research_local/PhotochemPy/photochem_clima_data", &
      err &
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    rad%surface_emissivity = surface_emissivity

    block
      real(dp) :: stellar_radiation
      integer :: i
      stellar_radiation = 0.0_dp
      do i = 1,rad%sol%nw
        stellar_radiation = stellar_radiation + rad%photons_sol(i)*(rad%sol%freq(i) - rad%sol%freq(i+1))
      enddo
      stellar_radiation = stellar_radiation/1.0e3_dp

      ! print*,1361.0_dp/stellar_radiation

    endblock

    if (size(rad%species_names) /= 3) then
      stop 1
    endif
    if (any(rad%species_names /= ['N2 ','CO2','H2O'])) then
      stop 1
    endif

    call rad%radiate( &
      input%T_surface, &
      input%T_layer, &
      input%P_layer, &
      input%densities, &
      input%dz, &
      err=err &
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    call write_output_file(rad, input)

  end subroutine

  subroutine write_output_file(rad, input)
    type(Radtran), intent(in) :: rad
    type(InputFile), intent(in) :: input
    real(dp) :: OLR, ASR
    character(1000) :: tmp
    character(s_str_len), allocatable :: labels(:)
    character(:), allocatable :: fmt_label, fmt_number, fmt_number1
    integer :: io, i
    real(dp) :: mix(3)
    real(dp) :: mmw

    fmt_label = "(a14)"
    fmt_number = "(es14.4e3)"
    fmt_number1 = "(f14.4)"

    OLR = rad%wrk_ir%fup_n(rad%nz+1)/1.0e3_dp
    ASR = rad%wrk_sol%fdn_n(rad%nz+1)/1.0e3_dp - rad%wrk_sol%fup_n(rad%nz+1)/1.0e3_dp

    open(unit=2,file='../test1/output_ModernEarth_codaccra_photochemclima_test1.txt',status='replace',iostat=io)
    write(2,'(a)') "# Model name = photochem.clima"
    write(2,'(a)') "# Diameter: 12742                # Diameter of the planet [km]"
    write(2,'(a)') "# Gravity: 9.81                  # Gravity at the surface of the atmosphere [m/s2]"
    write(2,'(a)') "# Star-distance: 1               # Planet-star distance [ua]"
    write(2,'(a)') "# Star-type: Sun                 # Host-star type"
    write(2,'(a)') "# Star-temperature: 5772         # Host-star temperature [K]"
    write(2,'(a)') "# Wavelength-min: N/A            # Minimum wavelength [um]"
    write(2,'(a)') "# Wavelength-max: N/A            # Maximum wavelength [um]"
    write(2,'(a)') "# Radiance-unit: W/m^2           # Radiation output unit"
    write(2,'(a)') "# Surface-temperature: 288       # Surface temperature [K]"
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
    write(tmp,'(i12)') rad%nz
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

    write(tmp,fmt_number) input%P(1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number1) input%T_surface
    write(2,fmt_label,advance='no') adjustl(tmp)

    mmw = sum([input%N2(1), input%CO2(1), input%H2O(1)] * [28.0_dp, 44.0_dp, 18.0_dp])
    write(tmp,fmt_number1) mmw
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number1) input%alt(1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) input%N2(1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) input%CO2(1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) input%H2O(1)
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) -(rad%wrk_ir%fdn_n(1) - rad%wrk_ir%fup_n(1))/1.0e3_dp
    write(2,fmt_label,advance='no') adjustl(tmp)

    write(tmp,fmt_number) (rad%wrk_sol%fdn_n(1) - rad%wrk_sol%fup_n(1))/1.0e3_dp
    write(2,fmt_label,advance='no') adjustl(tmp)
    write(2,*)

    ! atmospheric layers

    do i = 1,rad%nz
      write(tmp,fmt_number) input%P_layer(i)*1.0e5_dp
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number1) input%T_layer(i)
      write(2,fmt_label,advance='no') adjustl(tmp)

      mmw = sum(input%densities(i,:)/sum(input%densities(i,:)) * [28.0_dp, 44.0_dp, 18.0_dp])
      write(tmp,fmt_number1) mmw
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number1) input%alt(i+1)
      write(2,fmt_label,advance='no') adjustl(tmp)

      mix = input%densities(i,:)/sum(input%densities(i,:))

      write(tmp,fmt_number) mix(1)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) mix(2)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) mix(3)
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) -(rad%wrk_ir%fdn_n(i+1) - rad%wrk_ir%fup_n(i+1))/1.0e3_dp
      write(2,fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt_number) (rad%wrk_sol%fdn_n(i+1) - rad%wrk_sol%fup_n(i+1))/1.0e3_dp
      write(2,fmt_label,advance='no') adjustl(tmp)
      write(2,*)
    enddo

    close(2)

  end subroutine

end program
