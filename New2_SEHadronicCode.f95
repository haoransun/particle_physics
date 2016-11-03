!save Oct/7/2014

!Removed Data from various authors. See main program subroutines for list.

!!!!!!!!!!!!!!!!!!!!!input files needed and their correct format throughout the program.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	gpep072211 is a current raw data file for experimental data points. It is preformatted from george washington
!	gpkl060311 is a current raw data file for experimental data points. It is preformatted from george washington
!
!	initial_PWAs.dat is a read file containing initial PWA values for each bin. It needs a column for each PW
!		-Format is colmn 1 is srE col 2 is siE multipole col 3 is p1rE (non physical) etc. going 
!		 up in order s p1 p3 ... and putting the real part then the imaginary part in adjacent columns
!		-It contains a column label header line that shows what each column is. 
!		-Read is using * so it ignores any amount of spaces between columns, but a space MUST exist
!
!	flag_status.dat is a read file containing flag on/off values. 
!		-This time columns are in order rE are first 10, iE are next 10, rM are next ten, iM are last 10
!		-Read is using * so it ignores any amount of spaces between columns, but a space MUST exist
!		-File contains two header lines. First as a brief synopsis of the file. second line contains the L value of ang mmtm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extact_file..._debug 9 got rid of total_array%central_energy(3:15,nE) so that it stores same bin energy for all observables.
!change nE parameter from 80 to 100 to accomodate as many bins as needed, although it looks like I only need 20

!Various things to remember:
!SAID plots variables sigma and Pol as just that. So not Sigma*Diff_Cross. Since Diff Cross function is so is sigma.
!PWA such as S11 are fitted in units mfm as of Apr 2014. This has changed in the past though. Output then changes it to mub or other units as needed including dimensionless. 



!Module holds all variables that are specific to the computer it is running on. Mostly that is pathnames for files.
!It allows you to quickly review where files are being read from and output to.
	module computer_specific
	implicit none

	Integer, PARAMETER :: NumSubDirectory=10

	!Main Pathname that most other pathnames are based off of.
	!If you change length of directory variable, other variables in program need changed as well
	character*20 ,parameter:: directory =  '/home/hsun9/Dropbox/'!Haoran1st added, the path length should match the number after character  

	!Output module path directories
	!PWA15_folderis location of energy Dependent code
	!internal folder is a subfolder inside PWA15 which is location of input files for Energy Dep Code.
	Character*13 ,parameter::PWA15_folder =         "PWA15Aug2013/"
	Character*13 ,parameter::PWA15_internal_folder= "PartialWaves/"!some folders not used in single Energy code, but used in Energy dependent code  Haoran added

	!testing subroutine
	character*4  :: pwa = "PWA/"

	end module computer_specific


!This module contains global variables that are used throughout the program. These variables are initialized in the subroutine Initialize().
!The goal is to make it easy to add new reactions. To do that, in the subroutine Initialize, assign these variables and give the reaction a number.
!That number is the number the user enters on entering the program.
	module SetupReaction
	!Folder for the given reaction to analyze. All Input/Output and needed files are found in these folders.
	character*7  :: subdirectory
	! none computer specific and filename variables for raw data and its location
	character*10 :: datafilename
	!Contains the masses of all the particles used in the photoproduction reaction. Initial nucleon plus masses of two final particles
	Double Precision :: mass_nucleon
	Double Precision :: mass_incident
	Double Precision :: mass_outgoing1
	Double Precision :: mass_outgoing2
	Double Precision :: Threshold_Energy_Wcm
	!See paper R Arndt Phys Rev C Vol 42 #5, Nov 1990 page 1854 for their definition for each reaction channel
	!H1C is constant for  Helicity 1/2, and H3C is for Helicity 3/2 coefficient. 
	!Values determined by Arndt paper. See Initialize subroutine as to what to assign the values.
	Double Precision H1C 
	Double Precision H3C

	!Isospin 0 and 1 coefficients.
	Double Precision C0
	Double Precision C1
	Double Precision C2  !Extra set of Isospin 1 amplitudes for when comparing channels that are related by 3 amplitdues such as KL-> KS vs K0 and K_bar0 reactions.

	integer :: user_reaction	!Used so program knows what reaction was used for selecting array value.

	end module SetupReaction


	module constants
	use SetupReaction
	implicit none

!PROGRAM WIDE VARIABLES WARNING DO NOT USE THESE VARIABLES IN AN INDIVIDUAL MODULE OR SUB PROGRAM, ESPECIALLY X
!WARNING USED IN ALL MODULES
	Double Precision :: X !Cos(Theta), This is for Helicity amplitudes as well.
	integer, parameter :: num_observables = 4 !16 photo amplitudes plus 3 DSG helicity 1/2 and 3/2 values plus extra observables multiplied by DSG eg. E*DSG or G*DSG.
	integer :: bin_num !This is the bin_num you are currently on- This is used in all Helicity Functions.
	integer :: num_bins !This is a counter than holds final count of Single Energy Bins.
	logical :: TESTING_VARIABLE = .false.
	!value to keep track of todays date for plots/fileInfo and other needed items
	integer :: datetoday(3)
	
	Double Precision, parameter :: PI=3.1415926535897932
	Double Precision, parameter :: h_barc = 197.33
	Double Precision, parameter :: h_barc2 = h_barc**2
	integer, 	  PARAMETER :: nE=100!!! maximum number of bins-same variable as num_bins but for variable declaration below, 
!        integer, 	  PARAMETER :: NumberOfData_Maximum=50

	!bin variables for intial final bin number and bin width used as a single energy.
	integer           :: end_bin = 9  !this variable sets max bins processed. used for forward/backward direction
	double precision  :: bin_width(nE) 	!this gives the range of energy values treated as 1 energy

	!Most unit numbers used throughout the program. 
	!some unit output variables in 3000 range are found in output module
	!table i is used with tablei+number
	!It is my intention to combine some of these so it has fewer numbers and they are reused for the same purpose.
	integer, parameter :: PWA_open = 25

	!UNIT 1000 is the main data input file for reactions.
	integer :: SGT1i=1001, SGT3i=1002 !SGT 1/2 and 3/2 Helicity contributions 
	integer	:: SGTi=1003,DSGi=1004,POLi=1005 
	integer :: PDSGi=1006,Fi=1007 

	integer :: tmp=998
	integer :: tablei=1030 !allows generation of total array for quality assurance
	INTEGER :: Unit_OUTPUT=3000, UNIT_num_bins = 3020 !found in output module
	Integer :: ChiPerPoint_file= 7000
	integer :: output_amplitude_modulus = 7001
	integer :: unit_pwa_input_initial = 5000
	integer :: unit_phase = 8000
	integer :: unit_flag_status_file = 5002
	integer :: unit_debug = 50
        integer :: unit_s001r=8000,unit_p001r=8001,unit_p003r=8002,unit_p103r=8003,&
		   unit_d003r=8004,unit_d103r=8005,unit_d005r=8006,unit_d105r=8007,&
		   unit_f005r=8008,unit_f105r=8009,unit_f007r=8010,unit_f107r=8011,&
		   unit_g007r=8012,unit_g107r=8013,unit_g009r=8014,unit_g109r=8015,&
		   unit_h009r=8016,unit_h109r=8017,unit_h011r=8018,unit_h111r=8019

        integer :: unit_s001i=8020,unit_p001i=8021,unit_p003i=8022,unit_p103i=8023,&
		   unit_d003i=8024,unit_d103i=8025,unit_d005i=8026,unit_d105i=8027,&
		   unit_f005i=8028,unit_f105i=8029,unit_f007i=8030,unit_f107i=8031,&
		   unit_g007i=8032,unit_g107i=8033,unit_g009i=8034,unit_g109i=8035,&
		   unit_h009i=8036,unit_h109i=8037,unit_h011i=8038,unit_h111i=8039

	integer :: unit_s101r=8041,unit_s101i=8042,unit_p101r=8043,unit_p101i=8044

        integer :: unit_s101pe=8070,unit_p101pm=8071,unit_p103pe=8072,unit_p103pm=8073,&
		   unit_d103pe=8074,unit_d103pm=8075,unit_d105pe=8076,unit_d105pm=8077,&
		   unit_f105pe=8078,unit_f105pm=8079,unit_f107pe=8080,unit_f107pm=8081,&
		   unit_g107pe=8082,unit_g107pm=8083,unit_g109pe=8084,unit_g109pm=8085,&
		   unit_h109pe=8086,unit_h109pm=8087,unit_h111pe=8088,unit_h111pm=8089

	integer :: unit_pathnamefileObservable=8040
	integer :: unit_all_amplitudes=8041
	!testing unit number
	integer :: gw_data_file=9504

	character*7, dimension(1:3+num_observables) :: observable_files= &
				(/"SGT1   ","SGT3   ","SGT    ","DSG    ","P      ","PDSG   ",&
				"F      "/)

	!headerlines and hence graphs contain a label unit which may change in the future. 
	!Observables will be in units microbarnes or mub^2
	!PWA's such as S11 will be in units mfm.
	character*6 :: units_mub  = "(#mub)"
	character*7 :: units_mb  = "(mb/sr)"
	character*5 :: units_mfm = "(mfm)"
	character*9 :: units_mfm2 = "(mfm)^{2}"

	!197.327 is h_bar*c constant. 10^6 changes mfm^2 to fm^2. 100 changes micro barnes to mfm^2-reminder only for mod squared. coeff does not need a 10 because it is used after s11 is converted to mfm
	double precision, parameter :: coeff  = 1./(197.327*10**3)	
	double precision :: coeff2 =coeff*coeff

	!Const that will =1 when PW is read as Mod Phase
	!Const that will= 1/cos(45 deg) when read as Rl Im. This means the phase variable is not applicable.
	double precision :: VCC(40)
	double precision :: SGTC(40,50)!40 pw's and 50 bins
	end module constants


!this module holds the main database variable for data storage. It holds all raw data, energy bin info, and other information needed
!to output information to files, analyze the data, and bin it.
	module types
	use constants	
	
!holds all angles, diff cross section and errors for data(1 of each per row), and a number related to author
		type datapoints  
		double precision :: energies_datapoints(1500,4)=0! 
		double precision :: exact_energy_val(1500)=0!exact energy
		integer :: author_number(1500) !holds a number related to each author for plotting author data with different symbols.
		character*15 :: data_points_author_name(1500)
		character*60 :: data_points_paper_date(1500)
		character*4 :: reactionID(1500)
		end type datapoints

		!row gives the different observables, column is the different energy bins.
		type energies
		double precision :: central_energy(nE)=0   !this holds the central value of the energy bin
		integer :: energy_num_points(1:3+num_observables, nE)=0        !holds the # of points in a bin for a given observable
		type (datapoints) :: energy_values(1:3+num_observables, nE)   !holds data (angles, OBS, OBS error, and actual energy)
		end type energies
		
		!variable of type energies declaration
		type(energies) :: total_array
!toatal_array%energy_values(4,1)%energy_datapoints(21,2),this is an example about access the data.
	end module types

!module of functions for various Legendre Polynomials and their first and second derivative.
!returns for a given order and X value the value of P, P' or P" at X
	module Legendre_Poly
	use constants, only:X !only statement means all other system wide variables from constants are UNKNOWN to this module
	implicit none

	Contains
		
		!Definition of Legendre Polynomials up to order 6. N is the desired order that is passed in.
		!X is a system wide variable that needs to be SET FIRST. 
		!It is defined in usual way X=cos(A) where A is angle in radians
		!function Legendre(n)!returns value of legendre poly for a given angle and angular momentum expansion l		
		!	integer :: n
		!	Double Precision, Dimension(0:6) :: P !definition variables for first 6 legendre polynomials
		!	Double Precision :: Legendre
			
		!	P(0)=1.
        	!	P(1)=X
        	!	P(2)=0.5*(3.*X*X - 1.)
        	!	P(3)=0.5*(5.*X*X*X - 3.*X)
        	!	P(4)=0.125*(35.*X*X*X*X - 30.*X*X + 3.)
        	!	P(5)=0.125*(63.*X*X*X*X*X - 70.*X*X*X + 15.*X)
        	!	P(6)=0.0625*(231.*x*x*x*x*x*x - 315.*x*x*x*x + 105.*x*x - 5.)
			
		!	Legendre = P(n)

		!end function Legendre

		!Definition of Legendre Polynomials up to order 6. N is the desired order that is passed in.
		!X is a system wide variable that needs to be SET FIRST. 
		!It is defined in usual way X=cos(A) where A is angle in radians
		function d_L(n)  !returns value of first derivative of legendre polys
			integer :: n
			Double Precision, Dimension(0:6) :: dP !definition variables for the first 6 derivatives of those polynomials
			Double Precision d_L
			
			dP(0)=0.
			dP(1)=1.
			dP(2)=0.5*(6.*X)
			dP(3)=0.5*(15.*X*X - 3.)
			dP(4)=0.125*(35.*4.*X*X*X - 30.*2.*X)
			dP(5)=0.125*(63.*5.*X*X*X*X - 70.*3.*X*X + 15.)
			dP(6)=.0625*(231.*6.*X*X*X*X*X - 315.*4.*X*X*X + 105.*2.*X)

			d_L=dP(n)

		end function d_L

		!Definition of 2nd Derivative of Legendre Polynomials up to order 6. N is the desired order that is passed in.
		!X is a system wide variable that needs to be SET FIRST. 
		!It is defined in usual way X=cos(A) where A is angle in radians
		function d2_L(n) !returns value of second derivative of legendre poly.
			integer :: n
			Double Precision, Dimension(0:6) :: d2P	!second derivative
			Double Precision d2_L
			
			d2P(0)=0.
			d2P(1)=0.
			d2P(2)=0.5*6.
			d2P(3)=0.5*(15.*2.*X)
			d2P(4)=0.125*(35.*4.*3.*X*X - 30.*2.)
			d2P(5)=0.125*(63.*5.*4.*X*X*X - 70.*3.*2.*X)
			d2P(6)=	.0625*(231.*6.*5.*X*X*X*X - 315.*4.*3.*X*X + 105.*2.)	

			d2_L=d2P(n)

		end function d2_L

	end module Legendre_Poly


!Module Devoted to the Helicity Amplitudes---------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------
	module Helicity_Amplitudes
	use constants

        	!Iso 0 Expansion Coefficients   
        	Double Precision :: W(ne),&
          	s01a(ne), s01b(ne), &
          	p01a(ne), p01b(ne), p03a(ne), p03b(ne),&
          	d03a(ne), d03b(ne), d05a(ne), d05b(ne),&
          	f05a(ne), f05b(ne), f07a(ne), f07b(ne),&
          	g07a(ne), g07b(ne), g09a(ne), g09b(ne),& 
          	h09a(ne), h09b(ne), h011a(ne),h011b(ne)   

		Data W(ne), &
          	s01a(ne), s01b(ne), &
          	p01a(ne), p01b(ne), p03a(ne), p03b(ne),&
          	d03a(ne), d03b(ne), d05a(ne), d05b(ne),&
          	f05a(ne), f05b(ne), f07a(ne), f07b(ne),&
          	g07a(ne), g07b(ne), g09a(ne), g09b(ne),& 
          	h09a(ne), h09b(ne), h011a(ne),h011b(ne) /23*0.00/

		!Iso 1 Expansion Coefficients
          	Double Precision :: &
          	s11a(ne), s11b(ne), &
          	p11a(ne), p11b(ne), p13a(ne), p13b(ne),&
          	d13a(ne), d13b(ne), d15a(ne), d15b(ne),&
          	f15a(ne), f15b(ne), f17a(ne), f17b(ne),&
          	g17a(ne), g17b(ne), g19a(ne), g19b(ne),& 
          	h19a(ne), h19b(ne), h111a(ne),h111b(ne)
	
		Data        &    	
                s11a(ne), s11b(ne), &
          	p11a(ne), p11b(ne), p13a(ne), p13b(ne),&
          	d13a(ne), d13b(ne), d15a(ne), d15b(ne),&
          	f15a(ne), f15b(ne), f17a(ne), f17b(ne),&
          	g17a(ne), g17b(ne), g19a(ne), g19b(ne),& 
          	h19a(ne), h19b(ne), h111a(ne),h111b(ne)/22*0.00/

		!Extra set of Iso 1 for KL->KS reaction which is a fuction of 3 amplitudes when comparing to other channels
          	Double Precision :: &
          	s11c(ne), s11d(ne), &
          	p11c(ne), p11d(ne), p13c(ne), p13d(ne),&
          	d13c(ne), d13d(ne), d15c(ne), d15d(ne),&
          	f15c(ne), f15d(ne), f17c(ne), f17d(ne),&
          	g17c(ne), g17d(ne), g19c(ne), g19d(ne),& 
          	h19c(ne), h19d(ne), h111c(ne),h111d(ne)
	
		Data        &    	
                s11c(ne), s11d(ne), &
          	p11c(ne), p11d(ne), p13c(ne), p13d(ne),&
          	d13c(ne), d13d(ne), d15c(ne), d15d(ne),&
          	f15c(ne), f15d(ne), f17c(ne), f17d(ne),&
          	g17c(ne), g17d(ne), g19c(ne), g19d(ne),& 
          	h19c(ne), h19d(ne), h111c(ne),h111d(ne)/22*0.00/

		!Iso 0 Mod/phase values
        	Double Precision :: &  
          	s01m(ne), s01p(ne), &
          	p01m(ne), p01p(ne), p03m(ne), p03p(ne),&
          	d03m(ne), d03p(ne), d05m(ne), d05p(ne),&
          	f05m(ne), f05p(ne), f07m(ne), f07p(ne),&
          	g07m(ne), g07p(ne), g09m(ne), g09p(ne),& 
          	h09m(ne), h09p(ne), h011m(ne),h011p(ne)   

		Data &
          	s01m(ne), s01p(ne), &
          	p01m(ne), p01p(ne), p03m(ne), p03p(ne),&
          	d03m(ne), d03p(ne), d05m(ne), d05p(ne),&
          	f05m(ne), f05p(ne), f07m(ne), f07p(ne),&
          	g07m(ne), g07p(ne), g09m(ne), g09p(ne),& 
          	h09m(ne), h09p(ne), h011m(ne),h011p(ne)   /22*0.00/

		!Iso 1 Mod/phase values
          	Double Precision :: &
          	s11m(ne), s11p(ne), &
          	p11m(ne), p11p(ne), p13m(ne), p13p(ne),&
          	d13m(ne), d13p(ne), d15m(ne), d15p(ne),&
          	f15m(ne), f15p(ne), f17m(ne), f17p(ne),&
          	g17m(ne), g17p(ne), g19m(ne), g19p(ne),& 
          	h19m(ne), h19p(ne), h111m(ne),h111p(ne)
	
		Data &
          	s11m(ne), s11p(ne), &
          	p11m(ne), p11p(ne), p13m(ne), p13p(ne),&
          	d13m(ne), d13p(ne), d15m(ne), d15p(ne),&
          	f15m(ne), f15p(ne), f17m(ne), f17p(ne),&
          	g17m(ne), g17p(ne), g19m(ne), g19p(ne),& 
          	h19m(ne), h19p(ne), h111m(ne),h111p(ne)/22*0.00/


	contains

	Function LambdaBar2()
	use SetupReaction
	use constants
	double precision :: lambdabar2
	double precision :: W2
	double precision :: mass_plus2, mass_minus2

	W2 = W(bin_num)**2
        mass_plus2 =(Mass_incident+Mass_nucleon)**2
        mass_minus2=(Mass_incident-Mass_nucleon)**2

        Q2=((W2-Mass_plus2)*(W2-Mass_minus2))/(4.*W2)
        lambdaBar2=(10.*h_barc2)/Q2

	return
	end

        FUNCTION f_Real(X)
!*-----------------------------------------------------------------------
!C       REAL part of Spin-Non-Flip amplitude f at C.M. energy E
!*-----------------------------------------------------------------------
        double precision :: f_Real

        double precision :: X 
        double precision :: P0, P1, P2, P3, P4, P5
        double precision :: f_Real_0, f_Real_1, f2_Real_1
!*-----------------------------------------------------------------------
        f_Real=0

        P0=1.
        P1=X
        P2=0.5*(3.*X*X-1.)
        P3=0.5*(5.*X*X*X-3.*X)
        P4=0.125*(35.*X*X*X*X-30.*X*X+3.)
        P5=0.125*(63.*X*X*X*X*X-70.*X*X*X+15.*X)
!*-----------------------------------
        f_Real_0=&
               (1.*s01a(bin_num)+0.*0.000)*P0+&
               (2.*p03a(bin_num)+1.*p01a(bin_num))*P1+&
               (3.*d05a(bin_num)+2.*d03a(bin_num))*P2+&
               (4.*f07a(bin_num)+3.*f05a(bin_num))*P3+&
               (5.*g09a(bin_num)+4.*g07a(bin_num))*P4+&
               (6.*0.000        +5.*h09a(bin_num))*P5
    
        f_Real_1=&
               (1.*s11a(bin_num)+0.*0.000)*P0+&
               (2.*p13a(bin_num)+1.*p11a(bin_num))*P1+&
               (3.*d15a(bin_num)+2.*d13a(bin_num))*P2+&
               (4.*f17a(bin_num)+3.*f15a(bin_num))*P3+&
               (5.*g19a(bin_num)+4.*g17a(bin_num))*P4+&
               (6.*0.000        +5.*h19a(bin_num))*P5
    
        f2_Real_1=&
               (1.*s11c(bin_num)+0.*0.000)*P0+&
               (2.*p13c(bin_num)+1.*p11c(bin_num))*P1+&
               (3.*d15c(bin_num)+2.*d13c(bin_num))*P2+&
               (4.*f17c(bin_num)+3.*f15c(bin_num))*P3+&
               (5.*g19c(bin_num)+4.*g17c(bin_num))*P4+&
               (6.*0.000        +5.*h19c(bin_num))*P5
   
        f_Real=&
               C0*f_Real_0+&
               C1*f_Real_1+&
               C2*f2_Real_1
           !print *,"freal:",f_real,c1
        RETURN
        END

        FUNCTION f_Imaginary(X)
!*-----------------------------------------------------------------------
!C       IMAGINARY part of Spin-Non-Flip amplitude f at C.M. energy E
!*-----------------------------------------------------------------------
        double precision :: f_Imaginary

        double precision :: X

        double precision :: P0, P1, P2, P3, P4, P5
        double precision :: f_Imaginary_0, f_Imaginary_1, f2_Imaginary_1
!*-----------------------------------------------------------------------
        f_Imaginary=0

        P0=1.
        P1=X
        P2=0.5*(3.*X*X-1.)
        P3=0.5*(5.*X*X*X-3.*X)
        P4=0.125*(35.*X*X*X*X-30.*X*X+3.)
        P5=0.125*(63.*X*X*X*X*X-70.*X*X*X+15.*X)
!*-----------------------------------
        f_Imaginary_0=&
               (1.*s01b(bin_num)+0.*0.000)*P0+&
               (2.*p03b(bin_num)+1.*p01b(bin_num))*P1+&
               (3.*d05b(bin_num)+2.*d03b(bin_num))*P2+&
               (4.*f07b(bin_num)+3.*f05b(bin_num))*P3+&
               (5.*g09b(bin_num)+4.*g07b(bin_num))*P4+&
               (6.*0.000        +5.*h09b(bin_num))*P5
     
        f_Imaginary_1=&
               (1.*s11b(bin_num)+0.*0.000)*P0+&
               (2.*p13b(bin_num)+1.*p11b(bin_num))*P1+&
               (3.*d15b(bin_num)+2.*d13b(bin_num))*P2+&
               (4.*f17b(bin_num)+3.*f15b(bin_num))*P3+&
               (5.*g19b(bin_num)+4.*g17b(bin_num))*P4+&
               (6.*0.000        +5.*h19b(bin_num))*P5
      
        f2_Imaginary_1=&
               (1.*s11d(bin_num)+0.*0.000)*P0+&
               (2.*p13d(bin_num)+1.*p11d(bin_num))*P1+&
               (3.*d15d(bin_num)+2.*d13d(bin_num))*P2+&
               (4.*f17d(bin_num)+3.*f15d(bin_num))*P3+&
               (5.*g19d(bin_num)+4.*g17d(bin_num))*P4+&
               (6.*0.000        +5.*h19d(bin_num))*P5
       
        f_Imaginary=&
               C0*f_Imaginary_0+&
               C1*f_Imaginary_1+&
               C2*f2_Imaginary_1
         !print *,"fimaginary:",f_Imaginary
        RETURN
        END


        FUNCTION g_Real(X)
!*-----------------------------------------------------------------------
!C       REAL part of Spin-Flip amplitude g at C.M. energy E
!*-----------------------------------------------------------------------
        double precision :: g_Real

        double precision :: X
        double precision :: B0, B1, B2, B3, B4, B5
        double precision :: g_Real_0, g_Real_1, g2_Real_1
!*-----------------------------------------------------------------------
        g_Real=0

        B0=0.
        B1=sqrt(1.-X*X)
        B2=3.*X*B1
        B3=0.5*(15.*X*X-3.)*B1
        B4=0.5*(35.*X*X*X-15.*X)*B1
        B5=0.125*(315.*X*X*X*X-210.*X*X+15.)*B1
!*-----------------------------------
        g_Real_0=&
               (p03a(bin_num)-p01a(bin_num))*B1+&
               (d05a(bin_num)-d03a(bin_num))*B2+&
               (f07a(bin_num)-f05a(bin_num))*B3+&
               (g09a(bin_num)-g07a(bin_num))*B4+&
               (0.000        -h09a(bin_num))*B5
    
        g_Real_1=&
               (p13a(bin_num)-p11a(bin_num))*B1+&
               (d15a(bin_num)-d13a(bin_num))*B2+&
               (f17a(bin_num)-f15a(bin_num))*B3+&
               (g19a(bin_num)-g17a(bin_num))*B4+&
               (0.000        -h19a(bin_num))*B5
    
        g2_Real_1=&
               (p13c(bin_num)-p11c(bin_num))*B1+&
               (d15c(bin_num)-d13c(bin_num))*B2+&
               (f17c(bin_num)-f15c(bin_num))*B3+&
               (g19c(bin_num)-g17c(bin_num))*B4+&
               (0.000        -h19c(bin_num))*B5
          
        g_Real=&
               C0*g_Real_0+&
               C1*g_Real_1+&
               C2*g2_Real_1
         ! print *,"greal:",g_real
        RETURN
        END



        FUNCTION g_Imaginary(X)
!*-----------------------------------------------------------------------
!C       IMAGINARY part of Spin-Flip amplitude g at C.M. energy E
!*-----------------------------------------------------------------------
        REAL g_Imaginary

        double precision :: X
        double precision :: B0, B1, B2, B3, B4, B5
        double precision :: g_Imaginary_0, g_Imaginary_1, g2_Imaginary_1
!*-----------------------------------------------------------------------
        g_Imaginary=0

        B0=0.
        B1=sqrt(1.-X*X)
        B2=3.*X*B1
        B3=0.5*(15.*X*X-3.)*B1
        B4=0.5*(35.*X*X*X-15.*X)*B1
        B5=0.125*(315.*X*X*X*X-210.*X*X+15.)*B1
!*-----------------------------------
        g_Imaginary_0=&
               (p03b(bin_num)-p01b(bin_num))*B1+&
               (d05b(bin_num)-d03b(bin_num))*B2+&
               (f07b(bin_num)-f05b(bin_num))*B3+&
               (g09b(bin_num)-g07b(bin_num))*B4+&
               (0.000        -h09b(bin_num))*B5
        
        g_Imaginary_1=&
               (p13b(bin_num)-p11b(bin_num))*B1+&
               (d15b(bin_num)-d13b(bin_num))*B2+&
               (f17b(bin_num)-f15b(bin_num))*B3+&
               (g19b(bin_num)-g17b(bin_num))*B4+&
               (0.000        -h19b(bin_num))*B5
         
        g2_Imaginary_1=&
               (p13d(bin_num)-p11d(bin_num))*B1+&
               (d15d(bin_num)-d13d(bin_num))*B2+&
               (f17d(bin_num)-f15d(bin_num))*B3+&
               (g19d(bin_num)-g17d(bin_num))*B4+&
               (0.000        -h19d(bin_num))*B5
      
        g_Imaginary=&
               C0*g_Imaginary_0+&
               C1*g_Imaginary_1+&
               C2*g2_Imaginary_1
     ! print *,"gimaginary:",g_Imaginary
        RETURN
        END


!
!
!Definition of the various observable functions in terms of the Helicity Amplitudes. These are the actual functions that are called.
!The current units of the functions, and hence the PWA's is micro barnes to match the units of the experimental data.
!this computes diff_cross(cos(Theta))
        FUNCTION FUNC_Diff_Cross_Section()
	use constants, only:X
!*-----------------------------------------------------------------------
!C       differential cross section at X=Cos(Theta)
!*-----------------------------------------------------------------------
        double precision :: FUNC_Diff_Cross_Section
!*-----------------------------------------------------------------------
        FUNC_Diff_Cross_Section=0

        FUNC_Diff_Cross_Section=lambdaBar2()*( &
     		f_Real(X)**2+f_Imaginary(X)**2+ &
     		g_Real(X)**2+g_Imaginary(X)**2)

        RETURN
        END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!this computes Polarization observable
	function FUNC_Polarization()
		Double Precision :: FUNC_Polarization
		Complex*16 :: temp_Helicity
			!CALL Interpolate(Ecm)

		FUNC_Polarization = 0

		FUNC_Polarization = FUNC_PDSG()/FUNC_Diff_Cross_Section()
	return
	end function FUNC_Polarization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!this computes Sigma observable
	function FUNC_PDSG()
	use constants, only:X
		Double Precision :: FUNC_PDSG 
		Complex*16 :: temp_PDSG
			
		FUNC_PDSG = 0
        	FUNC_PDSG=2.*lambdaBar2()*(&
     		          f_Imaginary(X)*g_Real(X)-&
     		          g_Imaginary(X)*f_Real(X))
	return
	end function FUNC_PDSG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!this computes T observable
	function FUNC_T()
		Double Precision :: FUNC_T
		Complex*16 :: temp_Helicity

		FUNC_T = 0
                       
	return
	end function FUNC_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!this computes G observable
	function FUNC_G()
		Double Precision :: FUNC_G 
		Complex*16 :: temp_Helicity

		FUNC_G = 0

		!FUNC_G = (-1.)*dimag(temp_Helicity)*outgoing_coeff_q()/(incident_coeff_k()*FUNC_Diff_Cross_Section())

	return
	end function FUNC_G
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!this computes H observable
		function FUNC_H()
			Double Precision :: FUNC_H
			Complex*16 :: temp_Helicity

			FUNC_H = 0.0

			!FUNC_H = (-1.)*dimag(temp_Helicity)*outgoing_coeff_q()/(incident_coeff_k()*FUNC_Diff_Cross_Section())
		return
		end function FUNC_H
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		function FUNC_F()
			Double Precision :: FUNC_F
			Complex*16 :: temp_Helicity

			FUNC_F = 0.0

			!FUNC_F = real(temp_Helicity)*outgoing_coeff_q()/(incident_coeff_k()*FUNC_Diff_Cross_Section())
		return
		end function FUNC_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		function FUNC_E()
			Double Precision :: FUNC_E
			Double Precision :: temp_Helicity

			FUNC_E = 0.0


			!FUNC_E = temp_Helicity*outgoing_coeff_q()/(2*incident_coeff_k()*FUNC_Diff_Cross_Section())
		return
		end function FUNC_E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end module Helicity_Amplitudes



	module Chi_Square_mod
		use computer_specific
		use constants
 		use Helicity_Amplitudes
		use types

        	integer, PARAMETER :: NumberOfParameter_Maximum=40	!when {s01a,...,h09b} are all free parameters as are {s11a,...,h19b}
		integer, PARAMETER :: NumberOfFlag=NumberOfParameter_Maximum
		integer, PARAMETER :: NumberOfStack_Maximum=700 !number of one observable in one bin max value.


     		DOUBLE PRECISION :: FreeParameter(NumberOfParameter_Maximum)
     		DOUBLE PRECISION :: Error_FreeParameter(NumberOfParameter_Maximum)
     		DOUBLE PRECISION :: Pool_1(NumberOfParameter_Maximum)
        	DOUBLE PRECISION :: Pool_1_Error(NumberOfParameter_Maximum)
        	DOUBLE PRECISION :: Pool_1_Original(NumberOfParameter_Maximum)
		DOUBLE PRECISION :: Pool_1_modphase(NumberOfParameter_Maximum)
		Double Precision :: Pool_1_MP_Error(NumberOfParameter_Maximum)
		Double Precision :: Pool_1_Mod_Original(NumberOfParameter_Maximum) !holds starting mod value from initial parameters
		Double Precision :: ED_Mod(ne,NumberOfParameter_Maximum)     !Holds the ED fit mod value in cases where I iterate multiple times. Only has mod values
        	DOUBLE PRECISION :: Error_Default
           
!Variables for ChiSquareFit and Output CHiSquare
	        CHARACTER*4 :: PartialWave_1(NumberOfFlag/2)	
	        DATA PartialWave_1 /'S01','P01','P03','D03','D05','F05','F07','G07','G09','H09',&
				    'S11','P11','P13','D13','D15','F15','F17','G17','G19','H19'/!,&
				    !'S31E','P31E','P33E','D33E','D35E','F35E','F37E','G37E','G39E','H39E',&
				    !'S31M','P31M','P33M','D33M','D35M','F35M','F37M','G37M','G39M','H39M'/
                CHARACTER*2 :: VaryCode(0:2)
                DATA VaryCode /' 0','90', '45'/	!0: fixed as Rl/Im  90: varied as Rl/Im by D. Mark Manley's convention 45:PhaseFixed Mod varied BCH convention
!*           DATA VaryCode /'0' , '1'/	!0: fixed   1: varied
!*           DATA VaryCode /'F' , 'V'/	!F: fixed   V: varied 
   
       		INTEGER :: Flag_1(NumberOfFlag), Flag_1_Original(nE,NumberOfFlag) !original handles turning on/off flags for printing
        	INTEGER, parameter :: Flag_Off=0, Flag_On=1, Flag_MP_On=2
		INTEGER            :: FLAG_PNLTY
		double precision :: Pnlty_Strength(NumberOfFlag) !Holds strength of penalty contribution for PW. Rl and Im Parts are read in as same

                
		!Chi Square Summation Variables
		INTEGER :: Flag_last_ChiSummation = 0   !when 1, chi square writes to file the chi per point values 
		Double Precision :: ChiSummation(3:3+num_observables) = 0  !array holds observables total chi value
		Double Precision :: GLOBAL_ChiSum(3:3+num_observables) = 0 !array holds the total chi square summed over all bins (each observable get its own sum)
		Integer :: IterNum = 0!bug testing variable
		!Fit values and Data Points.
     		DOUBLE PRECISION :: Y(NumberOfStack_Maximum)                !holds datapoints of diff cross section
     		DOUBLE PRECISION :: ERROR_Y(NumberOfStack_Maximum)	!hold datapoints of error of diff cross
     		DOUBLE PRECISION :: Y_FIT(NumberOfStack_Maximum)	!holds best fit data
		

     		INTEGER :: NumberOfData = 0
		INTEGER :: NumberOfParameter
                integer :: numberofparameter_original

		Integer :: dataForObs(3:3+num_observables)

     		DOUBLE PRECISION :: ARRAY(NumberOfParameter_Maximum, NumberOfParameter_Maximum)
     		DOUBLE PRECISION :: dF(NumberofStack_Maximum, NumberOfParameter_Maximum) = 0


	contains

	FUNCTION FUNCTION_FittingFunction(fit_bin_number,fit_ObsNum, COSTHETA,data_number)
		!*-----------------------------------------------------------------------
		!*       fitting function for observables
		!*-----------------------------------------------------------------------

		!function variable
        	DOUBLE PRECISION FUNCTION_FittingFunction
	
		!passed variables
        	INTEGER fit_ObsNum, fit_bin_number
        	DOUBLE PRECISION  COSTHETA
		Integer data_number !number of a data in the bin
		!local variables
        	DOUBLE PRECISION fit_obs
		!*-----------------------------------------------------------------------
        	FUNCTION_FittingFunction=0.

		X=COSTHETA
		bin_num=fit_bin_number 

		call setisospincoeff(total_Array%energy_values(fit_ObsNum,bin_num)%reactionID(data_number))

        	SELECT CASE(fit_ObsNum)
		CASE (3)
			!fit_obs = func_Total_Cross_Section()
			!FUNCTION_FittingFunction=fit_obs
		
        	CASE (4) 
            		fit_obs  = FUNC_Diff_Cross_Section()
            		FUNCTION_FittingFunction=fit_obs

		CASE (5)
			fit_obs = FUNC_Polarization()
			FUNCTION_FittingFunction=fit_obs

		CASE (6)
			fit_obs = FUNC_PDSG()
			FUNCTION_FittingFunction=fit_obs

    		CASE (9)
			fit_obs = FUNC_F()
			FUNCTION_FittingFunction=fit_obs
        	CASE DEFAULT
            		PRINT *, 'Observable type= ', fit_Obsnum
            		PRINT *, 'New observable? see function_FittingFunction'
            		RETURN
        	END SELECT

        RETURN
        END FUNCTION FUNCTION_FittingFunction

!assigns the Isospin 1/2 and 3/2 coefficients. This is a combination of the isospin coefficient needed to go from Energy Dependent to Single Energy, as well
!as the Isospin coefficient needed for the charge channel in the case of mixed isospin reactions.
!See Niboh dissertation pg 23 for The Energy Dep -> Single Energy conversion (dissertation gives SE->ED conversion).
!See Arndt et.al. Phys Rev C Vol:42 #5(Nov) 42(1990) for the decomposition into Isospin 1/2 and 3/2 
!So Iso1 coeff from Niboh (ED->SE) is -1/Sqrt(3) and GP-> Pi0 p pArndt paper says H1/2p is 1 so total H1C is -1/Sqrt(3) for that reaction
!exmaple2: Iso3/2: Niboh (ED->SE)	  is Sqrt(3/2) and GP->Pi+ n arndt paper says H3/2 is -Sqrt(2)/3 so together gives -1/Sqrt(3) for H3C on GP->pi+ n
	SUBROUTINE setIsoSpinCoeff(reaction_ID)
	character *4 :: reaction_ID
	!print *,"reaction id:",reaction_id
		if(reaction_ID .eq. "GPEP") then
		  H1C = 1
		  H3C = 0
		else if(reaction_ID .eq. "K+L0") then
		  H1C = 1
		  H3C = 0
		else if(reaction_ID .eq. "PI+") then
		  H1C = -DSQRT(DBLE(2.0/3.0))
		  H3C = -DSQRT(DBLE(1./3.0))
		else if(reaction_ID .eq. "PI0") then
		  H1C = -DSQRT(DBLE(1.0/3.0))
		  H3C =  DSQRT(DBLE(2.0/3.0))
		else if(reaction_ID .eq. "GNEN") then
		  H1C = 1
		  H3C = 0
		else if(reaction_ID .eq. "PI-") then
		  H1C = -DSQRT(DBLE(2.0/3.0))
		  H3C =  DSQRT(DBLE(1./3.0))
		else if(reaction_ID .eq. "PI0N") then
		  H1C =  DSQRT(DBLE(1.0/3.0))
		  H3C =  DSQRT(DBLE(2.0/3.0))
		else if(reaction_ID .eq. "K-pE") then
		  C0 =  0.5
		  C1 =  0.5
		  C2 =  0.0
		else if(reaction_ID .eq. "KLKS") then
		  C0 =  0.25
		  C1 =  0.25
		  C2 =  -0.5!HaoranMultiple
		else if(reaction_ID .eq. "K-Kb") then
		  C0 =  -0.5
		  C1 =  0.5
		  C2 =  0!HaoranMultiple
		else
		  print *,"Reaction Helicity Coefficients not setup in subroutine: setIsoSpinCoeff"
		  stop "Program Terminates to allow reaction to be added"
		end if
	RETURN
	END

	FUNCTION ChiSquare(local_Y, local_Error_Y, local_Y_Fit)!-J was passed in to give datapoint number to look at.
		!*-----------------------------------------------------------------------
		!C       calculate Chi-Square of each fitting data
		!*-----------------------------------------------------------------------
        	DOUBLE PRECISION ChiSquare

		!passed variables
        	DOUBLE PRECISION local_Y, local_Error_Y, local_Y_Fit
        	DOUBLE PRECISION Chi
		!*-----------------------------------------------------------------------
        	ChiSquare=0.
        	IF (local_Error_Y.NE.0.) THEN
                Chi=(local_Y-local_Y_Fit)/local_Error_Y
                ChiSquare=Chi*Chi
        	END IF

		if(Flag_last_ChiSummation .eq. 1) Then
		write(ChiperPoint_file,'(f8.3,2x,f8.3,2x,f8.3,2x, f12.2)')&
			, local_Y_Fit, local_Y, local_Error_Y, ChiSquare
		End If
        RETURN
        END FUNCTION ChiSquare

!Definition of Penalty term is simply a Chi Square term [(Starting_PW - Fit_PW)/Error]**2 reinterpreted
!A strong Error value was found to be around .01 or 1/Error**2 being 10,000
!A weak Error value was found to be .1 or 1/Error**2 around 100.
!By defining 1/Error**2 = Pnlty_Strength*100, we find A weak strength ~= Pnlty_strength=1
!A strong Pnlty_Strength is then 100.
!While a strength of 0 implies this PW has no penalty term which is intuitive.	
	FUNCTION ChiSquare_Penalty()
		!*-----------------------------------------------------------------------
		!C       compute Chi-Square penalty for amplitude set
		!*-----------------------------------------------------------------------
       		DOUBLE PRECISION ChiSquare_Penalty
								
		double precision Pnlty_dX

        	ChiSquare_Penalty= 0.0

		do I=1,40
		  ChiSquare_Penalty = ChiSquare_Penalty + &
                                     (Pool_1_Original(I)-Pool_1(I))**2*(Pnlty_Strength(I)*100.)
		end do

                ChiSquare_Penalty= ChiSquare_Penalty!*FLag_Pnlty !Flag Penalty read in from file:flag_status. 0=Off and 1=On

        RETURN 
        END FUNCTION ChiSquare_Penalty


!With definition of Beta_Penalty, it comes from ChiSquare term ((PW-PW_original)/Error)**2. 
!Beta is then is defined as the deriviative of this function. However, I write it as
!(PW_Orig - PW_current)*Strength*100. A strong value for Error is around .01 and a weak value around .1. 
!This means 100 to 10,000. By defining the Strength as Strength*100=1/Error**2, 
!A good constraining value is between 1 and 100 with 0 being OFF or no contribution to penalty.
	FUNCTION Beta_Penalty(Free_Param)
	double precision :: Beta_Penalty
	integer :: Free_Param

	  Beta_Penalty = (Pool_1_Original(Free_Param)-Pool_1(Free_Param))*abs(Pnlty_Strength(Free_Param)*100.)

	RETURN
	END FUNCTION Beta_Penalty


!See comment above Functoin Beta_penalty for explanation of penalty term to Alpha matrix. 
!This time Alpha is a second derivative of the parameter. So all the remains is an error(Strength) term.
	FUNCTION Alpha_Penalty(Free_Param)
        double precision :: Alpha_Penalty
 	integer :: Free_Param

	  Alpha_Penalty = abs(Pnlty_Strength(Free_Param)*100.)

	RETURN
	END FUNCTION Alpha_Penalty


	FUNCTION ChiSquare_Summation(bin_number)
	
		!*-----------------------------------------------------------------------
		!C       compute fitting data Y_FIT and total Chi-Square
		!*-----------------------------------------------------------------------
        	DOUBLE PRECISION :: ChiSquare_Summation!this function and hold the total summation of all observables.
        	DOUBLE PRECISION :: ERROR_Y_Original
		Integer :: bin_number, int_obs
		Integer :: data_tally !used to total datapoints of previous observables

		!takes care of the bins and observables that have 0 data so when it does not enter if statement below, eror_y is not 0 but Y is.
		Y=0.
		Y_Fit=0.
		Error_Y=.001
		ChiSummation = 0. !this summation is an individual observables total chivalue	
        	ChiSquare_Summation=0. !Total of all observables Chi Square value.
		data_tally=0	!holds a tally of number of datapoints from all previous obs.

		IterNum = IterNum + 1 !used for bug testing commented out

		!*-----------------------------------
		!Perform Chi Summation process
		Do int_obs=3,3+num_observables
		  if(dataForObs(int_obs) .ne. 0) then
		    DO I=1, dataForObs(int_obs) 
			!if(int_obs .eq. 5)print *,"I val is:",I+data_tally,I
			Y(I+data_tally)=total_array%energy_values(int_obs,bin_number)%energies_datapoints(I,2)
			Error_Y(I+data_tally)=total_array%energy_values(int_obs,bin_number)%energies_datapoints(I,3)
		   	Y_FIT(I+data_tally)=FUNCTION_FittingFunction(bin_number,int_obs,&
					COS(total_array%energy_values(int_obs,bin_number)%energies_datapoints(I,1)),I)
			!Chi of a given observable
			ChiSummation(int_obs)=ChiSummation(int_obs)+ChiSquare(Y(I+data_tally),Error_Y(I+data_tally),Y_FIT(I+data_tally))
			!if(iternum .eq. 2) then
			 !print *,"Y_Fit testing",int_obs,Y_FIT(I+data_tally),iternum, I+data_tally,Y(I+data_tally)
			 !print *,"ob values are:", p11rm(bin_num),p11im(bin_num),d13rm(bin_num),d13im(bin_num),	
			!end if
        	    END DO
		    data_tally = data_tally + dataForObs(int_obs)

		    !Total CHi Value of all observables
		    ChiSquare_Summation = ChiSquare_Summation + ChiSummation(int_obs)

		  end if
		End Do

        RETURN 
        END


!Fit called by main program.
	SUBROUTINE FIT(bin_number,fit_flag_bin_read,user_Reaction)
!*-----------------------------------------------------------------------
!C       fit data by partial wave analysis
!*-----------------------------------------------------------------------
	
	!passed variables
	Integer :: bin_number
	Integer :: fit_flag_bin_read,user_Reaction
	!local variables
        INTEGER :: Unit_FIT_a = 1
	character*40 :: date_run
        !if(bin_number .eq. 4) then
	  !print *, "Penalty Term ON/OFF:"
	  !print *, "0-OFF"
	  !print *, "1-ON"
          !read *, Flag_Pnlty
          !if((Flag_PNLTY .ne. 0) .and. (Flag_Pnlty .ne. 1)) Then
	  !  print *,"Flag Penalty set to off for invalid entry"
	    !Flag_Pnlty = 1
	  !end if 
        !end if

	!initialize a few variables that are aliases for a longer variable name
	!after loop, dataforobs set to 0 for observables that I don't want to use data for.
	Do int_obs = 3, 3+num_observables
	DataForObs(int_obs) = total_array%energy_num_points(int_obs,bin_number)
	End Do
	!Don't currently want to analyze SGT and Dx13(difference between Helicity 1/2 and 3/2 which is not independent) observables.
	DataForObs(3) = 0 
!        open(unit=35,file="hntestRIP.dat",status="unknown",position="append")

!       SET DEFAULT ERRORS
        Error_Default=0.025 

!*-----------------------------------------------------------------------
!*       1. read patial wave amplitudes: 'FIT_DataRead'
!*       2. fit new amplitudes by Chi-Square goodness-of-fit algorithm: 'ChiSquareFit'
!*       3. calculate observables with new amplitudes
!*	 4. Main Program is in charge of calling output.
!*-----------------------------------------------------------------------
	!print *,"fit_flag_bin_read in call fit is:", fit_flag_bin_read,bin_number
       if(bin_number .eq. 1) CALL FIT_DataRead(bin_number,fit_flag_bin_read,user_reaction)
       CALL ChiSquareFit(bin_number)
!*-----------------------------------
	
!	call system('mkdir -p ' //FILE_path(1:44))
!*-----------------------------------

        !WRITE(*, 101)
	!101     FORMAT(60('-'))
!	close(unit=35)
        RETURN
        END

	

	SUBROUTINE FIT_DataRead(bin_number,fit_flag_bin_read,user_reaction)
!*-----------------------------------------------------------------------
!C       read data for FIT from previous work by JOHN TULPAN
!*-----------------------------------------------------------------------
	implicit none

	!passed variables
	Integer :: bin_number 
	Integer :: fit_flag_bin_read
	Integer :: user_Reaction
	!Local Variables
	Integer :: turn_off_flags,j,k,last_bin!last_bin will determine how many bins of data exist in file.

	!Isospin 1/2
        character*4 :: Amplitudes(18)=(/'S11E','P11M','P13E','P13M','D13E','D13M','D15E','D15M',&
					'F15E','F15M','F17E','F17M','G17E','G17M','G19E','G19M','H19E','H19M'/)
	character*3 :: Amp_Filename(10)=(/'S11','P11','P13','D13','D15','F15','F17','G17','G19','H19'/)

	!Isospin 3/2
        character*4 :: Amplitudes3(18)=(/'S31E','P31M','P33E','P33M','D33E','D33M','D35E','D35M',&
					'F35E','F35M','F37E','F37M','G37E','G37M','G39E','G39M','H39E','H39M'/)
	character*3 :: Amp_Filename3(10)=(/'S31','P31','P33','D33','D35','F35','F37','G37','G39','H39'/)

	CHARACTER*1 :: FileName(2)=(/'E','M'/)!for patial wave amplitudes E and M monopoles
	DOUBLE PRECISION :: DUMP	
	double precision :: energy_temp,Erl_temp,Eim_temp,Mrl_temp,Mim_temp
!variable iostatus TURNS NEGATIVE IF RECEIVES AN END OF FILE NOTICE. PROGRAM CONTINUES AS IF PW VARYCODE WAS '0' OR '90'
	integer :: iostatus


	ED_Mod=0.0

	!Open files for read
	open(unit_pwa_input_initial,file=directory//subdirectory//"Input/initial_PWAs_small.dat",status="unknown")
	read(unit_pwa_input_initial,*) !read header line


	!obtain number of bins for output. Also occurs later, but needed here too.
	open(30,file=directory//subdirectory//"Input/flag_status_small.dat",status="unknown")
	  rewind(30)
	  read(30,*)!Header Line
	  read(30,*)!Header Line
	  read(30,*) end_bin
	close(unit=30)


	!read in Energy Dependent values and compute modulus for the starting value. File is of form 'Electri Multi header(except P11) data' THEN 'M multipole header(except S11) data'
	do j=1,10 !do loops runs through different files

	!open(unit=31,file=directory//subdirectory//"Input/NewInput/"//Amp_filename(j)&
	!	//"ED.dat",status="unknown")
	!rewind(unit=31)
	!read(31,*)!header for ED Multipole Dimensioned AMplitude File
	 do k=1,end_bin
	  
	  if(j .eq. 1) then
	    !read(31,*,end=1000) energy_temp,Erl_temp,Eim_temp
	    !ED_Mod(k,j) = dsqrt(Erl_temp**2+Eim_temp**2)
	  else if(j .eq. 2) then
	    !read(31,*,end=1000) energy_temp,Mrl_temp,Mim_temp
	    !ED_Mod(k,j+20) = dsqrt(Mrl_temp**2+Mim_temp**2)
	  else
	    !read(31,*,end=1000), energy_temp,Erl_temp,Eim_temp,Mrl_temp,Mim_temp
	    !ED_Mod(k,j) = dsqrt(Erl_temp**2+Eim_temp**2)
	    !ED_Mod(k,j+20) = dsqrt(Mrl_temp**2+Mim_temp**2)
	  end if
	 end do
	!close(unit=31)
	end do
1000	continue



        !read starting PWA values as either mod(a and b values same now or later) or Rl/Im
	!format for file is Header containing PW name followed by 1 line per bin for the PW values.

          !print *, "Inside fit_dataread case 0 with bin:", bin_number
	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,s01a(k),s01b(k)
	!print *,"S01 values aRE:", dump,s01a(k),s01b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)	!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p01a(k),p01b(k)
	!print *,"P01 values aRE:", dump,p01a(k),p01b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p03a(k),p03b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d03a(k),d03b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d05a(k),d05b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f05a(k),f05b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f07a(k),f07b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g07a(k),g07b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here3"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g09a(k),g09b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here4"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h09a(k),h09b(k)
	  end do
	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h011a(k),h011b(k)

	  end do

!ISOSPIN 1 amplitudes
	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,s11a(k),s11b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p11a(k),p11b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p13a(k),p13b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d13a(k),d13b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d15a(k),d15b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f15a(k),f15b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f17a(k),f17b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g17a(k),g17b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here3"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g19a(k),g19b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here4"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h19a(k),h19b(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h111a(k),h111b(k)
	  end do


	print *,"User reaction:",user_Reaction

	if(user_reaction .eq. 11 .or. user_reaction .eq. 12) then
!ISOSPIN 1 amplitudes
	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,s11c(k),s11d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p11c(k),p11d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,p13c(k),p13d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d13c(k),d13d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,d15c(k),d15d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f15c(k),f15d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,f17c(k),f17d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g17c(k),g17d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here3"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,g19c(k),g19d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	!print *,"Reads unit_pwa_input_intial ok here4"
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h19c(k),h19d(k)
	  end do

	  read(unit_pwa_input_initial,*,iostat=iostatus)!PW NAME HEADERLINE
	  do k=1,end_bin
	  read(unit_pwa_input_initial,*,iostat=iostatus), dump,h111c(k),h111d(k)
	  end do

	end if !IF reaction is number 11, it needs an extra set of amplitudes to read in.
	close(unit=unit_pwa_input_initial)

	do k=1,end_bin
	W(k) = total_array%central_energy(k)
	end do

	do j=70,87
	close(8000+j)
        end do


        RETURN
        END

!read initial flag values and makes sure that unphysical flag values are turned and kept off.  	
SUBROUTINE PartialWaveVaryCode(bin_number)
Implicit None
!*-----------------------------------------------------------------------
!*-----------------------------------------------------------------------
!flags(I) with I=2,12,21,31 are all non physical- they correspond to S11M real and imaginary and P11E real and imag

        INTEGER bin_number
	integer :: iostatus
	integer :: I,J
!*-----------------------------------------------------------------------
        DO I=1, NumberOfFlag
            Flag_1(I)=Flag_Off      !for isospin=1/2
        END DO
!*-----------------------------------

	open(unit_flag_status_file,file=directory//subdirectory//"Input/flag_status_small.dat",status="unknown")
	if(bin_number .eq. 1) then
	  rewind(unit_flag_status_file)
	  read(unit_flag_status_file,*)
	  read(unit_flag_status_file,*)
	  read(unit_flag_status_file,*) end_bin!,Flag_Pnlty
	  read(unit_flag_status_file,*) (Pnlty_Strength(I),I=1,10)
	  read(unit_flag_status_file,*) (Pnlty_Strength(I),I=21,30)

	!Make sure P11E and S11M have Strength=0  (off)
	        Pnlty_Strength(21) = 0.0 !S11Mr
		Pnlty_Strength( 2) = 0.0 !P11Er 

	!Set Real and Imaginary strength equal.
	  Do I=1,10
		Pnlty_Strength(I+10) = Pnlty_Strength(I)
		Pnlty_Strength(I+30) = Pnlty_Strength(I+20)
	  end do
                
	  !flags 2,12,21,31 should always be off. their intial wave amplitudes set to 0 on physical grounds.
	  do I=1,end_bin
	    read(unit_flag_status_file,*,iostat=iostatus) (flag_1_Original(I,J), J=1,NumberOfFlag)
	  end do

	  close(unit=unit_flag_status_file)
	end if

	  !Set flag status for bin.
	  do J=1, NumberOfFlag
	    flag_1(J) = Flag_1_Original(bin_number,J)
	  end do
	!Terminate program if an invalid combination of flag values is found. Does not check for flag values other than 1 or 2 which is done later.
	!Valid combinations are (Rl=1,IM=1), (Rl=1,Im=0), (Rl=0,Im=1), (Rl=2,Im=0 :This is mod phase fitting) 
	!F15a varied as mod while F15b varied as im or vice versa as example.
	!Also do not allow phase to vary as flag value of 2 since program fails if phase varies.
	  do J=1,10
            if((flag_1(J+10) .eq. 2)  .or. (flag_1(J+30) .eq. 2) .or.&
	       ((flag_1(J)   .eq. 2)  .and. (flag_1(J+10) .eq. 1)) .or.&
	       ((flag_1(J+20).eq. 2)  .and. (flag_1(J+30) .eq. 1))&
	      )Then
		print *, "Illegal Flag combo found for a Partial wave. Bin:",bin_number
	        print *, "Flag for Im/Phase must be 0 or 1. If Real/Mod flag is 2, Im/Phase must be 0."
	        stop "Fix input file for bin listed above to continue fit."
	     end if
	  end do

        RETURN
        END


        SUBROUTINE ChiSquareFit(bin_number)
	IMPLICIT NONE
!*-----------------------------------------------------------------------
!C       perform a least square fit using Chi-Square goodness-of-fit criterion
!*       ref.: Data Reduction and Error Analysis for the Physical Sciences,
!*             P. R. Bevington, D. K. Robinson, McGraw-Hill, 2003
!*-----------------------------------------------------------------------
!*-----------------------------------
        INTEGER :: Counter_Parameter, I, J, K, I1,kk !counter variables
	Integer :: Num_Constraints
	DOUBLE PRECISION :: A(NumberOfParameter_Maximum), B(NumberOfParameter_Maximum)
        DOUBLE PRECISION :: ChiSquare_A, ChiSquare_B
	DOUBLE PRECISION :: ChiSquare_Original, ChiSquare_swap 
	DOUBLE PRECISION :: Alpha(NumberOfParameter_Maximum, NumberOfParameter_Maximum), Beta(NumberOfParameter_Maximum)
!*-----------------------------------
        Double Precision, PARAMETER :: F_lambda_Maximum=1.0E+10
        Double Precision, PARAMETER :: F_lambda_Minimum=1.0E-8
        DOUBLE PRECISION :: F_lambda
        INTEGER :: Flag_F_lambda
!*-----------------------------------
        INTEGER :: DegreeOfFreedom
        DOUBLE PRECISION :: ChiSquare_Average, ChiSquare_Average_Original 
!*-----------------------------------
	!Termination variables -currently in use is Iteration_number, not time fit
	Integer :: Iteration_number
!        DOUBLE PRECISION :: Time_FIT_Initial, Time_FIT_Final, Time_FIT_Ellipse

	Integer :: bin_number
	Character*2 :: char_bin_number

	!print *,VCC
	!print *,s11pE(bin_number),p13pE(bin_number),d13pE(bin_number),"Phase testE"
	!print *,p11pM(bin_number),p13pM(bin_number),d13pM(bin_number),"Phase testM"
	!print *,VCC(1)*s11pE(bin_number),"Check for value of 1"
!*-----------------------------------------------------------------------------------------------------
	write(char_bin_number,'(i2)'),bin_number
	Open(ChiperPoint_file,file=Directory//subdirectory//"Output/observable_files/"//&
				"ChiPerPoint"//trim(adjustl(char_bin_number))//".dat",status='unknown')
	bin_num=bin_number
	write(ChiperPoint_file, '("Date MDY format:", i2,"/",i2,"/", i4)') datetoday(2), datetoday(1), datetoday(3)
	write(ChiperPoint_file,'("Y fit  Y actual  Y act error  Chi Value")')  
	write(ChiperPoint_file,'("Initial Chi Square Values")') 
	if(bin_number .eq. 1) then
	open(unit=931,file=Directory//subdirectory//"Output/ChiABvalues.dat",status="unknown")
		end if
!*-----------------------------------------------------------------------------------------------------
	 A=0 !Set array to 0

	!PartialWaveVaryCode determines which flags are on or off and new parameters assigned initial value.
        CALL PartialWaveVaryCode(bin_number) 
	!PWA's-> pool_1 variable
	Call Amplitude_Assign(1,bin_number)
	!changes initial parameter values to be constrained,especially flag 2 values to make sure they read the same.
	Call Protect(bin_number)
!	print *, num_constraints
	!if(bin_number .gt. 3) Num_Constraints = 2
!*-----------------------------------
	!set original error and amplitude values. If 0, give small non-0 (fix division by 0 error in program)
           DO I=1, NumberOfFlag
            if((abs(Pool_1(I)) .lt. .001) .and. (Flag_1(I).eq.1) )Pool_1(I) = .001

            Pool_1_Original(I)=Pool_1(I)
            Pool_1_Error(I)=Error_Default

           END DO

	!set original mod values.
           DO I=1,10

		!Electric Values
            SELECT CASE(Flag_1(I))
            CASE (0,1)
            Pool_1_Mod_Original(I)=dSqrt(Pool_1_Original(I)**2 + Pool_1_Original(I+10)**2)
            CASE (2)
            Pool_1_Mod_Original(I)=Pool_1_Original(I)
           END SELECT

		!Magnetic Values
            SELECT CASE(Flag_1(I+20))
            CASE (0,1)
            Pool_1_Mod_Original(I+20)=dSqrt(Pool_1_Original(I+20)**2 + Pool_1_Original(I+30)**2)
            CASE (2)
            Pool_1_Mod_Original(I+20)=Pool_1_Original(I+20)
           END SELECT
           END DO

	 !if(bin_number .eq. 4)print *,"Mod original values are:", Pool_1_Mod_Original
!*-----------------------------------------------------------------------
!C       Pool->A, read amplitudes into array A
!*-----------------------------------------------------------------------
        Counter_Parameter=0

        DO I1=1, NumberOfFlag    !Number of Flag = Number_Of_Flag_Maximum==40
            SELECT CASE(Flag_1(I1))
            CASE (0)
                CONTINUE
            CASE (1,2)
               Counter_Parameter=Counter_Parameter+1
                A(Counter_Parameter)=Pool_1(I1) !Read in amplitudes that are flagged as on.
             !print *, '1669 A(i): ', A(I1) 
           END SELECT
        END DO
 	!print *, "poolRIP:", pool_1
	Flag_last_ChiSummation = 1	
        NumberOfParameter=Counter_Parameter
	!print *,"values before first call:", a,"POOL:",Pool_1
	!include all observable datapoints expand 
	NumberOfData = 0
	do kk=3,3+num_observables
        NumberOfData= NumberOfData + dataForObs(kk)
	end do
	if(numberofData .gt. NumberOfStack_Maximum) then
	  print*,"NumberOfData has exceeded NumberOfStack_Maximum. (More data than array alloc for bin:", bin_number, ")"
	  print *, "NumberOfData in bin is:", NumberOfData
          stop "Program Terminates. Increase parameter NumberOfStack_Maximum to value higher than NumberOfData."
	end if
        DegreeOfFreedom=NumberOfData-NumberOfParameter !+ Num_Constraints
        Iteration_number=1
        ChiSquare_Original=ChiSquare_Summation(bin_number)
	Flag_last_ChiSummation = 0
!	print *,"FIRST CHISUM RIP:",ChiSquare_Original
	!print *, "Y is: " !, Y ---->Y is ok after test
	!print *, "Error Y is: " , Error_Y- is ok after test
        ChiSquare_Average_Original=ChiSquare_Original/DegreeOfFreedom
        IF ((DegreeOfFreedom.LE.0).OR.(NumberOfData.EQ.0)) THEN
            Flag_F_lambda=-1
	Else
!           Time_FIT_Initial=time()
           Flag_F_lambda=1
           F_lambda=0.001 !bevington suggests making this .001 at start of a fit
           FreeParameter=0 !sets array to 0
           Error_FreeParameter=0 !sets array to 0
	END IF
!*-----------------------------------------------------------------------
!*End of intialization of Initial Variables
!*-----------------------------------------------------------------------
!
!
1       continue
	!print *, "goes here for chi b less than chi a"
!
!	Amp_Assigns does the following:
!	A->Free_Parameter->Pool_1
!	Error_Free_Parameter->Pool_1_Error
!	Pool_1 => PWA's 
!	For flag=2 values, it also adjusts fixed parameters that need to match varied parameters (a and b values F15a and F15b etc.)
!	Amp_Assign(5,...) is different then Amp_Assign(4,...) because it sets FreeParameter to be array A() first.(resets to previous values)
	Call Amplitude_Assign(5,bin_number,A)
	Call Amplitude_Assign(2,bin_number)
!        
!	Compute ChiSquare+penalty terms.
        ChiSquare_A=ChiSquare_Summation(bin_number)+ChiSquare_Penalty()  !line 16/17 in Bevington book code page 229	
	write(931,'(a2,f9.2,2x,f9.2,3x,i3)') "A:", ChiSquare_A,chisquare_penalty(),bin_number

!*-----------------------------------------------------------------------
!*       iteration process to get the parameter matrices
!         only purpose of do while loop is to skip lines IF Diff between chi_A and chi_B is small
!*-----------------------------------------------------------------------
	Do While(Flag_F_Lambda .Gt. 0)

!*------------------------------------------------------------------------
!COMPUTE THE DERIVATIVE OF FITTING FUNCTIONS 'assigns dF for ALL points/parameters' (alters PWA's) 
         CALL FittingFunctionDerivative(bin_number)

!*-----------------------------------
!*       COMPUTE 'Alpha' and 'Beta', see Bevington[1969] (11-16), (11-25)
            Beta=0. !vector to 0
            Alpha=0.!matrix to 0
!*-----------------------------------
!*       RECOVER THE ORIGINAL AMPLITUDES
!*       Pool_Original => s01r(bin_number)..., read amplitudes to parameter set
!*-----------------------------------------------------------------------
	 Call Amplitude_Assign(3,bin_number)
         !CALL PartialWaveVaryCode(bin_number)   !currently not needed because changes are made below.

	!calculation of Beta penalty contribution.
         Counter_Parameter=0
         DO I1=1, NumberOfFlag    !see Amplitude.h, for Isospin=1

            SELECT CASE(Flag_1(I1))
            CASE (0)
                CONTINUE
            CASE (1)
                Counter_Parameter=Counter_Parameter+1

                Beta(Counter_Parameter)=Beta_Penalty(I1)
		!c   Beta(Counter_Parameter)  = (Pool_1_Original(I1) - Pool_1(I1))/(Pool_1_Error(I1)**2)
                Alpha(Counter_Parameter,Counter_Parameter)=Alpha_Penalty(I1)
		!c   Alpha(Counter_Parameter,Counter_Parameter)  =  1./Pool_1_Error(Counter_Parameter)**2
            END SELECT
         END DO  

!*       LINEARIZATION OF FITTING FUNCTION
         DO I=1, NumberOfParameter
            IF (ERROR_Y(I).LT..0001) THEN
                PRINT *, 'Found an error that is less than .0001(i.e. 0), '
                PRINT *, '1:Check original database to fix bad error point OR'
		Print *, "2:If errors are less than .0001, reduce lower limit on error in code above"
                RETURN
            ELSE
                DO J=1, NumberOfData
		    Beta(I)=Beta(I)+dF(J, I)*(Y(J)-Y_FIT(J))/ERROR_Y(J)**2
                END DO

                DO K=1, I
                  DO J=1, NumberOfData
			Alpha(I, K)=Alpha(I, K)+dF(J, I)*dF(J, K)/ERROR_Y(J)**2
                  END DO
                    Alpha(K, I)=Alpha(I, K)
		!print *, "K,I,Alpha matrix(K,I)", K,I,Alpha(K,I)
                END DO
            END IF
         END DO


!*-----------------
2        continue            !goto statement below-goes here if Chi_A<Chi_B
!	print *, "goes here for chi a less than chi b"
!*-----------------------------------------------------------------------
!C       B->Pool, read comparable amplitudes from B to pools
!*-----------------------------------------------------------------------
!c       COMPUTE CURVATURE MATRIX 'ARRAY' from 'Alpha'
         DO I=1, NumberOfParameter
            DO J=1, I
                ARRAY(I, J)=Alpha(I, J)/DSQRT(Alpha(I, I)*Alpha(J, J))
                ARRAY(J, I)=ARRAY(I, J)
            END DO
            ARRAY(I, I)=1.+F_lambda
         END DO
!*-----------------

!*       COMPUTE THE INVERSE MATRIX OF 'ARRAY'
         CALL MatrixInverse


!*-----------------
!*       COMPUTE NEW PARAMETER ESTIMATES
         DO I=1, NumberOfParameter
            B(I)=A(I)
             !print *, '1669 b(i): ', B(I)                       
            DO J=1, NumberOfParameter
	     !B(I) is new parameter vector(what I want to know), Array: inverse matrix, and Beta: Y matrix in Chi Squared. 
	     B(I)=B(I)+ARRAY(I, J)*Beta(J)/DSQRT(Alpha(I, I)*Alpha(J, J))
            END DO
         END DO
!*----------------------------------            
         DO I=1, NumberOfParameter
            FreeParameter(I)=B(I)
         END DO

	!Update PWA's with B array values. These may or may not be better than array A's values.
	!Free_Parameter->Pool_1
	!For Flag=2 parameters, also adjusts fixed parameters that need to match varied parameters. (F15b=F15a etc.)
	!Pool_1 => PWA(bin_number) ...
	 Call Amplitude_Assign(4,bin_number)
         Call Amplitude_Assign(2,bin_number)

!
!	 Physical constraints on data.
!	 Call Protect_GPEP(f_lambda) 

!
!        Determine new Chi Square value.
         ChiSquare_B=ChiSquare_Summation(bin_number)+ChiSquare_Penalty()
	write(931,'(a2,f9.2,2x,f9.2,3x,i3)') "B:", ChiSquare_B,ChiSquare_Penalty(),bin_number
!	print *, "chi a and b are:", chisquare_a, chisquare_b
!*-----------------------------------------------------------------------
!*       Test for convergence. If change of ChiSquare is less than .5 terminate
!*-----------------------------------------------------------------------

	if(abs(ChiSquare_A -ChiSquare_B) .lt. .5) then
	ChiSquare_A = Min(ChiSquare_A,ChiSquare_B)
	ChiSquare_B = Chisquare_A
	exit
	end if

         IF (ChiSquare_A.LT.ChiSquare_B) THEN
!         IF(CHISQUARE_A .LT. CHISQUARE_C) CHISQUARE_C=CHISQUARE_A
!         if((chisquare_c .lt. chisquare_a) .and. (chisquare_c .lt. chisquare_b)) exit
         	F_lambda=F_lambda*10. 

	!If maximum iterations is met AND A is less than B, then use A values instead of B values as PW values before final iteration.
              IF ((Iteration_number .eq. 30) .or. (F_Lambda .GT. F_Lambda_Maximum)) Then
	!	Amp_Assigns does the following:
	!	A->Free_Parameter->Pool_1
	!	Error_Free_Parameter->Pool_1_Error
	!	Pool_1 => PWA's 
	!	For flag=2 values, it also adjusts fixed parameters that need to match varied parameters (a and b values F15a and F15b etc.)
	!	Amp_Assign(5,...) is different then Amp_Assign(4,...) because it sets FreeParameter to be array A() first.(resets to previous values)
		Call Amplitude_Assign(5,bin_number,A)
		Call Amplitude_Assign(2,bin_number)

		EXIT

	      End if

	      Iteration_number=Iteration_number+1
              GOTO 2

         ELSE

		F_lambda=F_lambda/10.

            DO I=1, NumberOfParameter
                A(I)=B(I)
                IF (Alpha(I, I).NE.0.) THEN
                    !write(*,*) Error_FreeParameter(i), "before",i
                    Error_FreeParameter(I)=DSQRT(ARRAY(I, I)/Alpha(I, I))
                    !write(*,*) Error_FreeParameter(i), "after",i
                END IF
            END DO
 
            IF ((Iteration_number .eq. 30) .or. (F_Lambda .LT. F_Lambda_Minimum)) Exit

            Iteration_number=Iteration_number+1
            GOTO 1
         END IF

	End Do	!ENd DO while loop





!*--------------------------------------------------------------------------------------------------------------
!*
!*       Entering final out of the ChiSquare Fit
!*	   -Computes Final errors and outputs PW fit values.
!*
!*--------------------------------------------------------------------------------------------------------------


	if (bin_number .eq. 1 ) then	
	open(unit=993,file=directory//subdirectory//"Output/Final_chi_results.dat",status='unknown')
	    write(993, '(a27,a6)'),"REACTION being analyzed is:", subdirectory
	else
        open(unit=993,file=directory//subdirectory//&
					"Output/Final_chi_results.dat",status='unknown',position='append')
	end if


	    write(993, '(90("="),/,90("="))')
	    write(993, '("Date MDY format:", i2,"/",i2,"/", i4)') datetoday(2), datetoday(1), datetoday(3)
	    write(993, '(a11,f6.1,5x,a3,f4.1,a5)'), "BIN ENERGY:",total_array%central_energy(bin_number),"+/-",&
					                     bin_width(bin_number)/2., "(MeV)"
	    write(993, '("CURRENT BIN:",i2,/)') bin_number
            WRITE(993, '(t3,a3,t8,a3,t14,a4,t21,a2,t27,a8,t36,a5,t44,a4,t50,a4,t57,a8,t67,a2,t70,a5,t80,a3,t86,a4,t92,a5,t99,a3)'),&
			'PW#',     'VC*',        'RL**',     'RL',    'ERR(+/-)',&
                        '||VC*',   'IM**', 'IM',  'ERR(+/-)' , "||", "M_ED?", "M**", 'MOD^', "M_err", 'PHS'
!350	Format (t3,a4,t8,a2,t11,f7.2,t18,f7.2,t26,f6.2,t38,a2,t42,f7.2,t53,f7.2,t62,f6.2,t73,a2,t76,f6.2,t83,f6.2,t90,f6.2,t97,f5.2)
            ChiSquare_Average=1
!
!	 call Protect_GPEP(Num_Constraints)
!
            IF ((DegreeOfFreedom.GT.0).AND.(NumberOfData.NE.0)) &
                ChiSquare_Average=(ChiSquare_Summation(bin_number)+ChiSquare_Penalty())/DegreeOfFreedom




	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!
	!           CALCULATE FINAL ERRORS FOR EACH AMPLITUDE
	!           
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
        NumberOfParameter_Original=NumberOfParameter
        NumberOfParameter=0

          DO I=1, NumberOfFlag

            IF ((Flag_1(I).EQ.Flag_On) .or. (Flag_1(I).EQ.Flag_MP_On)) THEN
               NumberOfParameter=NumberOfParameter+1
                FreeParameter(NumberOfParameter)=Pool_1(I)
            END IF
        END DO
	!*-----------------

	!COMPUTE THE DERIVATIVE OF FITTING FUNCTIONS 'dF'

        CALL FittingFunctionDerivative(bin_number)

	!Make Beta and Alpha vectors 0
        Beta=0.
        Alpha=0.
 
	!*-----------------
	!c       LINEARIZATION OF FITTING FUNCTION
         Counter_Parameter=0
         DO I1=1, NumberOfFlag    !see Amplitude.h, for Isospin=1

            SELECT CASE(Flag_1(I1))
            CASE (0)
                CONTINUE
            CASE (1)
                Counter_Parameter=Counter_Parameter+1

                Beta(Counter_Parameter)=Beta_Penalty(I1)
		!c   Beta(Counter_Parameter)  = (Pool_1_Original(I1) - Pool_1(I1))/(Pool_1_Error(I1)**2)
                Alpha(Counter_Parameter,Counter_Parameter)=Alpha_Penalty(I1)
		!c   Alpha(Counter_Parameter,Counter_Parameter)  =  1./Pool_1_Error(Counter_Parameter)**2
            END SELECT
         END DO 

        DO I=1, NumberOfParameter
            DO J=1, NumberOfData
		Beta(I)=Beta(I)+dF(J, I)*(Y(J)-Y_FIT(J))/ERROR_Y(J)**2
            END DO

            DO K=1, I
                DO J=1, NumberOfData
		    Alpha(I, K)=Alpha(I, K)+dF(J, I)*dF(J, K)/ERROR_Y(J)**2
                END DO
                Alpha(K, I)=Alpha(I, K)
            END DO
        END DO

	!*-----------------
	!c       COMPUTE CURVATURE MATRIX 'ARRAY' from 'Alpha'
        DO I=1, NumberOfParameter
            DO J=1, I
                IF ((Alpha(I, I)*Alpha(J, J)).NE.0) THEN
	            ARRAY(I, J)=Alpha(I, J)/DSQRT(Alpha(I, I)*Alpha(J, J))
                    ARRAY(J, I)=ARRAY(I, J)
                END IF
            END DO

            ARRAY(I, I)=1.0
        END DO

	!*-----------------
	!c       COMPUTE THE INVERSE MATRIX OF 'ARRAY' and set equal to 'ARRAY'

	CALL MatrixInverse

	!*-----------------------------------
	Do I=1,NumberOfFlag
	  Error_FreeParameter(I) = 0.
	End DO

        DO I=1, NumberOfParameter
            IF (NumberOfParameter.GE.NumberOfData) THEN
                Error_FreeParameter(I)=Error_Default
            ELSE IF (Alpha(I, I).NE.0) THEN
	        Error_FreeParameter(I)=DSQRT(ARRAY(I, I)/Alpha(I, I))
                !print *, "final error is:", error_FreeParameter(I), i
            END IF
        END DO

	!*-----------------
        Counter_Parameter=NumberOfParameter

        DO I1=NumberOfFlag, 1, -1 
            SELECT CASE(Flag_1(I1))
            CASE (0)
                CONTINUE
            CASE (1,2)
                Pool_1_Error(I1)=Error_FreeParameter(Counter_Parameter)
                Counter_Parameter=Counter_Parameter-1
            END SELECT
        END DO

	!*-----------------
            DO I=1, NumberOfFlag
                ChiSquare_Average=MAX(ChiSquare_Average, 1.0)
                IF (Pool_1_Error(I).EQ.Error_Default) then
		  Pool_1_Error(I)=0.0
		else
	          Pool_1_Error(I)=Pool_1_Error(I)!*SQRT(ChiSquare_Summation(bin_number))
                END IF
            END DO

	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!
	!           END CALCULATE FINAL ERRORS FOR EACH AMPLITUDE
	!           
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------
	!*------------------------------------------------------------------------------------------------------


	!Pool_1->PWA's  Called before mod_phase to prevent bugs in mod/phase amplitudes in Helicity Definitions.
	call Amplitude_assign(2,bin_number)

	!Changes Pool_1 into a completely Rl Im value holder(takes all mod values and changes them to RL IM)
	!Also determines the Pool_1_Modphase values.
	!Call Mod_Phase(bin_number)
	!Pool_Mod->PW_mod
        !call Amplitude_assign(6,bin_number) 
!---------------------------------------------------------------------------------------------------------------------
!Case 0 is flag off 
!Case 1,2 is flag on as either rl/im OR MOD phaes
!Outside case is for the real part , Inside cases then decide if the imaginary part is on or not.
!ELectric Multipoles
!Units are in mfm.
!---------------------------------------------------------------------------------------------------------------------

		!print *,"Chi penalty right be write statements:", chisquare_penalty()
	    DO I1=1, NumberOfFlag/4
		if(I1.ne.2) then
		if((abs(Pool_1_Original(I1)) .ge. .001) .or. (abs(Pool_1_Original(I1+10)) .ge. .001)) then
                SELECT CASE(Flag_1_Original(bin_number,I1))

                CASE (0)
                        SELECT CASE(Flag_1_Original(bin_number,I1+10))
                        CASE (0)
				WRITE(993,350),&
     				PartialWave_1(I1),&
     				VaryCode(Flag_1(I1)),Pool_1_Original(I1), Pool_1(I1), Pool_1_Error(I1), "||",&
     				VaryCode(Flag_1(I1+10)),Pool_1_Original(I1+10), Pool_1(I1+10), 0.00, "||",&
				ED_Mod(bin_number,I1),Pool_1_Mod_Original(I1),Pool_1_Modphase(I1),Pool_1_MP_Error(I1),Pool_1_modphase(I1+10)
                        CASE (1,2)
				WRITE(993,360),&
     				&PartialWave_1(I1),&
     				VaryCode(Flag_1(I1)),Pool_1_Original(I1), Pool_1(I1), 0.0, "||",&
     				VaryCode(Flag_1(I1+10)),Pool_1_Original(I1+10), Pool_1(I1+10), Pool_1_Error(I1+10), "||",&
				ED_Mod(bin_number,I1),Pool_1_Mod_Original(I1),Pool_1_Modphase(I1),Pool_1_MP_Error(I1),Pool_1_modphase(I1+10)

                        END SELECT
                CASE (1,2)
                        SELECT CASE(Flag_1_Original(bin_number,I1+10))
                        CASE (0)
				WRITE(993,350),&
     				PartialWave_1(I1),&
     				VaryCode(Flag_1(I1)),Pool_1_Original(I1), Pool_1(I1), Pool_1_Error(I1), "||",&
     				VaryCode(Flag_1(I1+10)),Pool_1_Original(I1+10), Pool_1(I1+10), Pool_1_Error(I1+10), "||",&
				ED_Mod(bin_number,I1),Pool_1_Mod_Original(I1),Pool_1_Modphase(I1),Pool_1_MP_Error(I1),Pool_1_modphase(I1+10)

                        CASE (1,2)
				WRITE(993,360),&
     				PartialWave_1(I1),&
     				VaryCode(Flag_1(I1)),Pool_1_Original(I1), Pool_1(I1), Pool_1_Error(I1), "||",&
     				VaryCode(Flag_1(I1+10)),Pool_1_Original(I1+10), Pool_1(I1+10), Pool_1_Error(I1+10), "||",&
				ED_Mod(bin_number,I1),Pool_1_Mod_Original(I1),Pool_1_Modphase(I1),Pool_1_MP_Error(I1),Pool_1_modphase(I1+10)

                        END SELECT

                END SELECT
		end if
		END IF

!Magnetic multipoles---------------------------------------------------------------------------------------------------------------------
		if(I1.ne.1) then
		if((abs(Pool_1_Original(I1+20)).ge. .001) .or. (abs(Pool_1_Original(I1+30)).ge. .001)) then
                SELECT CASE(Flag_1_Original(bin_number,I1+20))
                CASE (0)
                        SELECT CASE(Flag_1_Original(bin_number,I1+10))
                        CASE (0)
				WRITE(993,350)&
     				PartialWave_1(I1+10),&
     				VaryCode(Flag_1(I1+20)),Pool_1_Original(I1+20), Pool_1(I1+20), Pool_1_Error(I1+20), "||",&
     				VaryCode(Flag_1(I1+30)),Pool_1_Original(I1+30), Pool_1(I1+30), 0.0, "||",&
				ED_Mod(bin_number,I1+20),Pool_1_Mod_Original(I1+20),&
					Pool_1_Modphase(I1+20),Pool_1_MP_Error(I1+20),Pool_1_modphase(I1+30)
                        CASE (1,2)
				WRITE(993,360)&
     				&PartialWave_1(I1+10),&
     				VaryCode(Flag_1(I1+20)),Pool_1_Original(I1+20), Pool_1(I1+20), 0.0, "||",&
     				VaryCode(Flag_1(I1+30)),Pool_1_Original(I1+30), Pool_1(I1+30), Pool_1_Error(I1+30), "||",&
				ED_Mod(bin_number,I1+20),Pool_1_Mod_Original(I1+20),&
				Pool_1_Modphase(I1+20),Pool_1_MP_Error(I1+20),Pool_1_modphase(I1+30)

                        END SELECT
                CASE (1,2)
                        SELECT CASE(Flag_1_Original(bin_number,I1+30))
                        CASE (0)
				WRITE(993,350)&
     				PartialWave_1(I1+10),&
     				VaryCode(Flag_1(I1+20)),Pool_1_Original(I1+20), Pool_1(I1+20), Pool_1_Error(I1+20), "||",&
     				VaryCode(Flag_1(I1+30)),Pool_1_Original(I1+30), Pool_1(I1+30), Pool_1_Error(I1+30), "||",&
				ED_Mod(bin_number,I1+20),Pool_1_Mod_Original(I1+20),&
				Pool_1_Modphase(I1+20),Pool_1_MP_Error(I1+20),Pool_1_modphase(I1+30)

                        CASE (1,2)
				WRITE(993,360)&
     				PartialWave_1(I1+10),&
     				VaryCode(Flag_1(I1+20)),Pool_1_Original(I1+20), Pool_1(I1+20), Pool_1_Error(I1+20), "||",&
     				VaryCode(Flag_1(I1+30)),Pool_1_Original(I1+30), Pool_1(I1+30), Pool_1_Error(I1+30), "||",&
				ED_Mod(bin_number,I1+20),Pool_1_Mod_Original(I1+20),&
				Pool_1_Modphase(I1+20),Pool_1_MP_Error(I1+20),Pool_1_modphase(I1+30)

                        END SELECT
                END SELECT
		end if
		end if
            END DO

!---------------------------------------------------------------------------------------------------------------------
!Case 0 is flag off 
!Case 1 is flag on
!Outside case: real part , Inside cases: imaginary part.
!Magnetic Multipoles
!---------------------------------------------------------------------------------------------------------------------

	!print *,"HD3 and other 3 helicities are+bin_num:", HD_3(), HN_3(), Hsp_3(), Hsa_3(), bin_num
	    write(993, "('PWA units: (mfm)',/)")
            write(993, "(/,'*: 0=fixed, 90=Varied as Rl/Im, 45=Varied as Mod/Phs')")  !!! **: from John Tulpan'
	    write(993, "('**: Starting amplitude values')")
	    write(993, "('^: For VC=45, Pool original still is given in Rl/Im format')")
	    write(993, "('^: For VC=90, Mod is always Positive')") 
 	    write(993, "('?: Care Should be taken in determining if this is predicted or actual.')") 
 	    write(993, "('?: Predicted=not Analyzed in ED code??possible Bug; Actual=Best Current Value.',/,/)")
	!*-----------------------------------
            write(993, '(a25,i2)'), 'Number of parameters is: ', NumberOfParameter
            !write(993, '(a25,i2)'), 'Number of Constraints is:', Num_Constraints
	    write(993, '(a24,i4)'), 'Degree of freedom is:    ', DegreeOfFreedom
	    write(993, '("Fit finished Running in:",i2, " Iterations")'), Iteration_number
            IF ((NumberOfParameter .GE. NumberOfData) .OR. (NumberOfData .EQ. 0)) THEN
                write(993, '(a28,t35,f10.1)'), 'Total CHI-SQUARE(Sum+Penal):', ChiSquare_Summation(bin_number)+ChiSquare_Penalty()
                write(993, '(a19,t35,f10.1)'), 'Penalty CHI-SQUARE:', ChiSquare_Penalty()
            ELSE
                !print *, "Final call to penalty:"
                write(993, '(a29,t35,f10.2)'), 'Initial CHI-SQUARE(Sum only):',ChiSquare_Original
                write(993, '(a28,t35,f10.2)'), 'Final CHI-SQUARE(Sum+Pnlty):',(ChiSquare_Average)*DegreeOfFreedom
                write(993, '(a21,t35,f10.2,/,40("-"))'), 'Penalty Contribution:', ChiSquare_Penalty()
                write(993, '(a32,t35,f10.2)'),  '             Initial Chi per df:',ChiSquare_Average_Original
                write(993, '(a32,t35,f10.2)'),  'Final Chi per df  (Sum Only)/dF:',ChiSquare_Summation(bin_number)/DegreeOfFreedom
                write(993, '(a32,t35,f10.2)'),  'Total Chi per df (Sum+Pnlty)/dF:',ChiSquare_Average
               				                 

		write(993,'(/,a35)'),"Avg Chi Per Point from Observables:"
		!Prints to file an individual observables chi value and number of datapoints.
		Do kk=3,3+num_observables
		  if (dataForObs(kk) .ne. 0) then
		    write(993,370), observable_files(kk),ChiSummation(kk)/dataForObs(kk),"Number Of Data:",dataForObs(kk)
		    GLOBAL_ChiSum(kk) = ChiSummation(kk) + GLOBAL_ChiSum(kk)
		  !else
		    !write(993,380), observable_files(kk),                                "Number Of Data:",dataForObs(kk)
		  end if
		End DO

            END IF

!(t3,a3,t10,a3,t16,a4,t23,a2,t29,a8,t42,a3,t48,a4,t60,a2,t72,a8,t83,a3,t88,a3,t93,a3)')
350	Format(t3,a4,t8,a2,t11,f7.2,t18,f7.2,t26,f6.2,t36,a2,t38,a2,&
		&t40,f7.2,t49,f6.2,t56,f6.2,t67,a2,t70,f6.2,t78,f6.2,t84,f6.2,t91,f6.2,t98,f5.2)
360	Format(t3,a4,t8,a2,t11,f7.2,t18,f7.2,t26,f6.2,t36,a2,t38,a2,&
		&t40,f7.2,t49,f6.2,t56,f6.2,t67,a2,t70,f6.2,t78,f6.2,t84,f6.2,t91,f6.2,t98,f5.2)
370	Format(a7,T10,F8.2,5x,a15,i3)
380	Format(a7,T23,a15,i3)

	close(unit=993)

	!Compute the final chi square per point values and put into file.
	Flag_last_ChiSummation = 1
	write(ChiperPoint_file,'("Final Chi Square Values")')
	ChiSquare_Average = ChiSquare_Summation(bin_number)+ChiSquare_Penalty()
	write(*, '(20("-"),"Bin ",i2, " Of ", i2, " Completed",20("-"))') bin_number, end_bin

	close(unit=ChiperPoint_File)
        RETURN
        END











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE RESTRICTS VALUES OF PARAMETERS FOR GPEP REACTION
	SUBROUTINE PROTECT(bin_protect)
	!passed variables
	!integer, intent(out) :: Num_Constraints
	!local variables
	integer :: i_protect,bin_protect

	!Protect electric Multipoles - make them modulus values if different and flag=2
	do i_protect=1,10
	if(Flag_1(i_protect) .eq. 2) then
	  if(pool_1(i_protect).ne. pool_1(i_protect+10)) then
		Pool_1(i_protect) = dsqrt(pool_1(i_protect+10)**2 + Pool_1(i_protect)**2)
		Pool_1(i_protect+10) = Pool_1(i_protect)
	  end if
	end if
	end do

	!Protect Magnetic Multipoles - make them modulus values if different and flag=2
	do i_protect=21,30
	if(Flag_1(i_protect) .eq. 2) then
	  if(pool_1(i_protect).ne. pool_1(i_protect+10)) then
		Pool_1(i_protect) = dsqrt(pool_1(i_protect+10)**2 + Pool_1(i_protect)**2)
		Pool_1(i_protect+10) = Pool_1(i_protect)
	  end if
	end if
	end do

	!Call amplitude assign so that PW's get assigned correctly as well before first summation call.
	Call Amplitude_Assign(2,bin_protect)
	return
	END
!

!
!	Unused currently.
	SUBROUTINE Interpolate(E)
!*-----------------------------------------------------------------------
!C       interpolate patial-wave amplitudes at C.M. energy E
!*-----------------------------------------------------------------------
!C       for wide-width resonances whose partial-wave amplitude varies smoothly, energy-independent
!*		Taylor expansion: X(E)=X(E0)+X'(E0)*(E-E0) 
!*		where X'(E0)=(X(Emax)-X(Emin))/(Emax-Emin)
!C       for narrow-width resonances whose partial-wave amplitude varies greatly, energy-dependent
!*		Linear interpolation: (Y-Y0)/(Y1-Y0)=(X-X0)/(X1-X0)
!*		Generic interpolation via Lagrange polynomial
!*			ref.: http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
!*-----------------------------------------------------------------------
        DOUBLE PRECISION E
        DOUBLE PRECISION dE, E0
        INTEGER bin_number


!C       Linear Interpolation: (Y-Y0)/(Y1-Y0)=(X-X0)/(X1-X0)
!*-----------------------------------------------------------------------
!C       interpolate partial wave amplitudes for wide-width resonance, energy-independent
!*-----------------------------------------------------------------------
        !E0=total_array%central_energy(bin_number)
        !dE=1.*(W(bin_number+1)-W(bin_number-1))

        RETURN   
        END


	SUBROUTINE MatrixInverse
!*-----------------------------------------------------------------------
!C       invert a real symmetric matrix ARRAY
!*-----------------------------------------------------------------------

        INTEGER NumberOfRow(NumberOfParameter_Maximum)
        INTEGER NumberOfColumn(NumberOfParameter_Maximum)
 
        DOUBLE PRECISION A_Maximum, Stack
!*-----------------------------------------------------------------------
        DO K=1, NumberOfParameter
            A_Maximum=0.
        
!c           LOCATE THE LARGEST ELEMENT A_Maximum OF MATRIX ARRAY
            DO I=K, NumberOfParameter
                DO J=K, NumberOfParameter
                    IF(A_Maximum-ARRAY(I, J).LE.0.) THEN
                        A_Maximum=ARRAY(I, J)
                        NumberOfRow(K)=I
                        NumberOfColumn(K)=J
                    END IF
                END DO
            END DO
            IF(A_Maximum.EQ.0.) RETURN
  
!c           INTERCHANGE ROWS AND COLUMNS TO PUT A_Maximum INTO ARRAY(K, K)
            IF(NumberOfRow(K).GT.K) THEN
                DO J=1, NumberOfParameter
                    Stack=ARRAY(K, J)
                    ARRAY(K, J)=ARRAY(NumberOfRow(K), J)
                    ARRAY(NumberOfRow(K), J)=-Stack
                END DO
            END IF
      
            IF(NumberOfColumn(K).GT.K) THEN
                DO I=1, NumberOfParameter
                    Stack=ARRAY(I, K)
                    ARRAY(I, K)=ARRAY(I, NumberOfColumn(K))
                    ARRAY(I, NumberOfColumn(K))=-Stack
                END DO
            END IF

            DO I=1, NumberOfParameter
                IF(I.NE.K) THEN
                    ARRAY(I, K)=-ARRAY(I, K)/A_Maximum
                END IF
            END DO

            DO I=1, NumberOfParameter
                DO J=1, NumberOfParameter
                    IF((I.NE.K).AND.(J.NE.K)) THEN
                        ARRAY(I, J)=ARRAY(I, J)+ARRAY(I, K)*ARRAY(K, J)
                    END IF
                END DO
            END DO

            DO J=1, NumberOfParameter
                IF(J.NE.K) THEN
                    ARRAY(K, J)=ARRAY(K, J)/A_Maximum
                ELSE
                    ARRAY(K, J)=1./A_Maximum
                END IF
            END DO
        END DO
                                       
!c       RESTORE ORDERING OF MATRIX
        DO K=NumberOfParameter, 1, -1
            IF(NumberOfRow(K).GT.K) THEN
                DO I=1, NumberOfParameter
                    Stack=ARRAY(I, K)
                    ARRAY(I, K)=-ARRAY(I, NumberOfRow(K))
                    ARRAY(I, NumberOfRow(K))=Stack
                END DO
            END IF
    
            IF(NumberOfColumn(K).GT.K) THEN
                DO J=1, NumberOfParameter
                    Stack=ARRAY(K, J)
                    ARRAY(K, J)=-ARRAY(NumberOfColumn(K), J)
                    ARRAY(NumberOfColumn(K), J)=Stack
                END DO
            END IF
        END DO
      
        RETURN
        END




!Determines the value of the fitting function at a value slightly higher than the current PWA value and at a value slightly lower
!then it computes a slope based on the parameter. Each free parameter is treated as having a different slope.
!On a given call it will go through all observables that have data, all Data for that observable, AND all varied parameters.
	SUBROUTINE FittingFunctionDerivative(bin_number)
!*-----------------------------------------------------------------------
!C       calculate the derivatives of the fitting function with respect to free parameters
!*-----------------------------------------------------------------------

        DOUBLE PRECISION FreeParameter_original(NumberOfParameter)
        DOUBLE PRECISION FreeParameter_plus(NumberOfParameter)
        DOUBLE PRECISION FreeParameter_minus(NumberOfParameter)

        DOUBLE PRECISION Y_FIT_plus(NumberOfStack_Maximum)
        DOUBLE PRECISION Y_FIT_minus(NumberOfStack_Maximum)
	!epsilon needs to be small enough that a change in the smaller parameters does not produce a drastically different solution.
        Double Precision,PARAMETER :: epsilon=.002  !Parameter increment change from .01 to .002 Oct 2015. Small partial waves are on order of .1 or less
	integer bin_number, local_intObservable,NumberOfData_Obs
	integer :: data_tally = 0
!*-----------------------------------------------------------------------

          DO I=1, NumberOfParameter     !initialization by default
            FreeParameter_original(I)=FreeParameter(I)
            FreeParameter_plus(I)=0.
            FreeParameter_minus(I)=0.
	  End do

	data_tally = 0

	!loop over all observables that have data, but not SGT observables.
	Do local_intObservable = 3,3+num_observables
	If(dataForObs(local_intObservable).ne.0) Then
	  NumberOfData_Obs=dataForObs(local_intObservable) !added BCH 
	!print *, "df...dataforOBS val and loca_int val:", dataforobs(local_intObservable), local_intObservable
	  !Determine the deterivative of all parameters that are varied
	  !Determine the amplitude at a central value, then vary by either + or - epsilon/2.
          DO I=1, NumberOfParameter
          IF (FreeParameter(I).EQ.0.) THEN !avoid division by 0 in dF = ... line.
                CONTINUE
            ELSE

	        CALL FreeParameter_Amplitude(bin_number,I, FreeParameter_original(I))
                     FreeParameter_plus(I)=FreeParameter(I)*(1.+epsilon/2.)
                CALL FreeParameter_Amplitude(bin_number,I, FreeParameter_plus(I))

		!because FreeParameter_Amplitude call changes the values of Ith Parameter, Y plus and Y minus are different
                DO J=1, NumberOfData_Obs
                    Y_FIT_plus(J) =FUNCTION_FittingFunction(bin_number,local_intObservable,&
                                  Cos(total_array%energy_values(local_intObservable,bin_number)%energies_datapoints(J,1)),J)
		!print *, "J,local_intObservable, s11rE(1), Yplus: ", J,local_intObservable,s11rE(1),  Y_Fit_Plus(J)
                END DO
	

	        CALL FreeParameter_Amplitude(bin_number,I, FreeParameter_original(I))
	             FreeParameter_minus(I)=FreeParameter(I)*(1.-epsilon/2.)
	        CALL FreeParameter_Amplitude(bin_number,I, FreeParameter_minus(I))

		!Now determine fit value based on a slight change of parameters the other way.
                DO J=1, NumberOfData_Obs
                    Y_FIT_minus(J)=FUNCTION_FittingFunction(bin_number,local_intObservable,&
                                  Cos(total_array%energy_values(local_intObservable,bin_number)%energies_datapoints(J,1)),J)
		!print *, "J,local_intObservable,  Yminus: ", J,local_intObservable,s11rE(1), Y_Fit_Minus(J)
                END DO

		CALL FreeParameter_Amplitude(bin_number,I, FreeParameter_original(I))

                DO J=1, NumberOfData_Obs
		    IF (FreeParameter(I).NE.0.) THEN
			dF(J+data_tally, I)=(Y_FIT_plus(J)-Y_FIT_minus(J))/(FreeParameter(I)*epsilon)
			!if(local_intobservable .eq. 5) then
		!if(bin_number .eq. 19 .and.  j.eq. 1) print *,Freeparameter(i),i,df(j+data_tally,I),data_tally,j
			!end if
                    END IF
                END DO

               !if(bin_num.eq.6)write(unit_debug,*),y_fit_plus,"ZZZZZZ",y_fit_minus,"Obs is:", local_intObservable,"XXXX"
          END IF
          END DO

		!tallys number of data cycled through from previous observables.
		data_tally = data_tally + NumberOfData_Obs
		!print *,"data_tally in fittderiv is:", data_tally,local_intobservable
	End if
	end do

	!print *,"TEST. It leaves fittingfunctionderivative."
        RETURN
        END


!gets called by FittingFunctionderivative 5 times per free parameter
!Assigns the new Freeparameter value to Pool_1 and then PW. On final call for a given parameter, it assigns the 
!Free_Parameter_Original value back to the PW. This means all changes to the PW are temporary in this section and only for computing dF()
	SUBROUTINE FreeParameter_Amplitude(bin_number,CurrentFreeParameter, ValueOfFreeParameter)
!*-----------------------------------------------------------------------
!C       FreeParameter => Pool => s01r(bin_number)..., read single amplitude to parameter set
!*-----------------------------------------------------------------------
        INTEGER CurrentFreeParameter
        DOUBLE PRECISION ValueOfFreeParameter
        INTEGER Counter_Parameter
        INTEGER :: bin_number, I1
	
	!sets counter_parameter to the current number of free parameters. Decreased later on.
        Counter_Parameter=NumberOfParameter

	!print *, "Called FPA line 2878:NumParam=", numberofparameter
        DO I1=NumberofFlag,1,-1 !see Amplitude.h, for Isospin=1

            SELECT CASE(Flag_1(I1))
            CASE (0)
                CONTINUE
            CASE (1)
		!print *, "counter, numberofparameter, and I1 is line 2311: " ,CurrentfreeParameter, counter_parameter, I1
!	print *, '2047: C_P',Counter_Parameter
!c     +                    NumberOfFreeParameter
	        IF (Counter_Parameter.EQ.CurrentFreeParameter) THEN
	           FreeParameter(Counter_Parameter)=ValueOfFreeParameter
!print *, "line 2123",freeparameter(counter_parameter), pool_1(I1)
                END IF

		!this line is needed because I1 does not equal Counter_Parameter.
                Pool_1(I1)=FreeParameter(Counter_Parameter)
                Counter_Parameter=Counter_Parameter-1
            CASE (2)!This case is for flag value of 2 or mod with fixed phase, since the mod part flag value is 0, this works ok cycling over all parameters
	        IF (Counter_Parameter.EQ.CurrentFreeParameter) THEN
	           FreeParameter(Counter_Parameter)=ValueOfFreeParameter
                END IF

		!this line is needed because I1 does not equal Counter_Parameter.
                !For Flag=2 value, the PW for a and b must be same. F15aE=F15bE etc.
                Pool_1(I1)=FreeParameter(Counter_Parameter)
                Pool_1(I1+10)=FreeParameter(Counter_Parameter)
                Counter_Parameter=Counter_Parameter-1
		CASE DEFAULT
		print *, "flag value not correct,erroring out at line 2331"
            END SELECT
        END DO
	!Pool_1->PWA value(S11,P11 etc.)
        Call Amplitude_assign(2,bin_number)

        RETURN
        END


!This just handles Directional handling of variables. A=B type work where there is a long list of variables.
!Ie pool values have been changed, then assigned back into the s11/p11 amplitudes for calculations.
!Or once calculations have run, puts the original amp values back into s11 for final results.
!Direction values mean assigning follows below: 
!1 is PWA->Pool_1
!2 is PWA->Pool_1
!3 is Pool_Original->PWA
!4 is Free_Parameter ->Pool_1
!5 is multifunctioned... temp_a->free_parameter; free and free_error-> pool_1 and pool_1_error

!the hope was this slims down Chi_Fit section so it is easier to follow.
	Subroutine Amplitude_assign(direction,bin_number,optional_A)
		integer :: direction, bin_number
		integer :: local_I1,Counter_Parameter
		double precision, optional :: optional_A(NumberOfParameter_Maximum)
		double precision :: temp_A(NumberOfParameter_Maximum)
		
		if(present(optional_A)) temp_A=optional_A

		select case(direction)
			case(1)
				!print *, "case 1 amplitude assign",pool_1(1)
				!Electric Multipoles real and imaginary
        			Pool_1(1) =s01a(bin_number)
        			Pool_1(2) =p01a(bin_number)
			        Pool_1(3) =p03a(bin_number)
        			Pool_1(4) =d03a(bin_number)
        			Pool_1(5) =d05a(bin_number)
        			Pool_1(6) =f05a(bin_number)
        			Pool_1(7) =f07a(bin_number)
        			Pool_1(8) =g07a(bin_number)
        			Pool_1(9) =g09a(bin_number)
        			Pool_1(10)=h09a(bin_number)

        			Pool_1(11)=s01b(bin_number)
        			Pool_1(12)=p01b(bin_number)
        			Pool_1(13)=p03b(bin_number)
        			Pool_1(14)=d03b(bin_number)
        			Pool_1(15)=d05b(bin_number)
        			Pool_1(16)=f05b(bin_number)
        			Pool_1(17)=f07b(bin_number)
        			Pool_1(18)=g07b(bin_number)
        			Pool_1(19)=g09b(bin_number)
        			Pool_1(20)=h09b(bin_number)

				!Magnetic Pultipoles real and imaginary
				Pool_1(21) =s11a(bin_number)
        			Pool_1(22) =p11a(bin_number)
        			Pool_1(23) =p13a(bin_number)
        			Pool_1(24) =d13a(bin_number)
        			Pool_1(25) =d15a(bin_number)
        			Pool_1(26) =f15a(bin_number)
        			Pool_1(27) =f17a(bin_number)
        			Pool_1(28) =g17a(bin_number)
        			Pool_1(29) =g19a(bin_number)
        			Pool_1(30) =h19a(bin_number)

        			Pool_1(31)=s11b(bin_number)
        			Pool_1(32)=p11b(bin_number)
        			Pool_1(33)=p13b(bin_number)
        			Pool_1(34)=d13b(bin_number)
        			Pool_1(35)=d15b(bin_number)
        			Pool_1(36)=f15b(bin_number)
        			Pool_1(37)=f17b(bin_number)
        			Pool_1(38)=g17b(bin_number)
        			Pool_1(39)=g19b(bin_number)
        			Pool_1(40)=h19b(bin_number)
			case(2)
				!print *, "case 2 amplitude assign"
        			s01a(bin_number)=Pool_1(1)
        			p01a(bin_number)=Pool_1(2)
        			p03a(bin_number)=Pool_1(3)
        			d03a(bin_number)=Pool_1(4)
        			d05a(bin_number)=Pool_1(5)
        			f05a(bin_number)=Pool_1(6)
        			f07a(bin_number)=Pool_1(7)
        			g07a(bin_number)=Pool_1(8)
        			g09a(bin_number)=Pool_1(9)
        			h09a(bin_number)=Pool_1(10)

        			s01b(bin_number)=Pool_1(11)
        			p01b(bin_number)=Pool_1(12)
        			p03b(bin_number)=Pool_1(13)
        			d03b(bin_number)=Pool_1(14)
        			d05b(bin_number)=Pool_1(15)
        			f05b(bin_number)=Pool_1(16)
        			f07b(bin_number)=Pool_1(17)
        			g07b(bin_number)=Pool_1(18)
        			g09b(bin_number)=Pool_1(19)
        			h09b(bin_number)=Pool_1(20)

				s11a(bin_number)=Pool_1(21)
        			p11a(bin_number)=Pool_1(22)
        			p13a(bin_number)=Pool_1(23)
        			d13a(bin_number)=Pool_1(24)
        			d15a(bin_number)=Pool_1(25)
        			f15a(bin_number)=Pool_1(26)
        			f17a(bin_number)=Pool_1(27)
        			g17a(bin_number)=Pool_1(28)
        			g19a(bin_number)=Pool_1(29)
        			h19a(bin_number)=Pool_1(30)

        			s11b(bin_number)=Pool_1(31)
        			p11b(bin_number)=Pool_1(32)
        			p13b(bin_number)=Pool_1(33)
        			d13b(bin_number)=Pool_1(34)
        			d15b(bin_number)=Pool_1(35)
        			f15b(bin_number)=Pool_1(36)
        			f17b(bin_number)=Pool_1(37)
        			g17b(bin_number)=Pool_1(38)
        			g19b(bin_number)=Pool_1(39)
        			h19b(bin_number)=Pool_1(40)
			case(3)
        			s01a(bin_number)=Pool_1_Original(1)
        			p01a(bin_number)=Pool_1_Original(2)
        			p03a(bin_number)=Pool_1_Original(3)
        			d03a(bin_number)=Pool_1_Original(4)
        			d05a(bin_number)=Pool_1_Original(5)
        			f05a(bin_number)=Pool_1_Original(6)
        			f07a(bin_number)=Pool_1_Original(7)
        			g07a(bin_number)=Pool_1_Original(8)
        			g09a(bin_number)=Pool_1_Original(9)
        			h09a(bin_number)=Pool_1_Original(10)

        			s01b(bin_number)=Pool_1_Original(11)
        			p01b(bin_number)=Pool_1_Original(12)
        			p03b(bin_number)=Pool_1_Original(13)
        			d03b(bin_number)=Pool_1_Original(14)
        			d05b(bin_number)=Pool_1_Original(15)
        			f05b(bin_number)=Pool_1_Original(16)
        			f07b(bin_number)=Pool_1_Original(17)
        			g07b(bin_number)=Pool_1_Original(18)
        			g09b(bin_number)=Pool_1_Original(19)
        			h09b(bin_number)=Pool_1_Original(20)

        			s11a(bin_number)=Pool_1_Original(21)
        			p11a(bin_number)=Pool_1_Original(22)
        			p13a(bin_number)=Pool_1_Original(23)
        			d13a(bin_number)=Pool_1_Original(24)
        			d15a(bin_number)=Pool_1_Original(25)
        			f15a(bin_number)=Pool_1_Original(26)
        			f17a(bin_number)=Pool_1_Original(27)
        			g17a(bin_number)=Pool_1_Original(28)
        			g19a(bin_number)=Pool_1_Original(29)
        			h19a(bin_number)=Pool_1_Original(30)

        			s11b(bin_number)=Pool_1_Original(31)
        			p11b(bin_number)=Pool_1_Original(32)
        			p13b(bin_number)=Pool_1_Original(33)
        			d13b(bin_number)=Pool_1_Original(34)
        			d15b(bin_number)=Pool_1_Original(35)
        			f15b(bin_number)=Pool_1_Original(36)
        			f17b(bin_number)=Pool_1_Original(37)
        			g17b(bin_number)=Pool_1_Original(38)
        			g19b(bin_number)=Pool_1_Original(39)
        			h19b(bin_number)=Pool_1_Original(40)
			
			Case (4)

				Counter_Parameter=NumberOfParameter
         			DO local_I1=NumberOfFlag, 1, -1 !see Amplitude.h, for Isospin=1
            			  SELECT CASE(Flag_1(local_I1))
				  !flag off
            			  CASE (0)
                		  CONTINUE
				  !flag on
            			  CASE (1)
                		  Pool_1(local_I1)=FreeParameter(Counter_Parameter)
                		  Counter_Parameter=Counter_Parameter-1
            			  CASE (2)
                		  Pool_1(local_I1)=FreeParameter(Counter_Parameter)
                		  Pool_1(local_I1+10)=FreeParameter(Counter_Parameter)
                		  Counter_Parameter=Counter_Parameter-1
            			  END SELECT
         			END DO

			Case (5)
				DO local_I1=1, NumberOfParameter
            			  FreeParameter(local_I1)=temp_A(local_I1)
        			END DO

        			Counter_Parameter=NumberOfParameter  !necesssary because of GOTO 1 loops later on
        			DO local_I1=NumberOfFlag, 1, -1 
            			  SELECT CASE(Flag_1(local_I1))
				  !flag off
            			  CASE (0)
                		    CONTINUE
				  !flag on
            			  CASE (1)
                		    Pool_1(local_I1)=FreeParameter(Counter_Parameter)
                		    Pool_1_Error(local_I1)=Error_FreeParameter(Counter_Parameter)
                		    Counter_Parameter=Counter_Parameter-1
            			  CASE (2)
                		    Pool_1(local_I1)=FreeParameter(Counter_Parameter)
                		    Pool_1(local_I1+10)=FreeParameter(Counter_Parameter)
                		    Pool_1_Error(local_I1)=Error_FreeParameter(Counter_Parameter)
                		    Pool_1_Error(local_I1+10)=Error_FreeParameter(Counter_Parameter)
                		    Counter_Parameter=Counter_Parameter-1
            			  END SELECT
        			END DO
			case(6)
				!print *, "case 6 amplitude assign"
        			s01m(bin_number)=Pool_1_Modphase(1)
        			p01m(bin_number)=Pool_1_Modphase(2)
        			p03m(bin_number)=Pool_1_Modphase(3)
        			d03m(bin_number)=Pool_1_Modphase(4)
        			d05m(bin_number)=Pool_1_Modphase(5)
        			f05m(bin_number)=Pool_1_Modphase(6)
        			f07m(bin_number)=Pool_1_Modphase(7)
        			g07m(bin_number)=Pool_1_Modphase(8)
        			g09m(bin_number)=Pool_1_Modphase(9)
        			h09m(bin_number)=Pool_1_Modphase(10)

				s11m(bin_number)=Pool_1_Modphase(21)
        			p11m(bin_number)=Pool_1_Modphase(22)
        			p13m(bin_number)=Pool_1_Modphase(23)
        			d13m(bin_number)=Pool_1_Modphase(24)
        			d15m(bin_number)=Pool_1_Modphase(25)
        			f15m(bin_number)=Pool_1_Modphase(26)
        			f17m(bin_number)=Pool_1_Modphase(27)
        			g17m(bin_number)=Pool_1_Modphase(28)
        			g19m(bin_number)=Pool_1_Modphase(29)
        			h19m(bin_number)=Pool_1_Modphase(30)
			case(7)
				!print *, "case 6 amplitude assign"
        			s01p(bin_number)=Pool_1_Modphase(11)
        			p01p(bin_number)=Pool_1_Modphase(12)
        			p03p(bin_number)=Pool_1_Modphase(13)
        			d03p(bin_number)=Pool_1_Modphase(14)
        			d05p(bin_number)=Pool_1_Modphase(15)
        			f05p(bin_number)=Pool_1_Modphase(16)
        			f07p(bin_number)=Pool_1_Modphase(17)
        			g07p(bin_number)=Pool_1_Modphase(18)
        			g09p(bin_number)=Pool_1_Modphase(19)
        			h09p(bin_number)=Pool_1_Modphase(20)

				s11p(bin_number)=Pool_1_Modphase(31)
        			p11p(bin_number)=Pool_1_Modphase(32)
        			p13p(bin_number)=Pool_1_Modphase(33)
        			d13p(bin_number)=Pool_1_Modphase(34)
        			d15p(bin_number)=Pool_1_Modphase(35)
        			f15p(bin_number)=Pool_1_Modphase(36)
        			f17p(bin_number)=Pool_1_Modphase(37)
        			g17p(bin_number)=Pool_1_Modphase(38)
        			g19p(bin_number)=Pool_1_Modphase(39)
        			h19p(bin_number)=Pool_1_Modphase(40)

			Case Default
				print*, "Direction is not getting assigned properly in Assign_Amplitudes"
		end select
	return
	End
end module Chi_Square_mod





!Module Devoted to Output------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------
!WARNING READ THIS COMMENT FIRST 
!before making changes:
!Units of the PWA are mfm, as the observable data has been changed to mfm.
!So it is output's responsibillity to make the necessary unit changes to the PWAs and observables before outputting them to a file. 
!Conversion to (mub) is as follows:
!100 (mfm)^2 = 1mub    OR    10 mfm = 1 Sqrt[mub]
!------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------
module Output_mod

contains

	Subroutine OUTPUT(out_bin_number)
	use computer_specific
	use Helicity_Amplitudes
	use types
	!Character *7, parameter :: subD = ( / 'GPtoEP', 'GPtoKL', 'Pi0toP', 'Pi+toN' / )
	!*-----------------------------------------------------------------------
	!       output results for database-each run is for only 1 bin
	!*-----------------------------------------------------------------------

	!passed in variables
	integer :: out_bin_number, counter,int_func_to_compare
	character*30::testing123
	!local variables
	Integer :: i,test_counter
	!Directory variables for output subroutines

	!Character*10 :: File_a_ending =   "Output.dat"
	!Character*46 :: file_num_bins =    directory//subdirectory//"Input/num_bins"
!	INTEGER :: Time_Initial, Time_Final, Time_Ellipse	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!below is list of observable unit numbers. The counter value should be equal to the last 2 digits of the observable unit number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!integer	:: SGTi=1003,DSGi=1004,Ti=1005,Si=1006,POLi=1007 !Single Polarization observables
	!integer :: Gi=1008,Fi=1009,Hi=1010, Ei=1011 !Double Polarization observables 
	!integer :: OXi=1012,OZi=1013,CXi=1014,CZi=1015 !Linear and Circular Pol observables
	!integer :: Txi=1016,Tzi=1017,Lxi=1018,Lzi=1019 !Transversity Observables



	!call testing_output(out_bin_number)


	counter = 4
	call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/DSG",counter)
	counter = 5
	call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/P",counter)
	counter = 6
	call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/PDSG",counter)
	!counter = 7
	!call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/F",counter)
	!counter = 8
	!call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/H",counter)
	!counter = 9
	!call OUTPUT_obs(out_bin_number,directory//subdirectory//"Output/observable_files/E",counter)

	!call modified_obs(out_bin_number)
	call output_Amplitudes(out_bin_number)
	!call Output_Amplitude_Mod(out_bin_number) !output file used to graph modulus and error
	!CALL Output_Chi_per_point(out_bin_number)

	END



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Write to FIna_Chi_Results various mod phase and propagating solution items.
!At the end it will allow ease of moving solutions from 1 bin to the next.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine output_FinalChi(out_bin_number)

	use computer_specific
	use constants
	use types
	use Helicity_Amplitudes
	use Chi_Square_mod, only:GLOBAL_ChiSum
	integer :: out_bin_number, i
	integer :: sum_data(3:num_observables+3)=0
	double precision :: TOTAL_chisum =0.

	do i=4,3+num_observables
	do j=1,out_bin_number
	sum_data(i) = sum_data(i)+total_array%energy_num_points(i,j)
	end do
	end do

	open(unit=993,file=directory//subdirectory//"Output/Final_chi_results.dat",status='unknown',position="append")
	open(unit=994,file=directory//subdirectory//"Input/NewInput/PropagateSoln.dat",status='unknown')
	write(993,'(90("="),/,/)')
	write(993,'(90("="),/,/)')
	write(993,*) "GLOBAL CHI- OBS       Number Data for observable:"
	do i=4,3+num_observables
	if(i.eq. 22)cycle
	WRITE(993,*) observable_files(i),":", GLOBAL_ChiSum(i), "||         Num Data:",sum_data(i) ,&
			 "||            Global Chi Per Point:" ,GLOBAL_ChiSum(i)/(sum_data(i)+1)
	total_chisum = total_chisum + GLOBAL_ChiSum(i)
	end do
	write(993,*) "TOTAL CHI SQUARE VALUES OVER ALL OBERSVABLES AND BINS:"
	WRITE(993,*) "TOTAL SUM OF ALL CHI VALUES:", total_chisum
	write(993,'(90("="),/,/)')
	write(993,*) "** Global Chi per point is divided by the (# datapoints + 1) to avoid divide by 0. error."	
	write(993,'(90("="),/,/)')
	!print *, "outbin number signifies" , out_bin_number
	!write(993,'(/,/)')
	!!write(993,*)  "D13Em D13Ep D13Mm D13Mp F15Em F15Ep F15Mm F15Mp"
	!do k=1,out_bin_number
        !write(993,20)   d13mE(k), d13pE(k), d13mM(k), d13pM(k), f15mE(k), f15pE(k), f15mM(k), f15pM(k)
	!end do
	!Puts a list of modulus values.
	write(994,'(/,/)')
	write(994,*), "Modulus values of the different partial waves."
	write(994,*)  "W(MeV)  S11  P11  P13  P13  D13  D13  D15  D15"
	do k=1,out_bin_number
        write(994,25)  W(k), S01m(k), P01m(k), P03m(k), P03m(k),D03m(k), D03m(k), D05m(k), D05m(k)
	end do
	write(994,*)  "W(MeV)  F15  F15  F17  F17  G17  G17  G19  G19"
	do k=1,out_bin_number
        write(994,25)  W(k), F05m(k), F05m(k),F07m(k), F07m(k), G07m(k), G07m(k), G09m(k), G09m(k)
	end do

	write(994,'(/,/)')
	write(994,*) "Use below for Propagating solutions.-Copy to intial_pwas_small"
	if(out_bin_number .eq. end_bin-1) then
	!print *,"End bin is:", end_bin, "and out_bin is:", out_bin_number
	write(994,'(a3)'),'S11'
	do k=1,end_bin-1
        write(994,30)   W(k), s01a(k), s01b(k)
	end do
	write(994,'(a3)'),'P11'
	do k=1,end_bin-1
        write(994,30)   W(k), p01a(k), p01b(k)
	end do
	write(994,'(a3)'),'P13'
	do k=1,end_bin-1
        write(994,30) 	W(k), p03a(k), p03b(k), p13a(k), p13b(k)
	end do
	write(994,'(a3)'),'D13'
	do k=1,end_bin-1
        write(994,30) 	W(k), d03a(k), d03b(k), d13a(k), d13b(k)
	end do
	write(994,'(a3)'),'D15'
	do k=1,end_bin-1
        write(994,30) 	W(k), d05a(k), d05b(k), d15a(k), d15b(k)
	end do
	write(994,'(a3)'),'F15'
	do k=1,end_bin-1
        write(994,30) 	W(k), f05a(k), f05b(k), f15a(k), f15b(k)
	end do
	write(994,'(a3)'),'F17'
	do k=1,end_bin-1
        write(994,30) 	W(k), f07a(k), f07b(k), f17a(k), f17b(k)
	end do
	write(994,'(a3)'),'G17'
	do k=1,end_bin-1
        write(994,30) 	W(k), g07a(k), g07b(k), g17a(k), g17b(k)
	end do
	write(994,'(a3)'),'G19'
	do k=1,end_bin-1
        write(994,30) 	W(k), g09a(k), g09b(k), g19a(k), g19b(k)
	end do
	write(994,'(a3)'),'H19'
	do k=1,end_bin-1
        write(994,30) 	W(k), h09a(k), h09b(k), h19a(k), h19b(k)
	end do
	write(994,'(a4)'),'H111'
	do k=1,end_bin-1
        write(994,30) 	W(k), h111a(k), h111b(k), h111a(k), h111b(k)
	end do
	end if

	write(993,*) 
	write(994,*) 
	close(unit=993)
	close(unit=994)
25	format(f6.1,2x,8(f5.2,1x))
30	format(f7.1,3x,f9.5,1x,f9.5,1x,f9.5,1x,f9.5)
	end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!Amplitudes are originally in units of mfm. It is the outputs responsibility to change units of PWAs IF DESIRED
!Converts DSG observable data AND function return value from mfm to mub. 
!PWA's are still left as mfm units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE OUTPUT_obs(out_bin_number,file_pathname,int_observable_for_fit)
	use computer_specific
	use types
	use Helicity_Amplitudes

	!passed variables
	Character*2  :: char_bin_number
	Integer :: out_bin_number,int_observable_for_fit
	character :: file_pathname*(*)
	logical :: file_exist,file_exist2

	!local variables
	integer :: i
	double precision :: fac !conversion of DSG data back into mub.
	!Character*46 :: file_num_bins = directory//subdirectory//"Input/num_bins"


	!converts the bin number to character
	if(out_bin_number .lt. 9.5) write(char_bin_number,'(i1)'), out_bin_number
	if(out_bin_number .gt. 9.5) write(char_bin_number,'(i2)'), out_bin_number

	!files to open with a new name for each bin.
		
	  OPEN(UNIT=3002, FILE=file_pathname//trim(char_bin_number)//".dat", STATUS='replace', FORM='FORMATTED')
	if(testing_variable .neqv. .true.) then
	  OPEN(UNIT=3003, FILE=file_pathname//trim(char_bin_number)//"_raw.dat", STATUS='replace', FORM='FORMATTED')
	end if
	  OPEN(UNIT=UNIT_num_bins, FILE=directory//subdirectory//"Input/num_bins.dat", STATUS='unknown', FORM='FORMATTED')
	  rewind(Unit_num_bins)
	  write(unit_num_bins,'(i2)') out_bin_number
	!files to open on first bin and then reopen as append on all future bins.
	!if(out_bin_number.eq.1) then
	 ! OPEN(UNIT=Unit_OUTPUT, FILE=FILE_a//File_a_ending, STATUS='unknown', FORM='FORMATTED')
	 ! OPEN(UNIT=3001, FILE=FILE_a//"s11rE"//".dat", STATUS='unknown', FORM='FORMATTED')
	 ! WRITE(Unit_OUTPUT, 101) fdate()
!101       FORMAT('Date is ', A24, /)
	!else
	!  OPEN(UNIT=Unit_OUTPUT, FILE=FILE_a//File_a_ending, STATUS='unknown', position='append')
	!  OPEN(UNIT=3001, FILE=FILE_a//"s11rE"//".dat", STATUS='unknown', position='append')
	!end if
		

	!file 3000 which is to /PartialWaves/output/Output.dat    
	!Write(unit_output,'(2(F9.3))',advance='yes') , total_array%central_energy(out_bin_number)
	!Write(unit_output,'(a10,3x, i3)',advance='yes') , "bin number", out_bin_number
	!Write(Unit_Output,'(2(F8.3))'), Pool_1(1), Pool_1(11)
	!Write(Unit_Output,'(2(F8.3),/)'), Pool_1(22), Pool_1(32)
	!Write(Unit_Output,'(2(F8.3),/)'), Pool_1(24), Pool_1(34)

	!file 3001 which is to PartialWaves/output/s11rE.dat NEEDS TO CONTAIN INFO FOR ALL BINS
	!write(3001, '(2(f10.5,2x))'), s11rE(out_bin_number), W(out_bin_number)

	!print *,"various pw values are:,",W(out_bin_number), s11aE(out_bin_number), s11bE(out_bin_number)
	!print *,"other various pw values are:,",&
		!p13aE(out_bin_number), p13bE(out_bin_number), p13aM(out_bin_number), p13bM(out_bin_number)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!File 3002 for 4 columns at different angles and fixed energy ONLY CONTAINS INFO FOR 1 BIN
	fac = 1. !no conversion unless case 4 is selected.
	Select Case (int_observable_for_fit)
	!@ symbol in graphing program is a space. DO NOT include a space in the various titles, instead use a @
	!Ex. "cos (x)" to get the space between the 's' and '(' would be input as "cos@(x)"
	Case (4)
	  write(3002,100), "cos@#Theta","d#sigma/d#Omega@"//units_mb,-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_Diff_Cross_Section(), 0.01 , X, 0.01
	  end do

	Case (5)
	  write(3002,100), "cos@#Theta","P@@@@@@@@@@@@@@@",-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_Polarization(), 0.01 , X, 0.01
	  end do
	Case (6)
	  write(3002,100), "cos@#Theta","Pd#sigma/d#Omega",-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_PDSG(), 0.01 , X, 0.01
	  end do
	Case (7)
	  write(3002,100), "cos@#Theta","F@@@@@@@@@@@@@@@",-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_F(), 0.01 , X, 0.01
	  end do
	Case (8)
	  write(3002,100), "cos@#Theta","H@@@@@@@@@@@@@@@",-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_H(), 0.01 , X, 0.01
	  end do
	Case (9)
	  write(3002,100), "cos@#Theta","E@@@@@@@@@@@@@@@",-1.,1.,&
			   "W@=@", int(W(out_bin_number)),user_reaction
	  do i = 0,180,4
	  X=cos(i*Pi/180.)
	  write(3002,110), func_E(), 0.01 , X, 0.01
	  end do

	Case Default
		print *, "Reached default case in output_obs subroutine. observable has not been programmed in yet."
		print *, "Please add it into the module Output. Int number:",int_observable_for_fit
	End Select

100	Format(a11,2x,a27,2x,F3.0,2x,F3.0,2x,a4,i4,2x,I3)
110	Format(4(F9.4,2x))

	!print *, "fac for conversion is:", fac
	!File 3003 which will hold original raw data for the bin. ONLY CONTAINS INFO FOR 1 BIN
	if(TESTING_VARIABLE .neqv. .TRUE.) then
	!print *,"SHOULD NOT ENTER HERE ON TESTING"

	do i=1,total_Array%energy_num_points(int_observable_for_fit,out_bin_number)
	  write(3003,'(3(F9.4,2X),i2,2x,a15,2x,a60)'),&
				 	total_array%energy_values(int_observable_for_fit,out_bin_number)%energies_datapoints(i,2)/fac,&
					total_array%energy_values(int_observable_for_fit,out_bin_number)%energies_datapoints(i,3)/fac,&
				    cos(total_array%energy_values(int_observable_for_fit,out_bin_number)%energies_datapoints(i,1)),& 
					total_array%energy_values(int_observable_for_fit,out_bin_number)%author_number(i),& 
					total_array%energy_values(int_observable_for_fit,out_bin_number)%data_points_author_name(i),& 
					total_array%energy_values(int_observable_for_fit,out_bin_number)%data_points_paper_date(i) 					
	end do
	!in the case of 0 points, write 1 line of 0's and NONE
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILE GROUP THAT WILL BE GRAPHED TOGETHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !PRINT *, 'OUTPUT finished.'
	
        !WRITE(*, 201) fdate()
!201     FORMAT(60('-'), /, 'Date is ', A24)
	!write(6,'(//)')
	!CLOSE(UNIT=Unit_OUTPUT) 

	!if(out_bin_number.eq.1) then
	!	OPEN(UNIT=9001, FILE=pion_file, STATUS='unknown', FORM='FORMATTED')
	!	write(9001,*) total_array%central_energy(out_bin_number), s11rE(out_bin_number), s11iE(out_bin_number)
	!else
	!OPEN(UNIT=9001, FILE=pion_file, STATUS='unknown', position='append')
	!write(9001,*) total_array%central_energy(out_bin_number), s11rE(out_bin_number), s11iE(out_bin_number)
	!end if
	!rewind(9001)
	!Close(Unit=9001)
	!Close(3001)
	Close(unit=3002)
	Close(unit=3003)
	Close(UNIT_num_bins)
        RETURN

        END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generates output files for the individual amplitudes.
!Unit of Amplitudes is in mfm. 
!this is used for plotting option 2 selection1-8 of the cplus program main######.cxx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine output_Amplitudes(out_bin_number)
	use computer_specific
	use constants
	use types
	use Helicity_Amplitudes
	use Chi_Square_mod

	integer :: local_i, out_bin_number
	!character *48 :: amplitudes_output_file = directory//subdirectory//"Output/amplitudes"

	if(out_bin_number .eq. 1) then
	!real amplitude files separate
	open(unit=unit_s001r,file= directory//subdirectory//"Output/amplitudes/S001r.dat",status='unknown')
	open(unit=unit_s101r,file= directory//subdirectory//"Output/amplitudes/S101r.dat",status='unknown')
	open(unit=unit_p001r,file= directory//subdirectory//"Output/amplitudes/P001r.dat",status='unknown')
	open(unit=unit_p101r,file= directory//subdirectory//"Output/amplitudes/P101r.dat",status='unknown')
	open(unit=unit_p003r,file= directory//subdirectory//"Output/amplitudes/P003r.dat",status='unknown')
	open(unit=unit_p103r,file= directory//subdirectory//"Output/amplitudes/P103r.dat",status='unknown')
	open(unit=unit_d003r,file= directory//subdirectory//"Output/amplitudes/D003r.dat",status='unknown')
	open(unit=unit_d103r,file= directory//subdirectory//"Output/amplitudes/D103r.dat",status='unknown')
	open(unit=unit_d005r,file= directory//subdirectory//"Output/amplitudes/D005r.dat",status='unknown')
	open(unit=unit_d105r,file= directory//subdirectory//"Output/amplitudes/D105r.dat",status='unknown')
	open(unit=unit_f005r,file= directory//subdirectory//"Output/amplitudes/F005r.dat",status='unknown')
	open(unit=unit_f105r,file= directory//subdirectory//"Output/amplitudes/F105r.dat",status='unknown')
	open(unit=unit_f007r,file= directory//subdirectory//"Output/amplitudes/F007r.dat",status='unknown')
	open(unit=unit_f107r,file= directory//subdirectory//"Output/amplitudes/F107r.dat",status='unknown')
	open(unit=unit_g007r,file= directory//subdirectory//"Output/amplitudes/G007r.dat",status='unknown')
	open(unit=unit_g107r,file= directory//subdirectory//"Output/amplitudes/G107r.dat",status='unknown')
	open(unit=unit_g009r,file= directory//subdirectory//"Output/amplitudes/G009r.dat",status='unknown')
	open(unit=unit_g109r,file= directory//subdirectory//"Output/amplitudes/G109r.dat",status='unknown')
	open(unit=unit_h009r,file= directory//subdirectory//"Output/amplitudes/H009r.dat",status='unknown')
	open(unit=unit_h109r,file= directory//subdirectory//"Output/amplitudes/H109r.dat",status='unknown')
	!imaginary files separate
	open(unit=unit_s001i,file= directory//subdirectory//"Output/amplitudes/S001i.dat",status='unknown')
	open(unit=unit_s101i,file= directory//subdirectory//"Output/amplitudes/S101i.dat",status='unknown')
	open(unit=unit_p001i,file= directory//subdirectory//"Output/amplitudes/P001i.dat",status='unknown')
	open(unit=unit_p101i,file= directory//subdirectory//"Output/amplitudes/P101i.dat",status='unknown')
	open(unit=unit_p003i,file= directory//subdirectory//"Output/amplitudes/P003i.dat",status='unknown')
	open(unit=unit_p103i,file= directory//subdirectory//"Output/amplitudes/P103i.dat",status='unknown')
	open(unit=unit_d003i,file= directory//subdirectory//"Output/amplitudes/D003i.dat",status='unknown')
	open(unit=unit_d103i,file= directory//subdirectory//"Output/amplitudes/D103i.dat",status='unknown')
	open(unit=unit_d005i,file= directory//subdirectory//"Output/amplitudes/D005i.dat",status='unknown')
	open(unit=unit_d105i,file= directory//subdirectory//"Output/amplitudes/D105i.dat",status='unknown')
	open(unit=unit_f005i,file= directory//subdirectory//"Output/amplitudes/F005i.dat",status='unknown')
	open(unit=unit_f105i,file= directory//subdirectory//"Output/amplitudes/F105i.dat",status='unknown')
	open(unit=unit_f007i,file= directory//subdirectory//"Output/amplitudes/F007i.dat",status='unknown')
	open(unit=unit_f107i,file= directory//subdirectory//"Output/amplitudes/F107i.dat",status='unknown')
	open(unit=unit_g007i,file= directory//subdirectory//"Output/amplitudes/G007i.dat",status='unknown')
	open(unit=unit_g107i,file= directory//subdirectory//"Output/amplitudes/G107i.dat",status='unknown')
	open(unit=unit_g009i,file= directory//subdirectory//"Output/amplitudes/G009i.dat",status='unknown')
	open(unit=unit_g109i,file= directory//subdirectory//"Output/amplitudes/G109i.dat",status='unknown')
	open(unit=unit_h009i,file= directory//subdirectory//"Output/amplitudes/H009i.dat",status='unknown')
	open(unit=unit_h109i,file= directory//subdirectory//"Output/amplitudes/H109i.dat",status='unknown')
	end if

	if(out_bin_number .ne. 1) then
	!real amplitude files separate
	open(unit=unit_s001r,file= directory//subdirectory//"Output/amplitudes/S001r.dat",status='unknown',position='append')
	open(unit=unit_s101r,file= directory//subdirectory//"Output/amplitudes/S101r.dat",status='unknown',position='append')
	open(unit=unit_p001r,file= directory//subdirectory//"Output/amplitudes/P001r.dat",status='unknown',position='append')
	open(unit=unit_p101r,file= directory//subdirectory//"Output/amplitudes/P101r.dat",status='unknown',position='append')
	open(unit=unit_p003r,file= directory//subdirectory//"Output/amplitudes/P003r.dat",status='unknown',position='append')
	open(unit=unit_p103r,file= directory//subdirectory//"Output/amplitudes/P103r.dat",status='unknown',position='append')
	open(unit=unit_d003r,file= directory//subdirectory//"Output/amplitudes/D003r.dat",status='unknown',position='append')
	open(unit=unit_d103r,file= directory//subdirectory//"Output/amplitudes/D103r.dat",status='unknown',position='append')
	open(unit=unit_d005r,file= directory//subdirectory//"Output/amplitudes/D005r.dat",status='unknown',position='append')
	open(unit=unit_d105r,file= directory//subdirectory//"Output/amplitudes/D105r.dat",status='unknown',position='append')
	open(unit=unit_f005r,file= directory//subdirectory//"Output/amplitudes/F005r.dat",status='unknown',position='append')
	open(unit=unit_f105r,file= directory//subdirectory//"Output/amplitudes/F105r.dat",status='unknown',position='append')
	open(unit=unit_f007r,file= directory//subdirectory//"Output/amplitudes/F007r.dat",status='unknown',position='append')
	open(unit=unit_f107r,file= directory//subdirectory//"Output/amplitudes/F107r.dat",status='unknown',position='append')
	open(unit=unit_g007r,file= directory//subdirectory//"Output/amplitudes/G007r.dat",status='unknown',position='append')
	open(unit=unit_g107r,file= directory//subdirectory//"Output/amplitudes/G107r.dat",status='unknown',position='append')
	open(unit=unit_g009r,file= directory//subdirectory//"Output/amplitudes/G009r.dat",status='unknown',position='append')
	open(unit=unit_g109r,file= directory//subdirectory//"Output/amplitudes/G109r.dat",status='unknown',position='append')
	open(unit=unit_h009r,file= directory//subdirectory//"Output/amplitudes/H009r.dat",status='unknown',position='append')
	open(unit=unit_h109r,file= directory//subdirectory//"Output/amplitudes/H109r.dat",status='unknown',position='append')
	!imaginary files separate
	open(unit=unit_s001i,file= directory//subdirectory//"Output/amplitudes/S001i.dat",status='unknown',position='append')
	open(unit=unit_s101i,file= directory//subdirectory//"Output/amplitudes/S101i.dat",status='unknown',position='append')
	open(unit=unit_p001i,file= directory//subdirectory//"Output/amplitudes/P001i.dat",status='unknown',position='append')
	open(unit=unit_p101i,file= directory//subdirectory//"Output/amplitudes/P101i.dat",status='unknown',position='append')
	open(unit=unit_p003i,file= directory//subdirectory//"Output/amplitudes/P003i.dat",status='unknown',position='append')
	open(unit=unit_p103i,file= directory//subdirectory//"Output/amplitudes/P103i.dat",status='unknown',position='append')
	open(unit=unit_d003i,file= directory//subdirectory//"Output/amplitudes/D003i.dat",status='unknown',position='append')
	open(unit=unit_d103i,file= directory//subdirectory//"Output/amplitudes/D103i.dat",status='unknown',position='append')
	open(unit=unit_d005i,file= directory//subdirectory//"Output/amplitudes/D005i.dat",status='unknown',position='append')
	open(unit=unit_d105i,file= directory//subdirectory//"Output/amplitudes/D105i.dat",status='unknown',position='append')
	open(unit=unit_f005i,file= directory//subdirectory//"Output/amplitudes/F005i.dat",status='unknown',position='append')
	open(unit=unit_f105i,file= directory//subdirectory//"Output/amplitudes/F105i.dat",status='unknown',position='append')
	open(unit=unit_f007i,file= directory//subdirectory//"Output/amplitudes/F007i.dat",status='unknown',position='append')
	open(unit=unit_f107i,file= directory//subdirectory//"Output/amplitudes/F107i.dat",status='unknown',position='append')
	open(unit=unit_g007i,file= directory//subdirectory//"Output/amplitudes/G007i.dat",status='unknown',position='append')
	open(unit=unit_g107i,file= directory//subdirectory//"Output/amplitudes/G107i.dat",status='unknown',position='append')
	open(unit=unit_g009i,file= directory//subdirectory//"Output/amplitudes/G009i.dat",status='unknown',position='append')
	open(unit=unit_g109i,file= directory//subdirectory//"Output/amplitudes/G109i.dat",status='unknown',position='append')
	open(unit=unit_h009i,file= directory//subdirectory//"Output/amplitudes/H009i.dat",status='unknown',position='append')
	open(unit=unit_h109i,file= directory//subdirectory//"Output/amplitudes/H109i.dat",status='unknown',position='append')


	end if

	if(out_bin_number .eq.1) then		
	!print *, "end bin w value is: " , total_array%central_energy(end_bin), " " , end_bin,total_array%central_energy(1)
	  write(unit_s001r,100), "W@(MeV)","S01@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "S01",user_reaction
	  write(unit_s101r,100), "W@(MeV)","S11@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "S11",user_reaction
	  write(unit_p001r,100), "W@(MeV)","P01@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "P01",user_reaction
	  write(unit_p101r,100), "W@(MeV)","P11@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "P11",user_reaction
	  write(unit_p003r,100), "W@(MeV)","P03@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "P03",user_reaction
	  write(unit_p103r,100), "W@(MeV)","P13@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "P13",user_reaction
	  write(unit_d003r,100), "W@(MeV)","D03@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "D03",user_reaction
	  write(unit_d103r,100), "W@(MeV)","D13@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "D13",user_reaction
	  write(unit_d005r,100), "W@(MeV)","D05@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "D05",user_reaction
	  write(unit_d105r,100), "W@(MeV)","D15@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "D15",user_reaction   
	  write(unit_f005r,100), "W@(MeV)","F05@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "F05",user_reaction
	  write(unit_f105r,100), "W@(MeV)","F15@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "F15",user_reaction   
	  write(unit_f007r,100), "W@(MeV)","F07@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "F07",user_reaction
	  write(unit_f107r,100), "W@(MeV)","F17@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "F17",user_reaction     
	  write(unit_g007r,100), "W@(MeV)","G07@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "G07",user_reaction
	  write(unit_g107r,100), "W@(MeV)","G17@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "G17",user_reaction  
	  write(unit_g009r,100), "W@(MeV)","G09@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "G09",user_reaction
	  write(unit_g109r,100), "W@(MeV)","G19@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "G19",user_reaction     
	  write(unit_h009r,100), "W@(MeV)","H09@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "H09",user_reaction
	  write(unit_h109r,100), "W@(MeV)","H19@"//units_mfm,W(1),total_array%central_energy(end_bin)+1.,&
		   "H19",user_reaction        
	end if

	do local_i=1,40
	  if(Pool_1_Error(local_i) .gt. 5) then
          print *,"Files that ED code usese errors scaled down to 2 temporarily for graphing"
          Pool_1_Error(local_i) = 1.0
          end if
	end do

	write(unit_s001r , 105) Pool_1(1),  Pool_1_Error(1) , W(out_bin_number), out_bin_number
	write(unit_s001i , 110) Pool_1(11), Pool_1_Error(11), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_s101r , 105) Pool_1(21), Pool_1_Error(21) , W(out_bin_number), out_bin_number
	write(unit_s101i , 110) Pool_1(31), Pool_1_Error(31), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_p001r , 105) Pool_1( 2), Pool_1_Error( 2), W(out_bin_number), out_bin_number
	write(unit_p001i , 110) Pool_1(12), Pool_1_Error(12), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_p101r , 105) Pool_1(22), Pool_1_Error(22), W(out_bin_number), out_bin_number
	write(unit_p101i , 110) Pool_1(32), Pool_1_Error(32), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_p003r , 105) Pool_1(3),  Pool_1_Error(3) , W(out_bin_number), out_bin_number
	write(unit_p003i , 110) Pool_1(13), Pool_1_Error(13), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_p103r , 105) Pool_1(23), Pool_1_Error(23), W(out_bin_number), out_bin_number
	write(unit_p103i , 110) Pool_1(33), Pool_1_Error(33), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_d003r , 105) Pool_1(4),  Pool_1_Error(4) , W(out_bin_number), out_bin_number
	write(unit_d003i , 110) Pool_1(14), Pool_1_Error(14), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_d103r , 105) Pool_1(24), Pool_1_Error(24), W(out_bin_number), out_bin_number
	write(unit_d103i , 110) Pool_1(34), Pool_1_Error(34), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_d005r , 105) Pool_1(5),  Pool_1_Error(5) , W(out_bin_number), out_bin_number
	write(unit_d005i , 110) Pool_1(15), Pool_1_Error(15), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_d105r , 105) Pool_1(25), Pool_1_Error(25), W(out_bin_number), out_bin_number
	write(unit_d105i , 110) Pool_1(35), Pool_1_Error(35), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_f005r , 105) Pool_1(6),  Pool_1_Error(6) , W(out_bin_number), out_bin_number
	write(unit_f005i , 110) Pool_1(16), Pool_1_Error(16), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_f105r , 105) Pool_1(26), Pool_1_Error(26), W(out_bin_number), out_bin_number
	write(unit_f105i , 110) Pool_1(36), Pool_1_Error(36), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_f007r , 105) Pool_1(7),  Pool_1_Error(7) , W(out_bin_number), out_bin_number
	write(unit_f007i , 110) Pool_1(17), Pool_1_Error(17), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_f107r , 105) Pool_1(27), Pool_1_Error(27), W(out_bin_number), out_bin_number
	write(unit_f107i , 110) Pool_1(37), Pool_1_Error(37), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_g007r , 105) Pool_1(8),  Pool_1_Error(8) , W(out_bin_number), out_bin_number
	write(unit_g007i , 110) Pool_1(18), Pool_1_Error(18), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_g107r , 105) Pool_1(28), Pool_1_Error(28), W(out_bin_number), out_bin_number
	write(unit_g107i , 110) Pool_1(38), Pool_1_Error(38), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_g009r , 105) Pool_1(9),  Pool_1_Error(9) , W(out_bin_number), out_bin_number
	write(unit_g009i , 110) Pool_1(19), Pool_1_Error(19), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_g109r , 105) Pool_1(29), Pool_1_Error(29), W(out_bin_number), out_bin_number
	write(unit_g109i , 110) Pool_1(39), Pool_1_Error(39), W(out_bin_number), out_bin_number,"temp", "temp"

	write(unit_h009r , 105) Pool_1(10), Pool_1_Error(10), W(out_bin_number), out_bin_number
	write(unit_h009i , 110) Pool_1(20), Pool_1_Error(20), W(out_bin_number), out_bin_number,"temp", "temp"
	write(unit_h109r , 105) Pool_1(30), Pool_1_Error(30), W(out_bin_number), out_bin_number
	write(unit_h109i , 110) Pool_1(40), Pool_1_Error(40), W(out_bin_number), out_bin_number,"temp", "temp"

100	Format(a7,2x,a11,2x,f5.0,2x,f5.0,2x,a4,I3)
105	Format((3(f10.5,5x),i2))
110	FORMAT(3(f10.5,5x),i2,2x,2(a4,2x))
	do local_i = 8000,8039
	Close(unit=local_i)
	end do
	close(unit=8041) !unit_all_amplitudes
	return
	end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Outputs to amplitude_modulus.dat the energy of the bin and the modulus of the different partial waves(file will include all PW's)
!Location of file is Fortran_output_folder.
!Amplitudes are determined in units mfm. Then converted to Dl for final output.
!This is for the graphing program using main######.cxx for 'root'
!would select option 9 in program main of main.cxx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Output_Amplitude_Mod(out_bin_number)
	use types
	use Helicity_Amplitudes
	use Chi_Square_mod

	!passed variables
	integer :: out_bin_number
	integer :: I
!pool_1_mp_error
	!local variables
	!character *49 :: amplitudes_output_file=directory//subdirectory//"Output/"//amplitudes_folder
	double precision :: l_fac = 0
	double precision :: dl_fac = 0

	!print *,"amplitudes_output_file is:", amplitudes_output_file
	open(unit=output_amplitude_modulus,file=directory//subdirectory//&
					"Output/amplitudes/amplitude_modulus.dat",status="unknown") 
	open(unit=7002,file=directory//subdirectory//&
					"Output/amplitudes/ED_modulus.dat",status="unknown") 
	!print *,directory//subdirectory//&
	!				"Output/amplitudes/amplitude_modulus.dat", "Output modlus file name path"
	!outputs to file the PWA amplitude modulus (not squared)
	write(output_amplitude_modulus,'(f7.1,3x,i1)') W(out_bin_number), user_reaction
	write(7002,'(f7.1,3x,i1)') W(out_bin_number), user_reaction

	 do I=1, 10 !electric amps
          !Create correct dimensionless factor
	  if((I.eq.1) .or. (I .eq. 4)) l_fac = 2
	  if((I.eq.3) .or. (I .eq. 6)) l_fac = 6
	  if((I.eq.5) .or. (I .eq. 8)) l_fac =12
	  if((I.eq.7) .or. (I .eq. 10))l_fac =20
          if (I.eq.9)                  l_fac =30	  
	  if(I .eq. 2) cycle !unphysical amplitude

	  if(pool_1_mp_error(I) .gt. 1) then
            pool_1_mp_error(I) = 1
            print *, "In output_AMplitude_Mod, pool_1_mp_error(",I,") error size decreased to 1"
	  end if

          bin_num=out_bin_number

	  write(output_amplitude_modulus,'(2(e10.4,1x))', advance='no'), &
			(Abs(Pool_1_modphase(I)))**2,2*Abs(pool_1_mp_error(I)*Pool_1_modphase(I))
	  write(7002,'(2(e10.4,1x))', advance='no'), (ED_Mod(out_bin_number,I))**2,0.00
         end do

!MAGNETIC CONVERTED TO DIMENSIONELSS MODULUS VALUES
	  write(output_amplitude_modulus,'(" ")')!Advances to new line when all Electriv PW's have been output.
	  write(7002,'(" ")')!Advances to new line when all Electriv PW's have been output.

	 do I=1, 10 !magnetic amps
	  if(pool_1_mp_error(I+20) .gt. 5) then
            pool_1_mp_error(I+20) = 1
            print *, "In output_AMplitude_Mod, pool_1_mp_error(",I+20,") error size decreased to 1"
	  end if

	  !Select bin for kinematic variables and compute conversion factor to dimensionless
          bin_num=out_bin_number

	  !!for error below, DL is squared because one comes from error conversion and the other from mod conversion ie: (A +- err)**2 gives A +- 2A*err, but A and err both contain a DL_fac
	  write(output_amplitude_modulus,'(2(e10.3,1x))', advance='no'), &
			(Abs(Pool_1_modphase(I+20)))**2,2*Abs(pool_1_mp_error(I+20)*Pool_1_modphase(I+20))
	  write(7002,'(2(e10.3,1x))', advance='no'), (ED_Mod(out_bin_number,I+20)*DL_fac)**2,0.0
	!print *,"1/(Conversion to Dimensionless)", 1/(DL_fac) 
	 end do
	  write(output_amplitude_modulus,'(" ")')!Advances line once Magnetic PWs have been written on second line.
	  write(7002,'(" ")')!Advances line once Magnetic PWs have been written on second line.

100	Format(a7,2x,a14,2x,f5.0,2x,f5.0,2x,a4,I3)
	return
	end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OUtputs to file Chi_per_point for a given bin the chi per point values.
!file is located in Fortran_Output
!this is for a general picture of how the fit is working, seeing if points are bad, and need for scaling.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine Output_Chi_per_point(out_bin_number)
	use constants
	use types
	use Helicity_Amplitudes
	use Chi_Square_mod

	!passed variables
	integer :: out_bin_number
	
	!local variables
        integer :: user_input
	character *3 :: char_bin_number
	character *40:: Date_header
	character *40:: Column_Names_header
	character *15:: AuthorInfo(NumberofStack_Maximum)
	character*7, dimension(4:3+num_observables) :: observable_names= (/"DSG    ","P      ","PDSG   ",&
				"F      "/)
	integer :: i,row_j,III !cycles through all Chi points initial then final. 
	integer :: data_tally_chi
	double precision :: Chi_values(2,NumberofStack_Maximum) !stores initial chi values in (1,X) and final in (2,X)
	double precision :: Final_Y_Fit(NumberofStack_Maximum)
	double precision :: Y_Data(NumberofStack_Maximum)
	double precision :: Y_Data_error(NumberofStack_Maximum)
	double precision :: dump_variable, dump_variable2

	
	integer :: total_number_data

 	total_number_data = NumberOfData
	!print*, "FIX chiperpoint subroutine for new unit numbers of observables.Type 99 to exit. Any other number continue"
	!print*, "If you continue chiperpoint files will be incorrect."
        !read(5,*) user_input
	!if(user_input .eq. 99 ) stop "User term to avoid CHiperPointFile"
	!integer	:: SGTi=1003,DSGi=1004,Ti=1005,Si=1006,POLi=1007 !Single Polarization observables
	!integer :: Gi=1008,Fi=1009,Hi=1010, Ei=1011 !Double Polarization observables 
	!integer :: OXi=1012,OZi=1013,CXi=1014,CZi=1015 !Linear and Circular Pol observables
	!integer :: Txi=1016,Tzi=1017,Lxi=1018,Lzi=1019 !Transversity Observables

	!generate AuthorInfo in same form as Y_Data and Y_Data_Error so I can print an author name next to chi values
	row_j=0
	do III=4,num_observables
	  do i=1,total_array%energy_num_points(III,out_bin_number)
	    AuthorInfo(i+row_j) = total_array%energy_values(III,out_bin_number)%data_points_author_name(i)
	  !if(out_bin_number .eq. 9) print *,"row_j and i are:", row_j,i
	  end do
	  row_j=row_j + total_array%energy_num_points(III,out_bin_number)
        end do

	write(char_bin_number,'(i2)') out_bin_number
	Open(ChiperPoint_file,file=Directory//subdirectory//"Output/observable_files/"//&
				"ChiPerPoint"//trim(adjustl(char_bin_number))//".dat",status='unknown')


	rewind(ChiperPoint_file)
!reads the file ChiperPoint_file to get the initial and final chi values together.
	read(ChiperPoint_file, '(a40)',END=90), Date_header
	read(ChiperPoint_file, '(a40)',END=90), Column_Names_Header
	read(ChiperPoint_file, '(a40)',End=90) 

	!write(6,*) Date_header, Column_Names_Header,"///", test

	do row_j=1,total_number_data
	  read(ChiperPoint_file, '(t3,f7.3,t10,f9.3,t23,f8.2,t37,f9.3)',End=90,iostat=iostat_bad_read) &
		dump_variable,dump_variable, dump_variable2, Chi_values(1,row_j) 
	  !read(ChiperPoint_file, '(t11,f9.3,t23,f6.3,t38,f6.1)',advance='yes',End=90,iostat=iostat_bad_read)&
	!		 dump_variable, dump_variable2, Chi_values(1,row_j)  
!	  if(row_j .lt. 10)print *,"chi1 is:",chi_values(1,row_j)
	end do
	  read(ChiperPoint_file, '(a40)',End=90)
	!print *, "Column_Names_Header variable is: " , Column_Names_Header

	do row_j=1,total_number_data
	  read(ChiperPoint_file, '(t2,f7.3,t10,f9.3,t23,f8.2,t37,f9.3)',End=90,iostat=iostat_bad_read)&
		 Final_Y_Fit(row_j), Y_Data(row_j), Y_Data_error(row_j), Chi_values(2,row_j)  
!	  if(row_j .lt. 10)print *,"chi2 is:",chi_values(2,row_j)
	end do

90	continue
	rewind(ChiperPoint_file)

	Close(ChiperPoint_File)
	Open(ChiperPoint_File,file=Directory//subdirectory//"Output/ChiFiles/"//&
				"ChiPerPoint"//trim(adjustl(char_bin_number))//".dat",status='unknown')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start writing to best form Chi file. Include Some headers for informational purposes
!including Date,bin number, fixed non-zero PWAs, and varying PWAs
	write(ChiperPoint_file, '(a40,/,a11,i2,/)'),Date_Header, "Bin Number:", out_bin_number
        write(ChiperPoint_file, 20)
        Do row_j=1,10
        if((Pool_1(row_j) .gt. 0.001) .or. (flag_1(row_j) .eq. 1))&
		write(ChiperPoint_file, 30) PartialWave_1(row_j), VaryCode(Flag_1(row_j)), VaryCode(Flag_1(row_j+10))
        End DO
        Do row_j=21,30
        if((Pool_1(row_j) .gt. 0.001) .or. (flag_1(row_j) .eq. 1))&
		write(ChiperPoint_file, 30) PartialWave_1(row_j-10), VaryCode(Flag_1(row_j)), VaryCode(Flag_1(row_j+10))
        End DO

	data_tally_chi=0
		Do int_obs=4,3+num_observables-3
                    if(total_array%energy_num_points(int_obs,out_bin_number) .ne. 0) then
	            write(ChiperPoint_file, '(a7)'), observable_names(int_obs)
	            write(ChiperPoint_file, 50) Column_Names_Header(1:7),"Y Data  ","Y Data Er", "Init Chi", "Final Chi"
		    DO row_j=1, total_array%energy_num_points(int_obs,out_bin_number) 
	 		 write(ChiperPoint_file, 45) Final_Y_Fit(row_j+data_tally_chi),&
						     Y_Data(row_j+data_tally_chi),&
						     Y_Data_error(row_j+data_tally_chi),&
		                                     Chi_values(1,row_j+data_tally_chi),&
						     Chi_values(2,row_j+data_tally_chi),&
						     AuthorInfo(row_j+data_tally_chi),&
						     total_array%energy_values(int_obs,out_bin_number)%reactionID(row_j)
        	    END DO
                    end if
		    data_tally_chi = data_tally_chi + total_array%energy_num_points(int_obs,out_bin_number)
		End Do

!Finish Writing to final Chi file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
20    Format(1x,"PWA ",3x,"R_VC",3x,"I_VC",2x,"(All other waves: VC=0 and 0.0)")
30    Format(a5,2(5x,a2))
40    Format(3x,2(f7.2,2x),3x,f6.3,4x,f8.4,t45,f8.4,2x,a15)
45    Format(3x,2(f7.2,2x),f9.3,t34,f5.2,t45,f5.2,t56,a15)
50    Format(t5,a7,t13,a8,t24,a9,t34,a8,t45,a9,2x,"Author")
!30    Format(t5,a7,t13,a8,t24,a9,t35,a8,t46,a9)
	Close(ChiperPoint_file)
	end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This output section is for physica and plotting energy dependent information.
!It outputs all the amplitudes into a single folder in the correct manner to be read by the energy dep code.
!in this subroutine, amplitudes are input in unit mfm and get converted to dimensionless
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine create_DL_amps(out_bin_number)
	use computer_specific
	use constants
	use Helicity_Amplitudes

	!passed variables
	integer :: out_bin_number

	!local variables
	integer :: j,read_data
	integer :: iostatus
	integer :: l_fac
	double precision :: Rl, Rl_Err, Im, Im_Err, Energy
	character*40 :: header
	character*50 :: gpep_energy_dep_fit_file= directory//PWA15_folder//PWA15_internal_folder
	character*10 :: Amplitude_filename_ending(40) = &
			(/"S101rE.dat", "P101rM.dat", "P103rE.dat", "P103rM.dat",&
			  "D103rE.dat", "D103rM.dat", "D105rE.dat", "D105rM.dat",& 
			  "F105rE.dat", "F105rM.dat", "F107rE.dat", "F107rM.dat",&
			  "G107rE.dat", "G107rM.dat", "G109rE.dat", "G109rM.dat",&
			  "H109rE.dat", "H109rM.dat", "H111rE.dat", "H111rM.dat",&
										 
			  "S101iE.dat", "P101iM.dat", "P103iE.dat", "P103iM.dat",&
			  "D103iE.dat", "D103iM.dat", "D105iE.dat", "D105iM.dat",& 
			  "F105iE.dat", "F105iM.dat", "F107iE.dat", "F107iM.dat",&
			  "G107iE.dat", "G107iM.dat", "G109iE.dat", "G109iM.dat",&
			  "H109iE.dat", "H109iM.dat", "H111iE.dat", "H111iM.dat"  /)
	!month for filename variables(also local)
	integer     :: date_for_file(3)
	character*3 :: month
	character*4 :: year

	call idate(date_for_file)
	!Not used to make it so I don't have to update .com files in PWA15Aug2013 every month, instead every year.
	if(date_for_file(2) .eq. 1)  month = "Jan"
	if(date_for_file(2) .eq. 2)  month = "Feb"
	if(date_for_file(2) .eq. 3)  month = "Mar"
	if(date_for_file(2) .eq. 4)  month = "Apr"
	if(date_for_file(2) .eq. 5)  month = "May"
	if(date_for_file(2) .eq. 6)  month = "Jun"
	if(date_for_file(2) .eq. 7)  month = "Jul"
	if(date_for_file(2) .eq. 8)  month = "Aug"
	if(date_for_file(2) .eq. 9)  month = "Sep"
	if(date_for_file(2) .eq. 10) month = "Oct"
	if(date_for_file(2) .eq. 11) month = "Nov"
	if(date_for_file(2) .eq. 12) month = "Dec"

	write(year,'(i4)') date_for_file(3)

	!below writes to one of the input files for the energy dep code.

	write(unit_all_amplitudes,'("ENERGY   REAL   REAL_ERR   IMAGINARY   IMAGINARY_ERR")')
	!open all amplitudes files and generate one compound file for Energy Dependent fit program.
	do j = 8000, 8019
	!get l value for conversion factor from dimensioned to dimensionless. See for instance Dissertation by Neboh for equation 
	!2-61 and 2-62. formulas are either l*(l+1) or (l+1)*(l+2) and typed into code in that format.

	  !opens real and imaginary parts of the files.
	  open(unit=j,file=directory//subdirectory//"Output/amplitudes/"&
				//Amplitude_filename_ending(j-7999),status="unknown")
	  open(unit=j+20,file=directory//subdirectory//"Output/amplitudes/"&
				//Amplitude_filename_ending(j-7979),status="unknown")

 	  read(j,*,iostat=iostatus) header

	  !if the file opens and there exists a header line then the file has data. Otherwise skip the rest of the read
	  if(iostatus .eq. 0) then
	  write(unit_all_amplitudes,'(a5)') Amplitude_filename_ending(j-7999)(1:2)//&
			Amplitude_filename_ending(j-7999)(4:4)//Amplitude_filename_ending(j-7999)(6:6)
	  write(unit_all_amplitudes,*), out_bin_number !needed for final count of bins so physica knows how many datapoints to read.
          bin_num = 0
	  do read_data=1,50
	    read(j,100,end =10) Rl,Rl_Err, Energy 
	    read(j+20,100) Im, Im_Err, Energy
	    if(abs(Rl_err) .lt. .01) Rl_err = .06
 	    if(abs(Im_err) .lt. .01) Im_err = .06	  
            bin_num = bin_num + 1 !bin_num is the Helicity amplitude module bin_number of fit. makes sure q and k are correctly computed

	    write(unit_all_amplitudes,110) Energy, Rl, Rl_Err, Im, Im_Err
	  end do
	  else
	    print *, "The Amplitude:" , Amplitude_filename_ending(j-7999)(1:4)//&
			Amplitude_filename_ending(j-7999)(6:6), "contains no data!"
	  end if
10	    continue

	close(j)
	close(j+20)
	end do

100	Format(3(f10.5,5x))
110	Format(f7.1,2x,4(f11.5,3x))


	  print *, "PINGN.F reads #bin(#datapoints) from file."
	  print *, "User does NOT NEED TO ACCOUNT FOR MORE BINS"

	end

!Anything related to files that use the Total cross section for plotting is found here. Called at end of program after all bins finished.
	subroutine total_Cross_section()
	use computer_specific
	use constants
	use types
	use Helicity_Amplitudes
	use Chi_Square_Mod

	Double Precision :: C_SGT !SGT conversion factor

	end



end module output_mod



!GUIDE TO MAIN PROGRAM AND SUBROUTINES, FUNCTIONS, AND MODULES

!MODULES-------------------------------------------------------------------------------------------------------------------------------------
!mod CONSTANTS gives definitions to commonly used variables and physical constants.

!mod AMPLITUDE defines all the amplitude variables needed for the subprogram "DEFINED LATER " to compute partial wave analysis

!mod LEGENDRE_POLY defines the first 7 Legendre Polynomials and its first and second derivatives as functions to be called later on

!mod TYPES creates a user defined storage for all data in the data file. the datapoints type stores all angles, diff cross sections and its error,
!Datapoints is a table where each row holds each of its three values in column 1,2 and 3 respectively.
!energies stores the central energy value of a bin,number of points in a bin and connects those values to the datapoints. The first array element in energies refers
!to the observable(sigma, Polarization..., the second refers to different bins. 

!PROGRAM--------------------------------------------------------------------------------------------------------------------------------------
!extract_file program accepts user input for reaction type, energy range, and bin width and then makes the call to subroutines
!which create files for each observable, creates a table, or runs through the various data points

!SUBROUTINES----------------------------------------------------------------------------------------------------------------------------------
!CASE_SELECT takes the input and generates 1 file per observable in the following form:
!first line includes energy in CM frame, number of data points at that energy for observable and observable name
!next number of lines includes the actual data collected in 3 sets per line with each set including:angle, diff_cross section for angle, and error
!Igor Strakovsky generously provided most of the data and the datafiles.

!DATA_POINTS is a subroutine that is called from case_select which actually writes to the files generated the data points. this allows
!data points to be handled individually and changed easily if the file form is changed

!MAKE_TABLE is a temp subroutine that reads the new output files and analyzes how many data points are found within
!a given bin width. This allows us to get a feel for the best bin width to use for calculations later on. It also reads the files into 
!a large two dimensional array of form: each colum is a different observable, each row is the number of points in a given bin of width 
!bin_width(variable in program). It may be changed in the future to a user_defined Type to allow a different handling of the data. 

! this program is designed to obtain data from an data input file and then perform a chi square minimization process.
!Any new reaction not included below needs to be setup in subroutine Initialize(). Also, While Isospin 3/2 reactions
!  are programmed in, currently only isospin 1/2 parameters can be varied.
	program extract_file
	use Legendre_Poly
	use constants
	use types
	use Chi_Square_mod
	use Output_mod


	! Energy Variables
	double precision :: min_energy_lab, max_energy_lab
	double precision :: E_Minimum
	double precision :: E_Maximum
	integer :: current_bin
	integer :: bin_direction
	logical :: i_opened !checks for open files at end of fitting
	character *40 :: i_name
!this flag is set to 1 if a bin has been read before. 
!This means the fit will read the previous bin read as input rather than a guess random starting value.
	integer :: flag_bin_read(nE) 
	flag_bin_read = 0

	call idate(datetoday)

! this section pertains to user selected reaction type to analyze	
         print *, "What reaction would you like to look at:"
         print *, "1)GP to Eta P"
         print *, "2)GP to K+ Lambda"
	 print *, "3)GP to Pi+ Neutron"
	 print *, "4)GP to Pi0 Proton" 
	 print *, "5)GN to Eta N"
	 print *, "6)GN to K0 Lambda"
	 print *, "7)GN to Pi- Proton" 
	 print *, "8)GN to Pi0 Neutron" 
	 print *, "9)K Long to Lambda Pi+" 
	 print *, "10)K Long to S0Pi+" 
	 print *, "11)K Long to K Short"
	 print *, "12)K- p to K- p"
	 print *, "13)K- p to Sigma- pi+"
	 print *, "14)K- p to Sigma+ pi-"
	 print *, "15)K- p to Sigma0 pi0"
	 print *, "16)K- p to kbr0 neutron"!Haora1st added
	 print *, "99)Exit Program"

         read *, user_reaction
	 if (user_reaction.eq.99) then
	        stop "Early User termination of program"
	 end if

	call Initialize(user_reaction)!Sets up the main program variables for a given reaction

	!opens original data file for reading all data for all observables
	open (unit=1000,file=directory//subdirectory//datafilename,blank='null',action="read",status="old")

	
	call case_select() 
	!only need min/max energies for this call

	min_energy_lab = (Threshold_Energy_Wcm**2-mass_nucleon**2)/(2*mass_nucleon)
	max_energy_lab = (2500**2-mass_nucleon**2)/(2*mass_nucleon)
	E_minimum = Threshold_Energy_Wcm-.5 !Set min energy at a small amount below threshold to pick up any points below threshold
	E_maximum = 2500.	!Max energy set to 2500. here and ALSO in case select. No need to analyze data above 2500MeV

	print *, min_energy_lab, max_energy_lab,   "Reaction Data Energy Range in lab frame"
	print *, E_minimum, E_maximum,  "Reaction Data Energy Range in CM frame"


	call make_table(E_Minimum,E_Maximum)

        close(unit=1000) !close raw input file of observable data


	!print *, "Number of bins: ", num_bins
!stop "EARLY TERMINATION FOR TESTING"
	bin_direction = 1
	current_bin = 1
	do while(bin_direction .ne.2)
	  if(current_bin .eq. 1) then
		!print *, "you are on the first bin, direction after fit will automatically be set to forward"
		call fit(current_bin,flag_bin_read(current_bin),user_reaction)
	  	call output(current_bin)
		flag_bin_read(current_bin) = 1
		current_bin = current_bin + 1
		if(current_bin .gt. end_bin) end_bin = current_bin
	  else if( current_bin .gt. 1) then


	  	call fit(current_bin,flag_bin_read(current_bin),user_reaction)
	  	call output(current_bin)
		flag_bin_read(current_bin) = 1
		!print *, "Enter 0 for reverse bin fitting."
		!print *, "Enter 1 for forward bin fitting."
		!print *, "Enter 2 for terminating fitting."
		!read(5,'(i1)') bin_direction
		bin_direction = 1
		if     (bin_direction .eq. 0) then 
		  current_bin = current_bin -1
		  print *, "WARNING!!! Certain output files for graphing are not setup &
			&    for backward bin_direction calculation including PWA(W)"
		!else if(bin_direction .eq. 1) then
		!  current_bin= current_bin  +1
		!  if(current_bin .gt. end_bin) end_bin = current_bin
		else if(bin_direction .eq. 1) then
		  current_bin= current_bin  +1
		  if(current_bin .gt. end_bin) then
		     !call create_DL_amps(current_bin-1)
			print *, "Final output files have been generated"
			exit
		  end if
		else
		  !call create_DL_amps(current_bin) 
		  print*, "Final output files have been generated" 
		  exit
		end if

	  else 
		!call fit(current_bin,flag_bin_read(current_bin),user_reaction)
	  	!call output(current_bin)
		flag_bin_read(current_bin) = 1
		print *, "Program did not receive a valid input, forward has been set as default."
		current_bin= current_bin +1
		if(current_bin .gt. end_bin) end_bin = current_bin
	  end if

	  !call output_debug(i)
	end do
        end_bin = current_bin
	!call total_Cross_Section() !output file to graph total cross section as function of energy bins.
	!call output_FinalChi(end_bin-1) !append to end of Final CHi Results file propagation info.
	  print *,"---------------------------------------------------------------------"
	  print *,"---------------------------------------------------------------------"
	  print *, "               Extract_File OUTPUT Terminated successfully"
	if(user_reaction .gt. 5) then
	print *, "User_reaction:", user_reaction, "has not been setup for DL graphing"
	print *, "Some features may not work/be setup for this reaction, especially with"
	print *, "Energy Dep code"
	end if
	  print *,"---------------------------------------------------------------------"
	  print *,"---------------------------------------------------------------------"
	print *,"**Run conversion.F95 to:"
	print *,"1:convert PWA's back into Dimensioned Amplitudes" 
	print *,"2:store new amplitudes values into INPUT FILE Input_PWAs.dat"
        print *,"-Program will be ready to run with new values."
	!do current_bin = 1,end_bin
	!  call output(current_bin)
	!end do
	do j=70,85
	close(unit=8000+j)
        end do
        close(unit=5000)
	close(unit=unit_debug)
        do j=10,10000

	!Make sure any files left open are closed.
	Inquire(j,opened=i_opened,name=i_name)
        if(i_opened) then
	  !print *,"Files left open are:", j,i_opened,i_name
          close(unit=j)
	end if
	end do

	return
	end program extract_file 
	


	subroutine initialize(IntReac)
	use SetupReaction
	use types
	!passed variables
	integer :: IntReac
	integer :: central_count

	double precision :: Plab_loc
	!Info on each global variable:
	!SubDirectory is the folder located on the desktop for the given reaction.
	!datafilename contains the name of the experimental data file. It must have 10 character so rename if necessary
	!H1C, and H3C are constants in front of the Helicity amplitudes. 
	!C0,C1,C2 are Isospin constants as well for Hadronic Reactions.
	!C0 is Isospin 0; C1 is primary Isospin 1 constant. C2 is an extra Isospin 1 constant.	
		!See paper R Arndt Phys Rev C Vol 42 #5, Nov 1990 page 1854 for their value for each reaction channel
		!Also, For isospin 1/2 only reaction, set H3C to 0.0 and then H1C to 1.0
	!Mass_outgoing 1 and 2 are masses of the two final state particles. If doesn't matter which is which since the difference is squared.
	!Mass_nucleon is initial target which is either proton or neutron in all current cases. (While name doesn't suggest it, program SHOULD work okay with other targets)
	!Mass_incident is rest Mass of incoming particle.
	!Threshold_Energy_Wcm is the threshold in CM frame of the reaction.

	if(IntReac .eq. 1) then
	SubDirectory = 'GPtoEP/'
	datafilename = 'gpep090714'
	mass_nucleon = 938.27 
	mass_outgoing1 = 938.27
	mass_outgoing2 = 547.51
	Threshold_Energy_Wcm = 1485.782
	H1C = 1.0
	H3C = 0.0
	end if

	!Reaction is Gamma proton -> K+ Lambda
	if(IntReac .eq. 2) then
	SubDirectory = 'GPtoKL/'
	datafilename = 'gpkl060311'
	mass_nucleon = 938.27 
	mass_outgoing1 = 1115.68
	mass_outgoing2 = 493.677
	Threshold_Energy_Wcm = 1613.33
	H1C = 1.0
	H3C = 0.0
	end if

	!Reaction is Gamma proton -> pi+ neutron
	if(IntReac .eq. 3) then
	SubDirectory = 'GP2PIN/'
	datafilename = 'pionPN.dat'
	mass_nucleon = 938.27 
	mass_outgoing1 = 939.565
	mass_outgoing2 = 139.57
	Threshold_Energy_Wcm = 1079.135
	H1C = 1.414213562
	H3C = -.4714045208
	end if

	!Reaction is Gamma proton -> pi0 proton
	if(IntReac .eq. 4) then
	SubDirectory = 'GP2PIP/'
	datafilename = 'pionPN.dat'
	mass_nucleon = 938.27 
	mass_outgoing1 = 938.27
	mass_outgoing2 = 134.977
	Threshold_Energy_Wcm = 1073.247
	H1C = 1.0
	H3C = (2.0/3.0)
	end if

	!Reaction is Gamma neutron -> eta neutron
	if(IntReac .eq. 5) then
	SubDirectory = 'GNtoEN/'
	datafilename = 'neutronDSG'
	mass_nucleon = 939.565 
	mass_outgoing1 = 939.565
	mass_outgoing2 = 547.51
	Threshold_Energy_Wcm = 1487.075
	H1C = 1.0
	H3C = 0.0
	end if

	!Reaction is Gamma neutron -> K0 Lambda
	if(IntReac .eq. 6) then
	SubDirectory = 'GNtoKL/'
	datafilename = 'gnkltest01'
	mass_nucleon = 939.565 
	mass_outgoing1 = 1115.68
	mass_outgoing2 = 497.65
	Threshold_Energy_Wcm = 1613.33
	H1C = 1.0
	H3C = 0.0
	end if

	!Reaction is Gamma neutron -> pi- proton
	if(IntReac .eq. 7) then
	SubDirectory = 'GN2PIP/'
	datafilename = 'pionPN.dat'
	mass_nucleon = 939.565 
	mass_outgoing1 = 938.27
	mass_outgoing2 = 139.57
	Threshold_Energy_Wcm = 1077.84
	H1C = 1.414213562
	H3C = .4714045208
	end if

	!Reaction is Gamma neutron -> pi0 neutron
	if(IntReac .eq. 8) then
	SubDirectory = 'GN2PIN/'
	datafilename = 'pionPN.dat'
	mass_nucleon = 939.565 
	mass_outgoing1 = 939.565
	mass_outgoing2 = 134.977
	Threshold_Energy_Wcm = 1074.542
	H1C = -1.0
	H3C = (2.0/3.0)
	end if

	!Reaction is KL p -> Pi+ Lambda HADRONIC    KLKS.dat 4 letter reaction ID: KLPL
	if(IntReac .eq. 9) then
	SubDirectory = 'KLPPiL/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 497.614
	mass_outgoing1 = 497.614
	mass_outgoing2 = 938.27
	Threshold_Energy_Wcm = 1480.00
	C0 = 0.0
	C1 = (0.70710678)
	C2=0
	end if

	!Reaction is KL p -> Sigma0 Pi+ HADRONIC    KLKS.dat 4 letter reaction ID: KLPS
	if(IntReac .eq. 10) then
	SubDirectory = 'KLPiS0/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 497.614
	mass_outgoing1 = 1192.64
	mass_outgoing2 = 139.57
	Threshold_Energy_Wcm = 1450.00
	C0 = 0.0
	C1 = (-0.70710678)
	C2= 0
	end if

	!Reaction is KL p -> KS p    KLKS.dat 4 letter reaction ID: KLKS
	if(IntReac .eq. 11) then
	SubDirectory = 'KLtoKS/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 497.614
	mass_outgoing1 = 938.27
	mass_outgoing2 = 497.614
	Threshold_Energy_Wcm = 1450.00
	C0 = 0.25
	C1 = 0.25
	C2= -0.5
	end if


	!Reaction is K- p -> K- p    KLKS.dat 4 letter reaction ID: K-pE
	!Haoranrevised here, changed the subdirectory name and the reaction notes
	if(IntReac .eq. 12) then
	SubDirectory = 'k-pk-p/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 493.677
	mass_outgoing1 = 938.27
	mass_outgoing2 = 493.677
	Threshold_Energy_Wcm = 1450.00
	C0 = 0.5
	C1 = 0.5
	C2=  0.0
	end if


	!Reaction is K- p -> Sigma- Pi+ HADRONIC    KLKS.dat 4 letter reaction ID: K-S-
	if(IntReac .eq. 13) then
	SubDirectory = 'K-p2S-/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 493.677
	mass_outgoing1 = 1197.449
	mass_outgoing2 = 139.57
	Threshold_Energy_Wcm = 1450.00
	C0 = -Sqrt(1.0/6.0)
	C1 = 0.5
	C2=  0.0
	end if


	!Reaction is K- p -> Sigma+ Pi- HADRONIC    KLKS.dat 4 letter reaction ID: K-S+
	if(IntReac .eq. 14) then
	SubDirectory = 'K-p2S+/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 493.677
	mass_outgoing1 = 1189.37
	mass_outgoing2 = 139.57
	Threshold_Energy_Wcm = 1450.00
	C0 = -Sqrt(1.0/6.0)
	C1 = -0.5
	C2=  0.0
	end if

	!Reaction is K- p -> Sigma0 Pi0 HADRONIC    KLKS.dat 4 letter reaction ID: K-S0
	if(IntReac .eq. 15) then
	SubDirectory = 'K-p2S0/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 493.677
	mass_outgoing1 = 1192.67
	mass_outgoing2 = 134.977
	Threshold_Energy_Wcm = 1450.00
	C0 = Sqrt(1.0/6.0)
	C1 = 0.0
	C2=  0.0
	end if
	!Haoran1st added
	!Reaction is K- p -> Kbr0 n HADRONIC    KLKS.dat 4 letter reaction ID: K-Kb
	if(IntReac .eq. 16) then
	SubDirectory = 'Kbr0_n/'
	datafilename = 'KLKS.dat'
	mass_nucleon = 938.27 
	mass_incident= 493.677
	mass_outgoing1 = 497.614
	mass_outgoing2 = 939.565
	Threshold_Energy_Wcm = 1431.947
	C0 = -0.5
	C1 = 0.5
	C2=  0.0
	end if
!--------------------------------------------------------------------------------------------------
!	THIS table gives ISOSpin coeff. for K-p reactions to final state listed in last column.
!			                C0
!     +			/
!     +			 0.5,			!Kaon- Proton 
!     +			-0.5,			!KaonBar0 Neutron
!     +			 0.,			!Lambda0 Pion0
!     +			-0.40824829,		!Sigma- Pion+
!     +			 0.40824829,		!Sigma0 Pion0
!     +			-0.40824829		!Sigma+ Pion-
!     +                  0.25		 	!Klong -> K short
!     +			/
	
!			                C1
!     +			/
!     +			 0.5,			!Kaon- Proton
!     +			 0.5,			!KaonBar0 Neutron
!     +			 0.70710678,		!Lambda0 Pion0
!     +			 0.5,			!Sigma- Pion+
!     +			 0.,			!Sigma0 Pion0
!     +			-0.5			!Sigma+ Pion-
!     +                  0.25                   !Klong -> K short

!					C2
!			-0.5			!Klong -> Kshort - extra isospin 1 const. Used for predictions from K-p and KN reactions.
!---------------------------------------------------------------------------------------------------


	!initialize the Central Energy and width of each bin for the reaction.
	!Warning-if this is changed, starting PW values will not be correct and will need an update using 0 iteration fit of ED code after this is run
	!Steps:Run this code knowing values are bad.
	!2:Run all PWs in ED code using 0 iteration fit.
	!3:Run conversion program to convert PW's (This has been automated so running a .com file (s11.com,p11.com etc.) will finish and run the conversion.
	!4:insert new PW values into starting PW values file.
	!*IF more bins are added, a dummy input line needs to be added to each input file for each PW.
	do central_count=1,nE

         if(user_reaction .eq. 1) then
	!Bin width of 5
	if(central_count .le. 11) then
	  total_array%central_energy(1)=1490.0
	  total_array%central_energy(2)=1500.0
	  total_array%central_energy(3)=1510.0
	  total_array%central_energy(4)=1520.0
	  total_array%central_energy(5)=1530.0
	  total_array%central_energy(6)=1540.0
	  total_array%central_energy(7)=1550.0
	  total_array%central_energy(8)=1560.0
	  total_array%central_energy(9)=1570.0
	  total_array%central_energy(10)=1580.0
	  total_array%central_energy(11)=1590.0
	  bin_width(central_count) = 10.
	else if(central_count .gt. 11 .and. (central_count.le.17)) then
	  bin_width(central_count) = 20
	  total_array%central_energy(central_count) = 1605.0+(central_count-12)*bin_width(central_count)
	else if(central_count .gt. 17) then
	!else if(central_count .gt. 24) then
	  bin_width(central_count) = 30
	  total_array%central_energy(central_count) = 1730.0+(central_count-18)*bin_width(central_count)
	  !bin_width(central_count) = 25
	  !total_array%central_energy(central_count) = 1903.0+(central_count-25)*bin_width(central_count)
	end if
	end if

	if(user_reaction .eq. 2) then
          if(central_count .lt. 15) then
	  bin_width(central_count) = 30!15
          total_array%central_energy(central_count) = 1625+bin_width(central_count)*(central_count-1)
          elseif(central_count .ge.15) then
	  bin_width(central_count) = 40
          total_array%central_energy(central_count) = 2050.0+bin_width(central_count)*(central_count-15)
          end if
	end if

	if(user_reaction .eq. 3 .or. user_reaction .eq. 4) then
	  bin_width(central_count) = 18.6
          total_array%central_energy(central_count) = 1096.7+bin_width(central_count)*(central_count-1)
	end if

	if(user_reaction .eq. 5) then
	  if(central_count .le. 8) then
	  bin_width(central_count) = 5
          total_array%central_energy(central_count) = 1492.5+bin_width(central_count)*(central_count-1)
          elseif(central_count .ge. 9) then
	  bin_width(central_count) = 20
          total_array%central_energy(central_count) = 1540.0+bin_width(central_count)*(central_count-9)
          elseif(central_count .ge. 26) then
	  bin_width(central_count) = 100
          total_array%central_energy(central_count) = 1940.+bin_width(central_count)*(central_count-26)
          end if
	end if

	if(user_reaction .eq. 6) then
	  bin_width(central_count) = 18.6
          total_array%central_energy(central_count) = 1617.6+bin_width(central_count)*(central_count-1)
		!1617.6
	end if

	!Reaction is KL p -> Pi+ Lambda HADRONIC    KLKS.dat 4 letter reaction ID: KLPL
	if(user_reaction .eq. 9) then
	  bin_width(central_count) = 20
          total_array%central_energy(central_count) = 1480+bin_width(central_count)*(central_count-1)
		!1617.6
	end if


	if(user_reaction .eq. 10) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1540+bin_width(central_count)*(central_count-1)
		!1617.6
	end if

	!Klong -> Kshort    KLKS.dat 4 letter reaction ID: KLKS
	if(user_reaction .eq. 11) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1450+bin_width(central_count)*(central_count-1)
		!1617.6
	end if

	!K-p Elastic    KLKS.dat 4 letter reaction ID: K-pE
	if(user_reaction .eq. 12) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1450+bin_width(central_count)*(central_count-1)
	  !Plab_loc = total_array%central_energy(central_count)**2
	  !Plab_loc = (Plab_loc - mass_nucleon**2 - mass_incident**2)/(2*mass_nucleon)
	  !Plab_loc = Plab_loc**2
	  !Plab_loc = Plab_loc - mass_incident**2
	  !Plab_loc = Sqrt(Plab_loc)
	  !print *,"Plab for bin",central_count," is:",Plab_loc
	end if

	!Reaction is K- p -> Sigma- Pi+ HADRONIC    KLKS.dat 4 letter reaction ID: K-S-
	if(user_reaction .eq. 13) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1540+bin_width(central_count)*(central_count-1)

	end if

	!Reaction is K- p -> Sigma+ Pi- HADRONIC    KLKS.dat 4 letter reaction ID: K-S+

	if(user_reaction .eq. 14) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1540+bin_width(central_count)*(central_count-1)

	end if

	!Reaction is K- p -> Sigma0 Pi0 HADRONIC    KLKS.dat 4 letter reaction ID: K-S0
	if(user_reaction .eq. 15) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1540+bin_width(central_count)*(central_count-1)

	end if

	!Haoran1st added
	!Reaction is K- p -> Kbr0 n HADRONIC    KLKS.dat 4 letter reaction ID: K-Kb
	if(user_reaction .eq. 16) then
	  bin_width(central_count) = 30
          total_array%central_energy(central_count) = 1450+bin_width(central_count)*(central_count-1)

	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	end do
	end subroutine



	!gathers from given files all data in specified range. call 1 is full range. call 2 is user specified range
	!returns it in both CM and lab frame
	subroutine  case_select() 
	use computer_specific
	use SetupReaction
	use constants
	implicit none

	integer	:: num_data_at_mmtm
	integer :: i,j
	double precision :: temp_mmtm,temp_mmtm_lab, momentum_val,W_cm
	double precision :: min_energy_cm,max_energy_cm !values in CM frame
	double precision :: min_energy_lab, max_energy_lab
	character :: obs*4  	   !obs=observable for given data
	character*12 :: author_name, page_num
	character*60 :: paper_date
	integer :: Counter_reactions
        character*4 :: PhotoReaction !determines if the reaction is KLambda or KSigma since both are in same file.

	character*6 :: DataFormatType
	integer :: Int_DataFormatType
	!***flag variable is needed to clear a file with no data. right now if I run reaction 2 which has extra
	! observables, then the extra files store data. Then when I run another reaction, those files will not clear and 
	!will store data into total array when it should not. flag will be set to 1 if data exists from raw datafile, then
	!it will clear files that do not have data.
	integer :: flag(1:3+num_observables) =0

	character*16 :: dump

	do i=1,3+num_observables
	  open (unit=1000+i,file=directory//subdirectory//observable_files(i),status="unknown")
	end do
	  open (11,file=directory//subdirectory//"paper_refs.dat",status="unknown")

	rewind (1000)
	rewind (11)

	read(1000,'(a16)')!headerline	
	!Haoran need and reactions here
	Do i=1,400000
		!Determines if data is in Igor format (3 sets per line) or Old format of 1 per line with File info.
		read(1000,'(a16,a6)',end=1000) dump,DataFormatType
	!if(i .gt. 100) print *,"Dataformattype:",dataformattype,"xx"
		if(DataFormatType .eq. "Sigma0")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date
			read(1000,*)obs
			Author_Name = "K-p_S0Pi0"
			PhotoReaction = "K-S0"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)

		else if(DataFormatType .eq. "KaonBa")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date

			read(1000,*)obs
			Author_Name = "K-p_Kbar0+n"
			PhotoReaction = "K-Kb"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)
		else if(DataFormatType .eq. "Sigma+")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date

			read(1000,*)obs
			Author_Name = "K-p_S+Pi-"
			PhotoReaction = "K-S+"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)
		else if(DataFormatType .eq. "Sigma-")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date

			read(1000,*)obs
			Author_Name = "K-p_S-Pi+"
			PhotoReaction = "K-S-"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)
		else if(DataFormatType .eq. "Lambda")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date

			read(1000,*)obs
			Author_Name = "K-p_LPi+"
			PhotoReaction = "K-pL"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)
		else if(DataFormatType .eq. "Kaon- ")Then
		        !print *,"XX",DataFormatType,"XX"
			!print *,"Old Data Type Found. Special Handling OKAY"
			read(1000,'(A)')Paper_date

			read(1000,*)obs
			Author_Name = "K^{-}p->K^{-}p"
			PhotoReaction = "K-pE"
			read(1000,*)momentum_val,num_data_at_mmtm

			Int_DataFormatType = 1   !Old data type flag
			W_CM = 2*mass_nucleon*Sqrt(mass_incident**2+momentum_val**2)
			W_CM = W_CM + mass_nucleon**2 + mass_incident**2
			W_CM = Sqrt(W_CM)
		else
		  backspace(1000)
			Int_DataFormatType = 2   !Igor File Type
		read(1000,100, END=1000), W_cm, &
					num_data_at_mmtm,PhotoReaction,&
					obs
		
		!W_cm=DSqrt(mass_nucleon**2+2*mass_nucleon*momentum_val) !in cm frame
 		!print *, momentum_val
		read(1000,*, END=1000) author_name, paper_date, page_num!, author line, may pick up a comma
		write(11,*) author_name,paper_date, page_num
		!write(20, '(a15,a15,a15)'),author_name, paper_date, page_num
		!this if statement goes to default case statement if energy is not in correct range
		end if

		if(W_cm.lt.Threshold_Energy_Wcm) THEN
		  !print *,"Found a momentum val below threshold, point(s) put into bin 1:",obs, momentum_val
		END if
		if(W_cm.gt. 2500.) then
		  !print *,"Found a momentum val above 2500MeV., point rejected from fits:", obs, momentum_val
		  obs = "NONE"
		end if


		select case(obs)
		
		case ("SGT1")
			flag(1)=1
			write (SGT1i, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,SGT1i)			  
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,SGT1i)		
		case ("SGT3")
			flag(2)=1
			write (SGT3i, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date			  
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,SGT3i)
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,SGT3i)		
		case ("SGT")
			flag(3)=1
			write (SGTi, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,SGT1i)		  
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,SGTi)
		case ("DSG")
			flag(4)=1
			write (DSGi, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,DSGi)
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,DSGi)
		case ("P")
			flag(5)=1
			write (Poli, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date		  
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,Poli)
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,Poli)
		case ("PDSG")
			flag(6)=1
			write (PDSGi, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date
			if(Int_DataFormatType .eq. 1) call old_data_points(num_data_at_mmtm,1000,PDSGi)
			if(Int_DataFormatType .eq. 2) call data_points(num_data_at_mmtm,1000,PDSGi)
		case ("F")
			flag(7)=1
			write (Fi, 150), W_cm,num_data_at_mmtm,PhotoReaction,obs,author_name, paper_date
			call data_points(num_data_at_mmtm,1000,Fi)
		case Default
			!print *, "Observable:", obs, "is not a valid observable value. Please fix input file."
			call data_points(num_data_at_mmtm,1000,tmp)
			!PRINT *,"nUM OF POINTS OF OBS IS:" , num_data_at_mmtm
		end select
	
	End DO
100	format(T2,F8.3,T12,I2,t43,a4,T48,a4)
150	format(F9.3,5x,I2,5x,a4,4x,a4,2x,a15,2x,a60)
1000	continue

	!close files -also delete files that have no data
	do j=1004,1000+num_observables+3 
	if(flag(j-1000).lt. .2) then
	  close(unit=j,status='DELETE')
		else
	    	  close (j)
		end if
	end do
	close(unit=1001)!sgt stays available for future call at end of program.
	close(unit=1002)!sgt stays available for future call at end of program.
	close(unit=1003)!sgt stays available for future call at end of program.
        close(unit=11)
	return
	end
	!this is for reading and writing to and from the observable files





	subroutine old_data_points(data_p, obs_fr,obs_fw)
	use constants
	!the characters a1,a2,a3 are the three real numbers in raw data file for actual data
	!points 1,4,7 are angles, 2,5,8 are diff cross section at angle, and 3,6,9 is diff_cx uncertainty
	!it allows me to read the whole line into separate variables
	integer data_p
	integer obs_fr
	integer n
	integer obs_fw
	double precision :: a1,a2,a3,b1,b2,b3,d1,d2,d3
	n=mod(data_p,3)  !gathers the number of datapoints on a line with <3 data_points
		!print *,"data_p,obsfr/fw:",data_p,obs_fr,obs_fw
	do i=1,data_p/3		!the do portion reads and writes to file the lines with 3 points
		read (obs_fr,*),a1,a2,a3
		read (obs_fr,*),b1,b2,b3
		read (obs_fr,*),d1,d2,d3
		!print *,"PI,",PI
		a1 = acos(a1)*180./Pi
		b1 = acos(b1)*180./Pi
		d1 = acos(d1)*180./Pi
		write (obs_fw,'(" ",3(f7.2,3x,f7.4,3x,f7.4,3x))')a1,a2,a3,b1,b2,b3,d1,d2,d3
	end do

	select case(n)	!this gets the last 1 or 2 points
	
	case(0)
		return		
	case (1)
		read (obs_fr,*) a1,a2,a3

		a1 = acos(a1)*180./Pi
		write (obs_fw,'(" ",3(f7.2,3x))') a1,a2,a3
	case (2)
		read (obs_fr,*) a1,a2,a3
		read (obs_fr,*) b1,b2,b3

		a1 = acos(a1)*180./Pi
		b1 = acos(b1)*180./Pi
		write (obs_fw,'(" ",6(f7.2,3x))') a1,a2,a3,b1,b2,b3
	case Default
		print *, "An unknown error has occurred!!!"
		stop
	end select

	
1000	return
	end 



	!this is for reading and writing to and from the observable files
	subroutine data_points(data_p, obs_fr,obs_fw)
	!the characters a1,a2,a3,b1... are the nine real numbers in raw data file for actual data
	!points 1,4,7 are angles, 2,5,8 are diff cross section at angle, and 3,6,9 is diff_cx uncertainty
	!it allows me to read the whole line into separate variables
	integer data_p
	integer obs_fr
	integer n
	integer obs_fw
	character*7 a1,a2,a3,b1,b2,b3,d1,d2,d3
		
	n=mod(data_p,3)  !gathers the number of datapoints on a line with <3 data_points
		
	do i=1,data_p/3		!the do portion reads and writes to file the lines with 3 points
		read (obs_fr,*, end=1000),a1,a2,a3,b1,b2,b3,d1,d2,d3
		write (obs_fw,'(" ",9(a7,3x))') &
			a1,a2,a3,b1,b2,b3,d1,d2,d3
	end do

	select case(n)	!this gets the last 1 or 2 points
	
	case(0)
		return		
	case (1)
		read (obs_fr,*,end=1000) a1,a2,a3
		write (obs_fw,'(" ",3(a7,3x))') a1,a2,a3
	case (2)
		read (obs_fr,*,end=1000) a1,a2,a3,b1,b2,b3
		write (obs_fw,'(" ",6(a7,3x))') a1,a2,a3,b1,b2,b3
	case Default
		print *, "An unknown error has occurred!!!"
		stop
	end select

1000	return
	end 
 
	!Duplicate of data_points, except it will store the data required to the internal variable total_array.
	!It is only called if the points fall within the energy range selected by the user.
	!first number is angle, second is OBS, and third is error OBS
	subroutine data_points_array(II,i,j, data_p, obs_fr,author_name, paper_date,energy_sgt,reaction_ID) !i=observable, j=bin, II=datapoint in bin
	use types

	!passed variables
	integer data_p
	integer obs_fr
	integer, intent(in) :: i,j !various passed variables. i j are passed for obs and bin correspondence.
	integer :: II
	character*12 :: author_name
	character*60 :: paper_date
	double precision :: energy_sgt
	character*4 :: reaction_ID
	!local variables
	integer :: n, counteri !n -case select, II dataarray counter,counteri=loop for lines with 3 data points
	real a1,a2,a3,b1,b2,b3,d1,d2,d3
	double precision :: fac !conversion factor-either 1 or 100.
	double precision :: radian= Pi/180.	
	n=mod(data_p,3)  !gathers the number of datapoints on a line with <3 data_points
	!if(i .eq. 4) fac = 100.
	fac = 1.
	

	do counteri=1,data_p/3		!the do portion reads and writes to file the lines with 3 points splitting data by observable type
		read (obs_fr,*, end=1000),a1,a2,a3,b1,b2,b3,d1,d2,d3

		!if(i .eq. 4 .and. j .eq. 15) print*,"Data bin 15 3 on line:",a1,a2,b1,b2,c1,c2

		if(a3.lt.0.0) a3=abs(a3)
		if(b3.lt.0.0) b3=abs(b3)
		if(d3.lt.0.0) d3=abs(d3)

		total_array%energy_values(i,j)%energies_datapoints(II,1)=a1
		total_array%energy_values(i,j)%energies_datapoints(II,2)=a2*fac
		total_array%energy_values(i,j)%energies_datapoints(II,3)=abs(a3)*fac+0.0005
		call Assign_author_name(i,j,II,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II)=reaction_ID
		total_array%energy_values(i,j)%energies_datapoints(II+1,1)=b1
		total_array%energy_values(i,j)%energies_datapoints(II+1,2)=b2*fac
		total_array%energy_values(i,j)%energies_datapoints(II+1,3)=abs(b3)*fac+0.0005
		call Assign_author_name(i,j,II+1,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II+1)=reaction_ID
		total_array%energy_values(i,j)%energies_datapoints(II+2,1)=d1
		total_array%energy_values(i,j)%energies_datapoints(II+2,2)=d2*fac
		total_array%energy_values(i,j)%energies_datapoints(II+2,3)=abs(d3)*fac+0.0005
		call Assign_author_name(i,j,II+2,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II+2)=reaction_ID


		!convert degree to radian- only needed for angle points
		total_array%energy_values(i,j)%energies_datapoints(II,1)=total_array%energy_values(i,j)%energies_datapoints(II,1)*radian
		total_array%energy_values(i,j)%energies_datapoints(II+1,1)=total_array%energy_values(i,j)%energies_datapoints(II+1,1)*radian
		total_array%energy_values(i,j)%energies_datapoints(II+2,1)=total_array%energy_values(i,j)%energies_datapoints(II+2,1)*radian

		II=II+3
	end do

	select case(n)	!this gets the last 1 or 2 points
	
	case(0)
		return !for lines where there are multiples of 3
	case (1)!case where we need to pick up 1 extra point
		read (obs_fr,*,blank="zero", end=1000) a1,a2,a3
		if(a3.lt.0) a3=abs(a3)

		!if(i .eq. 4 .and. j .eq. 15) print*,"Data bin 15 1extra:",a1,a2

		total_array%energy_values(i,j)%energies_datapoints(II,1)=a1
		total_array%energy_values(i,j)%energies_datapoints(II,2)=a2*fac
		total_array%energy_values(i,j)%energies_datapoints(II,3)=abs(a3)*fac+0.0005
		call Assign_author_name(i,j,II,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II)=reaction_ID
		total_array%energy_values(i,j)%energies_datapoints(II,1)=total_array%energy_values(i,j)%energies_datapoints(II,1)*radian!convert to radian
	!SGT points
		if(i .eq. 1) then
		  total_array%energy_values(i,j)%energies_datapoints(II,4) = energy_sgt
		end if
		if(i .eq. 2) then
		  total_array%energy_values(i,j)%energies_datapoints(II,4) = energy_sgt
		end if
		if(i .eq. 3) then
		  total_array%energy_values(i,j)%energies_datapoints(II,4) = energy_sgt
		end if

	if(i.eq. 11) then
	!print *,"Temp scaling down of E obs error bars for testing"
	!total_array%energy_values(i,j)%energies_datapoints(II,3)=total_array%energy_values(i,j)%energies_datapoints(II,3)/5.	
	end if
		II=II+1
	case (2)!case we need to pick up 2 extra points.
		read (obs_fr,*,end=1000) a1,a2,a3,b1,b2,b3

		!if(i .eq. 4 .and. j .eq. 15) print*,"Data bin 15 2extra:",a1,a2,b1,b2

	if (obs_fr .eq. 1005) then
	!print *, a1,a2,a3,b1,b2,b3
	end if
		if(a3.lt.0) a2=abs(a3)
		if(b3.lt.0) b2=abs(b3)
		
		total_array%energy_values(i,j)%energies_datapoints(II,1)=a1
		total_array%energy_values(i,j)%energies_datapoints(II,2)=a2*fac
		total_array%energy_values(i,j)%energies_datapoints(II,3)=abs(a3)*fac+0.0005
		call Assign_author_name(i,j,II,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II)=reaction_ID
		total_array%energy_values(i,j)%energies_datapoints(II+1,1)=b1
		total_array%energy_values(i,j)%energies_datapoints(II+1,2)=b2*fac
		total_array%energy_values(i,j)%energies_datapoints(II+1,3)=abs(b3)*fac+0.0005
		call Assign_author_name(i,j,II+1,author_name,paper_date)
		total_array%energy_values(i,j)%reactionID(II+1)=reaction_ID

		total_array%energy_values(i,j)%energies_datapoints(II,1)=total_array%energy_values(i,j)%energies_datapoints(II,1)*radian
		total_array%energy_values(i,j)%energies_datapoints(II+1,1)=total_array%energy_values(i,j)%energies_datapoints(II+1,1)*radian
	if(i.eq. 11) then
	!print *,"Temp scaling down of E obs error bars for testing"
	!print *,"Temp scaling down of E obs error bars for testing"	
	!total_array%energy_values(i,j)%energies_datapoints(II,3)=total_array%energy_values(i,j)%energies_datapoints(II,3)/5.
	!total_array%energy_values(i,j)%energies_datapoints(II+1,3)=total_array%energy_values(i,j)%energies_datapoints(II+1,3)/5.	
	end if
		II=II+2
	case Default
		print *, "Data points did not fit expected format of 3,6,or 9 points on a line!!!"
		stop
	end select

1000	return
	end

	!this subroutine creates table files for the given observables.
	!It is a little ancient in that most of its function was initially for testing.
	!However it still assigns energy values to the total array for each bin; As such it is required on each run.
	!It also counts the number of bins of data for user selected energy range.
	subroutine make_table(user_min_energy,user_max_energy)
	use computer_specific
	use constants
	use types
	implicit none
	integer :: bin_II=0 !counter for bins of width 25.
	integer :: i,j,k,num_data_points !i is the observable- j is the energy bin
	!min and max_energy here contain the total range(all bins) of energy values user has selected
	double precision :: user_min_energy, user_max_energy, energy
	double precision :: energy_temp
	integer :: II !II is the row within the given observable and energy bin. 
	!integer, Dimension(3:15,1000) :: counter  only needed for creating tables instead of array
        character*4 :: PhotoReaction !determines if the reaction is KLambda or KSigma
	character*18 :: author_name
	character*60 ::  paper_date
	character adjust_l*5,frame*3
	character*10 :: temp_holder

	open(30,file=directory//subdirectory//"Input/flag_status_small.dat",status="unknown")
	  rewind(30)
	  read(30,*)!Header Line
	  read(30,*)!Header Line
	  read(30,*) num_bins
	close(unit=30)

	open (unit=tablei+2001,file=directory//subdirectory//&     
				"DataTallyTable.dat",status="unknown") !needed for tables

	!print *, min_energy, max_energy
	do i=1,3+num_observables
	  open (unit=1000+i,file=directory//subdirectory//observable_files(i),status="unknown")
	  rewind(1000+i)
	end do
		
	do i=3,3+num_observables !cycle through observable files

	energy_temp = total_array%central_energy(1)
	bin_II=0
	  do j=1,num_bins !cycle through all energy bins
			read (i+1000,*, End=10000), energy 
                        !if(j.lt. 5)print *, "j and energy are:", energy, j,i,energy_temp
			backspace (i+1000)
			II=1
	!if(i .eq. 4)print *, "energy temp and width:", energy_temp,bin_width(j),j
!loop for k is arbitrary number big enough to make sure it finds an energy that is in user energy range. needed for files with no points in energy range and
!keeps j at 1 until it finds a value that is in the selected range.-->j=1 is first bin, regardless of starting energy
		  do k=1,100
		    if(energy.lt.energy_temp-bin_width(j)/2.) then
		    read (i+1000,*, End=10000) energy,num_data_points		
		    call data_points(num_data_points,1000+i,998)
		    read (i+1000,*, End=10000), energy 
		    backspace (i+1000)
		    end if
		    if(energy.ge.energy_temp-bin_width(j)/2.) exit 
		  end do

		!if(i.eq.6)print *, "before dowhile", energy
		if(energy.gt.user_max_energy) exit   !exits the second do loop of j meaning no additional energy in range found for observable
	!write(6, '(a34,3(f6.1,2x))'), "E, temp+width, and temp-width are:" ,&
			 !energy, energy_temp + bin_width(j)/2., energy_temp - bin_width(j)/2.
			do while (energy.lt.energy_temp + bin_width(j)/2. .and.&
				 energy .ge. energy_temp - bin_width(j)/2. &
				 .and. energy .le. user_max_energy)
				!compute #energies in bin range for files
				!if(energy .lt. 1600) print *, "various vals:", i,j,energy,energy_temp

				read (i+1000,55, End=10000) energy,num_data_points,&
						           PhotoReaction,temp_holder,author_name,paper_date
55	FOrmat(t2,f8.3,t15,i2,t22,a4,t30,a3,t39,a12,t53,a60)

!remove inconsistent data for gpep reaction
                                if(user_Reaction .eq. 1) then
				if((author_name(1:5) .eq. "SOEZU") .or. &
				   (author_name(1:5) .eq. "BACCI") .or. &
				   (author_name(1:5) .eq. "BLOOM") .or. &
				   !(author_name(1:4) .eq. "BOCK") .or. &
				   (author_name(1:6) .eq. "PSEUDO") .or. &
				   !(author_name(1:6) .eq. "PRELIM") .or. &
				   (author_name(1:5) .eq. "PREPO")) then
					call data_points(num_data_points,1000+i,998)
					!print *, "Author ", author_name(1:5), " being rejected"
					!print *, "Bin:" , j
					cycle
				end if
                                end if

                                if(user_Reaction .eq. 2) then
                                  if(PhotoReaction .ne. "K+L0") then
			           call data_points(num_data_points,1000+i,998)				   
				   cycle
				   end if
				  if ((author_name(1:4) .eq. "TRAN") .or. &
				    !(author_name(1:10) .eq. "PRELIM_CAS" .and. i .eq. 18) .or. &
				    !(author_name(1:10) .eq. "PRELIM_CAS" .and. i .eq. 19) .or. &
				    !(author_name(1:10) .eq. "PRELIM_WOL") .or. &
				    !(author_name(1:4) .eq. "HAAS") .or. &
				    (author_name(1:5) .eq. "HICKS") .or. &
				    (author_name(1:5) .eq. "GLAND"))&
				    !(author_name(1:5) .eq. "BOCKH") .or. &
				    !(author_name(1:5) .eq. "ALTHO")) 
					then
					call data_points(num_data_points,1000+i,998)
					!print *, "Author ", author_name(1:5), " being rejected"
					!print *, "Bin:" , j
					cycle
				  end if
				  !if((author_name(1:5) .eq. "BOCKH") .or. &
				   !(author_name(1:4) .eq. "THOM") .or. &
				   !(author_name(1:5) .eq. "FUJII") .or. &
				   !(author_name(1:4) .eq. "TRAN") .or. &
				   !(author_name(1:5) .eq. "GRILL")) then
			           !call data_points(num_data_points,1000+i,998)
				   !cycle
				  !end if
                                end if

                                if(user_Reaction .eq. 3) then
				  if(Photoreaction .ne. "PI+") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
				  if(Photoreaction .eq. 'PI0N') print *,"PI0N found for wrong reaction"
                                end if

                                if(user_Reaction .eq. 4) then
				  if(Photoreaction .ne. "PI0") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
				  if(Photoreaction .eq. 'PI0N') print *,"PI0N found for wrong reaction"
				!if((author_name(1:6) .eq. "HORNID") .or. &
				!    author_name(1:6) .eq. "HARTMA") then
				!   call data_points(num_data_points,1000+i,998)				   
				!   cycle
				!end if
                                end if

                                if(user_Reaction .eq. 7) then
				  if(Photoreaction .ne. "PI-") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
				  if(Photoreaction .eq. 'PI0N') print *,"PI0N found for wrong reaction"
                                end if

                                if(user_Reaction .eq. 8) then
				  if(Photoreaction .ne. "PI0N") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if

				!Klong p -> Pi+ Lambda
                                if(user_Reaction .eq. 9) then
				  if(Photoreaction .ne. "KLPL" .and. &
				     Photoreaction .ne. "K-pL") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if

				!Klong p -> Pi+ Sigma0
                                if(user_Reaction .eq. 10) then
				  if(Photoreaction .ne. "KLPS") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if

				!Klong p -> Kshort p
                                if(user_Reaction .eq. 11) then
				  if(Photoreaction .ne. "KLKS") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if


				!K-p Elastic Reaction
                                if(user_Reaction .eq. 12) then
				  if(Photoreaction .ne. "K-pE") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if

				!K-p to Sigma- Pi+
                                if(user_Reaction .eq. 13) then
				  if(Photoreaction .ne. "K-S-") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if



				!K-p to Sigma+ Pi-
                                if(user_Reaction .eq. 14) then
				  if(Photoreaction .ne. "K-S+") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if




				!K-p to Sigma+ Pi-
                                if(user_Reaction .eq. 15) then
				  if(Photoreaction .ne. "K-S0") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if

				!Haoran1st add
				
				!K-p to Kbr0 Neutron
                                if(user_Reaction .eq. 16) then
				  if(Photoreaction .ne. "K-Kb") then
				   call data_points(num_data_points,1000+i,998)				   
				   cycle
				  end if
                                end if




!Haoran1st added
				if(user_Reaction .gt. 16 .or. user_Reaction .eq. 6) stop "Reaction not setup in make_table subroutine"


		!if(j.lt. 4 .and. i .eq. 4) print *,energy,energy_temp-bin_width(j)/2.,energy_temp+bin_width(j)/2.,j
				!If the energy falls below mad bin energy, then assign datapoints to bin. Increment total number of points.
				if (energy.lt.(energy_temp + bin_width(j)/2.)) then
				call data_points_array(II,i,j,num_data_points,1000+i,author_name,paper_date,energy,Photoreaction)
				total_array%energy_num_points(i,j)= &
				total_array%energy_num_points(i,j) + num_data_points

				else
				backspace (i+1000)
				end if

			end do
		!Sets central point of the next bin for data binning
		   energy_temp = total_array%central_energy(j+1)
	end do




10000	continue
	if (j.gt.num_bins) then!this sets high limit on output values
		num_bins=j
	end if	
	end do

	!print *, "num bins is: " , num_bins, j
		!!!Haoran: originally was write(tablei+2001,'(a6,a64,a45)'),"W(GeV)"," DSG    P    G    F    H    E"

		write(tablei+2001,'(a6,a24,a45)'),"W(GeV)"," DSG    P    PDSG    F    H    E"
		write(tablei+2001,'(a15)'),"Bin width is 30"
		!write(tablei+2001,'(a76)'),"This reaction is K-p to Kbar0+n with DSG obs only,but I made 8 PDSG datapoints"
		do j=1,num_bins
		  	!write(tablei+2001,'(f10.7,2x)',advance='no'),pi
			write(tablei+2001,'(f5.0,2x)',advance='no'), total_array%central_energy(j)
!the do loop was  from 4 to 3+num_observables, to avoid 3 no-useful obs
	  	  do i=4,3+num_observables
			write (tablei+2001, '(i3,2x)',advance='no'), total_array%energy_num_points(i,j)
		  end do
			write(tablei+2001, '()')
		end do

		!do i=4,5   !observables
			  !write (tablei+2001,'(a7)'), observable_files(i)
			!do j=1,num_bins   !prints bins starting with first containing data, until last data filled bin
			  !write (tablei+2001,'(f8.2," ",f8.2, " ")',advance='no')&
			        !min_energy+bin_width*(j-1), min_energy+bin_width*j
			 ! write (tablei+2001,'((f8.3,3x))',advance='yes') total_array%central_energy(j)!,&
									!(total_array%central_energy(j)**2-mass_nucleon**2)&
									!/(2*mass_nucleon)
			!print *, "line 2834",total_array%central_energy(j)
			  !write (tablei+2001,'(i3,3x)',advance='no') total_array%energy_num_points(i,j) 
				!do k=1,total_array%energy_num_points(i,j) !k is the datapoints themselves
					!write (tablei+2001,'(f6.2,3x)',advance='no') &
					!      total_array%energy_values(i,j)%energies_datapoints(k,1)
					!write (tablei+2001,'(f9.4,3x)',advance='no') &
					!      total_array%energy_values(i,j)%energies_datapoints(k,2)
					!write (tablei+2001,'(f9.4,3x)') &
					!      total_array%energy_values(i,j)%energies_datapoints(k,3)
				!end do
				!write (tablei+2001,'(a7,3x,a5,3x,i2)') "NExt1  ", "obs: ",i
			!end do
		 !write (tablei+2001,'(" ")')	
		!end do
	close(unit=1001)
	close(unit=1002)	
	close(unit=1003)
	close(unit=1004) !close sgt and dsg but not delete
	do i=5,3+num_observables	
	close(unit=1000+i)
	end do
	!delete my temp file with data not needed.
	close(unit=tablei+2001)
	close(unit=tmp, status='delete')
	return

	end 
!program uses bin range calculated in CM frame.



!Routine that assigns a given data point an author name
!This is so when you graph observables in root, you can graph by author name and paper year.
	subroutine Assign_Author_name(int_i,int_j,int_II, author_name,paper_date)
	use constants
	use types

	implicit none

	integer :: int_i, int_j, int_II
	character :: author_name*(*)
	character :: paper_date*(*)

	!IF(INT_J .EQ. 6) PRINT *,AUTHOR_NAME
	!If paper date gets assigned '9999' because it is a new author, author name will get changed to "Not.Coded" at end of this routine
	total_array%energy_values(int_i,int_j)%data_points_author_name(int_II)=author_name
        if(user_reaction .eq. 1) then
	select case(author_name)

	case("PREPOST")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1967"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=2
	case("BACCI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=3
	case("BLOOM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=4
	case("ERBE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case("DELCOURT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=6
	case("HOLT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=7
	case("HEUSCH")
	  if(paper_date(4:5) .eq. '17') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	  end if
	  if(paper_date(4:5) .eq. '25') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=8
	  end if
	case("HONGOH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=9
	case("CHRIST")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1973"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=10
	case("BOOTH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1974"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=11
	case("VARTAPETYAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1980"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=12
	case("BOCK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1998"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=13
	case("AJAKA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1998"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=14
	case("KOUZNETSOV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1998"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=15
	case("HOMMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1988"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=16
	case("DYTMAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1995"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=17
	case("KRUSCHE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1995"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=18
	case("PRICE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1995"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=19
	case("DUGGER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2002"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=20
	case("AHRENS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2003"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=21
	case("CREDE")
	  if(paper_date(4:5) .eq. '94') then
	  total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2005"
	  total_array%energy_values(int_i,int_j)%author_number(int_II)=22
	  end if
	  if(paper_date(4:5) .eq. '80') then
	  total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	  total_array%energy_values(int_i,int_j)%author_number(int_II)=29
	  end if
	case("NAKABAYASHI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2006"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=23
	case("BARTALINI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=24
	case("ELSNER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=25
	case("SUMIHAMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=26
	case("WILLIAMS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=27
	case("SOEZUEER") !removed from fit because never published
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="unknown"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=28
	case("MCNICOLL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=30
	case("SIKORA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("AKONDI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=32
	case("PRELIM_HART")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=33
	case("PRELIM_MULLE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=34
	case("PSEUDODATA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=35
	case("SNDEROVICH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=36
	case("PRELIM_THIEL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=37
	case("PRELIM_WITTH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=38
	case("PRELIM_STRUB")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=40
	case Default
	print *, author_name, " author name found was not in list"
	end select
	end if

        if(user_reaction .eq. 2) then
	select case(author_name)
	case("DONOHO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1958"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BRODY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=2
	case("MCDANIEL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=3
	case("ANDERSON")
	  if(paper_date(4:4) .eq. '9') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1962"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=4
	  end if
	  if(paper_date(1:3) .eq. 'ISE') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1965"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=8
	  end if
	case("THOM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1963"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case("PECK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=6
	case("BORGIA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=7
	case("GRILLI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1965"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=9
	case("MORI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=10
	case("GROOM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1967"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=11
	case("ERBE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=12
	case("ABBHHM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=13
	case("FUJII")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=14
	case("DECAMP")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=15
	case("BLECKMANN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=16
	case("GEOING")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=17
	case("GOING")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=18
	case("GOEING")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=19
	case("FELLER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=20
	case("ALTHOFF")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=21
	case("HAAS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=22
	case("BOCKHORST")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1994"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=23
	case("TRAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1998"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=24
	case("GOERS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=25
	case("CARNAHAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2003"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=26
	case("ZEGERS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2003"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=27
	case("GLANDER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2004"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=28
	case("MCNABB")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2004"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=29
	case("SUMIHAMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2006"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=30
	case("KOHRI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2006"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("BRADFORD")
	  if(paper_date(4:5) .eq. '73') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2006"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=32
	  end if
	  if(paper_date(4:5) .eq. '75') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=33
	  end if
	case("LLERES")
	  if(paper_date(5:6) .eq. '31') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=34
	  end if
	  if(paper_date(5:6) .eq. '39') then
	   total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	   total_array%energy_values(int_i,int_j)%author_number(int_II)=36
	  end if
	case("HICKS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=35
	case("MCCRACKEN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=37
	case("DEY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=38
	case("JUDE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=39
	case("PSEUDODATA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="9999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=89
	case("PRELIM_WOLFO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=40
	case("PRELIM_PATER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=41
	case("PRELIM_CASEY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=46
	CASE dEFAULT
	PRINT *, "author:", AUTHOR_NAME, "NOT FOUND"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="9999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=99
	END SELECT
	end if

	if (user_reaction .eq. 3) then
	select case(author_name)
	case("DIXON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1958"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("TAYLOR")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("J.H.BOYDEN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1961"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("J.R.KILNER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1962"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("ALTHOFF")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1963"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("SMITH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1963"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BAZIN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1963"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("MCPHERSON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("LIU")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("FREYTAG")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1965"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("ADAMOVICH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("GRILLI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("FISHER")
	IF(paper_date(5:6) .eq. '70') then
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=18
	else if(paper_date(5:6) .eq. '72') then
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=22
	end if
	case("U.HAHN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=4
	case("KUZNETSOV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case("W.WALLROFF")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=6
	case("ALSPECTOR")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=14
	case("ZDARKO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=14
	case("KNIES")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1974"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=15
	case("ZABAEV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1975"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=16
	case("GANENKO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=18
	case("FUJII")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=22
	case("FUKUSHIMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=23
	case("ARAI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=25
	case("BUSSEY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1979"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=26
	case("ABRAHAMIAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1980"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=27
	case("P.HAMPE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1980"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=28
	case("EGAWA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1981"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=29
	case("BELYAEV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1984"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=30
	case("GETMAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1989"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("FISUM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=33
	case("MACCORMIC")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=34
	case("DUTZ")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=35
	case("KORKMAZ")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=38
	case("AJAKA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2000"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=40
	case("BECK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2000"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=60
	case("BLANPIED")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2001"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=65
	case("SCHMIT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2001"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=66
	case("BARTALINI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2002"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=67
	case("AHRENS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2004"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=70
	case("DUGGER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=75
	CASE dEFAULT
	PRINT *, "author:", AUTHOR_NAME, "NOT FOUND. PAPER DATE=9999; AUTHOR NUM=99,OBS=",INT_I
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="9999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=99
	END SELECT

 
	end if

	if (user_reaction .eq. 4) then
	select case(author_name) 
	case("STEIN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1959"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("JACKSON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("WORLOCK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("VASILKOV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("QUERZOLI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1960"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("TALMAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1962"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BERTANZA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1962"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("DRICKEY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("ALVAREZ")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1964"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("MALOY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1965"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("DE.STAEBLER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1965"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("ALTHOFF")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("G.L.HATCH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BRAUNSCHWEIG")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("WARD")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1967"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BLOOM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1967"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BACCI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1967"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("LUNDQUIST")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BUSCHORN","BUSCHHORN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("HAYAKAWA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("F.WOLVERTON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("GOVORKOV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1968"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("HITZEROTH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("DELCOURT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("BARBIELLINI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1969"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=1
	case("G.BOLOGNA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1970"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=22
	case("FISCHER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1971"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=23
	case("P.S.L.BOOTH","BOOTH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="ALL "
	total_array%energy_values(int_i,int_j)%author_number(int_II)=24
	case("PRENTICE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=25
	case("ALSPECTOR")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=26
	case("D.TRINES")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=4
	case("KABE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=7
	case("D.MENZE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=11
	case("ZDARKO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1972"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=27
	case("BECKS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1973"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=28
	case("GONCHAROV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1973"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=9
	case("HEMMI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1973"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=29
	case("TANAKA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1973"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=30
	case("KNIES")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1974"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=36
	case("HILGER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1974"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=12
	case("GENZEL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1975"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=13
	case("BREFELD")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1975"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=37
	case("BARTON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1975"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=38
	case("FELLER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=39
	case("DOUGAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=40
	case("ABRAHAMIAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=41
	case("DEREBCHINSKI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case("GANENKO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=42
	case("BUSSEY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=43
	case("BLUEM", "P.BLUEM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="ALL "
	total_array%energy_values(int_i,int_j)%author_number(int_II)=44
	case("ARAI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=45
	case("HUSMANN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=46
	case("HERR")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1977"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=50
	case("AVAKYAN", "R.O.AVAKYAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="ALL "
	total_array%energy_values(int_i,int_j)%author_number(int_II)=51
	case("GORBENKO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=52
	case("ZYBALOV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=8
	case("FUKUSHIMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=53
	case("SHUPE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1979"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=54
	case("YOSHIOKA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1980"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=55
	case("KATO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1980"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=56
	case("BRATASHEVSKI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="ALL "
	total_array%energy_values(int_i,int_j)%author_number(int_II)=57
	case("BELYAEV")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1983"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=58
	case("MAZZUCATO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1986"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=59
	case("ASATURYAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1986"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=60
	case("AGABABYAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1989"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=61
	case("BECK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1990"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=62
	case("BERGSTROM")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=63
	case("FUCHS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=74
	case("HAERTER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=75
	case("MACCORMIC")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1996"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=76
	case("SCHNEIDER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1997"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=13
	case("BOCK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1998"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=77
	case("KRUSCHE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=78
	case("REBREYEND")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=79
	case("SCHMIT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2001"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=80
	case("ADAMIAN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2001"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=81
	case("BLANPIED")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2001"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=82
	case("WIJESOORIYA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2002"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=86
	case("AHRENS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2004"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=87
	case("BARTHOLOMY")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2005"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=88
	case("BARTALINI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2005"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=89
	case("DUGGER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=2
	case("SUMIHAMA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2007"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=3
	case("ELSNER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=10
	case("SPARKS")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=6
	case("SCHUMANN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2010"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=14
	case("CREDE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2011"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=11
	case("THIEL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2012"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=90
	case("HORNIDGE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2013"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=91
	case("HARTMANN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=92
	case("GOTTSCHALL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=93
	case("SIKORA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=94
	case("LEUKEL")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="PC  "
	total_array%energy_values(int_i,int_j)%author_number(int_II)=95 
	case("ZENZ")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="19XX"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=96 
	case("S.ALMEHED")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="CONF"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=97 
	CASE dEFAULT
	!PRINT *, "author:", AUTHOR_NAME, "NOT FOUND. PAPER #=99. DATE=9999",INT_J,INT_I
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="9999"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=99
	END SELECT
	end if


	if (user_reaction .eq. 5) then
	select case(author_name)
	case("JAEGLE")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2008"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case("FANTINI")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2008"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=15
	case("WERTHMULLER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2014"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=35
	case("PRELIM_WITTH")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=45
	case("PRELIM_STRUB")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2015"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=65
	case default
	print *,"Author not spelled correctly:",author_name
	end select
	end if


	if (user_reaction .eq. 6) then
	select case(author_name)
	case("BANTAWA")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="2009"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=5
	case default
	print *,"Author not spelled correctly:",author_name
	end select
	end if


	if (user_reaction .eq. 9) then
	select case(author_name)
	case("K-p_LPi+")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="ALL"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=9
	case("KADYK")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1966"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31	
	case("BURKHARDT")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1975"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("CAMERON")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("CHO")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1976"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("ENGLER")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1978"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case("CORDEN")
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="1979"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=31
	case default
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K_Long Data"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=30
	end select
	end if


	if (user_reaction .eq. 10) then
	select case(author_name)

	case default
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K_Long Data"
	total_array%energy_values(int_i,int_j)%author_number(int_II)=15
	end select
	end if


	if (user_reaction .eq. 11) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="KLtoKS"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=12
	end select
	end if


	if (user_reaction .eq. 12) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K-PtoK-p"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=11
	end select
	end if


	if (user_reaction .eq. 13) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K-PS-Pi+"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=10
	end select
	end if



	if (user_reaction .eq. 14) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K-pS+Pi-"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=13
	end select
	end if



	if (user_reaction .eq. 15) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K-pS0Pi0"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=7
	end select
	end if

	!Haoran1st added
	if (user_reaction .eq. 16) then
	select case(author_name)

	case default
	!total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)="K-pKbr0n"
	total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II)=paper_date
	total_array%energy_values(int_i,int_j)%author_number(int_II)=7
	end select
	end if



	if(total_array%energy_values(int_i,int_j)%data_points_paper_date(int_II) .eq. '9999') then
	total_array%energy_values(int_i,int_j)%data_points_author_name(int_II)="Not.Coded"
	end if
	if(total_array%energy_values(int_i,int_j)%author_number(int_II) .EQ. 1) then
	total_array%energy_values(int_i,int_j)%data_points_author_name(int_II)="PRE-1970-DATA"
	end if
	return
	end
