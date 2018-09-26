ALL_TOOLS      += alpgen

ALL_TOOLS      += blackhat
blackhat_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/blackhat/0.9.9-giojec4/include
blackhat_EX_LIB := Ampl_eval BG BH BHcore CutPart Cut_wCI Cuteval Integrals Interface OLA RatPart Rateval Spinors assembly ratext
blackhat_EX_USE := qd

ALL_TOOLS      += boost_filesystem
boost_filesystem_EX_LIB := boost_filesystem
boost_filesystem_EX_USE := boost_system boost

ALL_TOOLS      += boost_header
boost_header_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-giojec4/include
boost_header_EX_USE := root_cxxdefaults

ALL_TOOLS      += boost_iostreams
boost_iostreams_EX_USE := boost

ALL_TOOLS      += boost
boost_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-giojec4/include
boost_EX_LIB := boost_thread boost_signals boost_date_time
boost_EX_USE := root_cxxdefaults sockets
boost_EX_FLAGS_CPPDEFINES  := -DBOOST_SPIRIT_THREADSAFE -DPHOENIX_THREADSAFE
boost_EX_FLAGS_CXXFLAGS  := -Wno-error=unused-variable

ALL_TOOLS      += boost_program_options
boost_program_options_EX_LIB := boost_program_options
boost_program_options_EX_USE := boost

ALL_TOOLS      += boost_python
boost_python_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-giojec4/include
boost_python_EX_LIB := boost_python
boost_python_EX_USE := root_cxxdefaults python

ALL_TOOLS      += boost_regex
boost_regex_EX_LIB := boost_regex
boost_regex_EX_USE := boost

ALL_TOOLS      += boost_serialization
boost_serialization_EX_LIB := boost_serialization
boost_serialization_EX_USE := boost

ALL_TOOLS      += boost_signals
boost_signals_EX_LIB := boost_signals
boost_signals_EX_USE := boost

ALL_TOOLS      += boost_system
boost_system_EX_LIB := boost_system
boost_system_EX_USE := boost

ALL_TOOLS      += boost_test
boost_test_EX_LIB := boost_unit_test_framework
boost_test_EX_USE := boost

ALL_TOOLS      += bz2lib
bz2lib_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/bz2lib/1.0.6/include
bz2lib_EX_LIB := bz2
bz2lib_EX_USE := root_cxxdefaults

ALL_TOOLS      += cascade_headers
cascade_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/cascade/2.2.04-giojec7/include
cascade_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += cascade
cascade_EX_LIB := cascade_merged
cascade_EX_USE := f77compiler cascade_headers

ALL_TOOLS      += castor_header
castor_header_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/castor/2.1.13.9/include /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/castor/2.1.13.9/include/shift
castor_header_EX_USE := root_cxxdefaults

ALL_TOOLS      += castor
castor_EX_LIB := shift castorrfio castorclient castorcommon
castor_EX_USE := castor_header libuuid

ALL_TOOLS      += ccache-ccompiler
ccache-ccompiler_EX_USE := gcc-ccompiler

ALL_TOOLS      += ccache-cxxcompiler
ccache-cxxcompiler_EX_USE := gcc-cxxcompiler

ALL_TOOLS      += ccache-f77compiler
ccache-f77compiler_EX_USE := gcc-f77compiler

ALL_TOOLS      += cgalimageio
cgalimageio_EX_LIB := CGAL_ImageIO
cgalimageio_EX_USE := zlib boost_system cgal

ALL_TOOLS      += cgal
cgal_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/cgal/4.2-giojec6/include
cgal_EX_LIB := CGAL_Core CGAL
cgal_EX_USE := root_cxxdefaults zlib boost_system gmp mpfr

ALL_TOOLS      += charybdis
charybdis_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/charybdis/1.003-giojec7/include
charybdis_EX_LIB := charybdis
charybdis_EX_USE := root_cxxdefaults f77compiler herwig pythia6

ALL_TOOLS      += classlib
classlib_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/classlib/3.1.3-ikhhed/include
classlib_EX_LIB := classlib
classlib_EX_USE := zlib bz2lib pcre openssl root_cxxdefaults
classlib_EX_FLAGS_CPPDEFINES  := -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS

ALL_TOOLS      += clhepheader
clhepheader_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/clhep/2.3.1.1-giojec4/include
clhepheader_EX_USE := root_cxxdefaults
clhepheader_EX_FLAGS_CXXFLAGS  := -Wno-error=unused-variable

ALL_TOOLS      += clhep
clhep_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/clhep/2.3.1.1-giojec4/include
clhep_EX_LIB := CLHEP
clhep_EX_USE := root_cxxdefaults
clhep_EX_FLAGS_CXXFLAGS  := -Wno-error=unused-variable

ALL_TOOLS      += cmsswdata
cmsswdata_EX_FLAGS_CMSSW_DATA_PACKAGE  := CalibCalorimetry/CaloMiscalibTools=V01-00-00 CalibTracker/SiPixelESProducers=V02-00-00 CalibTracker/SiStripDCS=V01-00-00 Calibration/Tools=V01-00-00 CondFormats/JetMETObjects=V01-00-03 Configuration/Generator=V01-00-01 DQM/PhysicsHWW=V01-00-00 DetectorDescription/Schema=V02-02-01 EventFilter/L1TRawToDigi=V01-00-00 FastSimulation/MaterialEffects=V04-02-08 FastSimulation/TrackingRecHitProducer=V01-00-02 Fireworks/Geometry=V07-05-01 GeneratorInterface/EvtGenInterface=V02-00-06 GeneratorInterface/ReggeGribovPartonMCInterface=V00-00-02 HLTrigger/JetMET=V01-00-00 L1Trigger/L1TCalorimeter=V01-00-16 L1Trigger/L1TGlobal=V00-00-06 L1Trigger/L1THGCal=V01-00-01 L1Trigger/L1TMuon=V01-00-07 L1Trigger/RPCTrigger=V00-15-00 MagneticField/Interpolation=V01-00-00 RecoBTag/CTagging=V01-00-02 RecoBTag/Combined=V01-00-02 RecoBTag/SecondaryVertex=V02-00-02 RecoBTag/SoftLepton=V01-00-01 RecoEgamma/ElectronIdentification=V01-00-07 RecoEgamma/PhotonIdentification=V01-00-03 RecoHI/HiJetAlgos=V01-00-01 RecoJets/JetProducers=V05-10-20 RecoLocalCalo/EcalDeadChannelRecoveryAlgos=V01-00-00 RecoMET/METPUSubtraction=V01-00-03 RecoMuon/MuonIdentification=V01-12-01 RecoParticleFlow/PFBlockProducer=V02-04-02 RecoParticleFlow/PFProducer=V15-00-00 RecoParticleFlow/PFTracking=V13-01-00 SLHCUpgradeSimulations/Geometry=V01-00-03 SimG4CMS/Calo=V03-01-00 SimG4CMS/Forward=V02-03-18 SimTracker/SiStripDigitizer=V01-00-00 Validation/Geometry=V00-07-00

ALL_TOOLS      += coral
ALL_SCRAM_PROJECTS += coral
coral_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/coral/CORAL_2_3_21-giojec12/include/LCG
coral_EX_USE := root_cxxdefaults
coral_ORDER := 98000

ALL_TOOLS      += cppunit
cppunit_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/cppunit/1.12.1/include
cppunit_EX_LIB := cppunit
cppunit_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += csctrackfinderemulation
csctrackfinderemulation_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/CSCTrackFinderEmulation/1.2/include
csctrackfinderemulation_EX_LIB := CSCTrackFinderEmulation

ALL_TOOLS      += curl
curl_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/curl/7.47.1/include
curl_EX_LIB := curl
curl_EX_USE := root_cxxdefaults

ALL_TOOLS      += cvs2git

ALL_TOOLS      += cython
cython_EX_USE := python

ALL_TOOLS      += das_client
das_client_EX_USE := python

ALL_TOOLS      += davix
davix_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/davix/0.6.4-giojec3/include/davix
davix_EX_LIB := davix
davix_EX_USE := boost_system openssl libxml2

ALL_TOOLS      += db6
db6_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/db6/6.0.30/include
db6_EX_LIB := db

ALL_TOOLS      += dbs-client

ALL_TOOLS      += dcap
dcap_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/dcap/2.47.8/include
dcap_EX_LIB := dcap
dcap_EX_USE := root_cxxdefaults

ALL_TOOLS      += distcc-ccompiler
distcc-ccompiler_EX_USE := gcc-ccompiler

ALL_TOOLS      += distcc-cxxcompiler
distcc-cxxcompiler_EX_USE := gcc-cxxcompiler

ALL_TOOLS      += distcc-f77compiler
distcc-f77compiler_EX_USE := gcc-f77compiler

ALL_TOOLS      += dmtcp

ALL_TOOLS      += doxygen

ALL_TOOLS      += dpm
dpm_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/dpm/1.8.0.1/include
dpm_EX_LIB := dpm
dpm_EX_LIBDIR := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/dpm/1.8.0.1/lib
dpm_EX_USE := root_cxxdefaults

ALL_TOOLS      += eigen
eigen_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2/include
eigen_EX_FLAGS_CPPDEFINES  := -DEIGEN_DONT_PARALLELIZE

ALL_TOOLS      += evtgen
evtgen_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/evtgen/1.6.0-giojec5/include
evtgen_EX_LIB := EvtGen EvtGenExternal
evtgen_EX_USE := hepmc pythia8 tauolapp photospp

ALL_TOOLS      += expat
expat_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/expat/2.1.0/include
expat_EX_LIB := expat
expat_EX_USE := root_cxxdefaults

ALL_TOOLS      += fastjet-contrib-archive
fastjet-contrib-archive_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fastjet-contrib/1.020/include
fastjet-contrib-archive_EX_LIB := EnergyCorrelator GenericSubtractor JetCleanser JetFFMoments Nsubjettiness ScJet SubjetCounting VariableR

ALL_TOOLS      += fastjet-contrib
fastjet-contrib_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fastjet-contrib/1.020/include
fastjet-contrib_EX_LIB := fastjetcontribfragile

ALL_TOOLS      += fastjet
fastjet_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fastjet/3.1.0/include
fastjet_EX_LIB := fastjetplugins fastjettools siscone siscone_spherical fastjet
fastjet_EX_USE := root_cxxdefaults

ALL_TOOLS      += fftjet
fftjet_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fftjet/1.5.0/include
fftjet_EX_LIB := fftjet
fftjet_EX_USE := root_cxxdefaults

ALL_TOOLS      += fftw3
fftw3_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fftw3/3.3.2/include
fftw3_EX_LIB := fftw3
fftw3_EX_USE := root_cxxdefaults

ALL_TOOLS      += file
file_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/file/5.18/include
file_EX_LIB := magic

ALL_TOOLS      += freetype
freetype_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/freetype/2.5.3/include
freetype_EX_LIB := freetype-cms
freetype_EX_USE := root_cxxdefaults

ALL_TOOLS      += frontier_client
frontier_client_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/frontier_client/2.8.20/include
frontier_client_EX_LIB := frontier_client
frontier_client_EX_USE := root_cxxdefaults zlib openssl expat python

ALL_TOOLS      += gcc-atomic
gcc-atomic_EX_LIB := atomic

ALL_TOOLS      += gcc-ccompiler
gcc-ccompiler_EX_FLAGS_CFLAGS  := -O2 -pthread
gcc-ccompiler_EX_FLAGS_CSHAREDOBJECTFLAGS  := -fPIC

ALL_TOOLS      += gcc-cxxcompiler
gcc-cxxcompiler_EX_FLAGS_CXXSHAREDFLAGS  := -shared -Wl,-E
gcc-cxxcompiler_EX_FLAGS_CXXSHAREDOBJECTFLAGS  := -fPIC
gcc-cxxcompiler_EX_FLAGS_LDFLAGS  := -Wl,-E -Wl,--hash-style=gnu
gcc-cxxcompiler_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -D_GNU_SOURCE
gcc-cxxcompiler_EX_FLAGS_CXXFLAGS  := -O2 -pthread -pipe -Werror=main -Werror=pointer-arith -Werror=overlength-strings -Wno-vla -Werror=overflow -std=c++14 -ftree-vectorize -Wstrict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 -fipa-pta -Wa,--compress-debug-sections -msse3 -felide-constructors -fmessage-length=0 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -Werror=delete-non-virtual-dtor -Werror=maybe-uninitialized -Werror=strict-aliasing -Werror=narrowing -Werror=uninitialized -Werror=unused-but-set-variable -Werror=reorder -Werror=unused-variable -Werror=conversion-null -Werror=return-local-addr -Werror=switch -fdiagnostics-show-option -Wno-unused-local-typedefs -Wno-attributes -Wno-psabi
gcc-cxxcompiler_EX_FLAGS_LD_UNIT  := -r -z muldefs

ALL_TOOLS      += gcc-f77compiler
gcc-f77compiler_EX_LIB := gfortran m
gcc-f77compiler_EX_FLAGS_FFLAGS  := -fno-second-underscore -Wunused -Wuninitialized -O2 -cpp
gcc-f77compiler_EX_FLAGS_FOPTIMISEDFLAGS  := -O2
gcc-f77compiler_EX_FLAGS_FSHAREDOBJECTFLAGS  := -fPIC

ALL_TOOLS      += gdb

ALL_TOOLS      += gdbm
gdbm_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gdbm/1.10/include
gdbm_EX_LIB := gdbm
gdbm_EX_USE := root_cxxdefaults

ALL_TOOLS      += geant4core
geant4core_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/geant4/10.02.p02-giojec5/include/Geant4
geant4core_EX_LIB := G4digits_hits G4error_propagation G4event G4geometry G4global G4graphics_reps G4intercoms G4interfaces G4materials G4parmodels G4particles G4persistency G4physicslists G4processes G4readout G4run G4tracking G4track G4analysis
geant4core_EX_USE := clhep root_cxxdefaults
geant4core_EX_FLAGS_CXXFLAGS  := -DG4MULTITHREADED -DG4USE_STD11 -ftls-model=global-dynamic -pthread
geant4core_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -DG4V9

ALL_TOOLS      += geant4data

ALL_TOOLS      += geant4
geant4_EX_USE := geant4core geant4vis xerces-c

ALL_TOOLS      += geant4-parfullcms

ALL_TOOLS      += geant4static
geant4static_EX_LIB := geant4-static
geant4static_EX_USE := clhep xerces-c
geant4static_EX_FLAGS_CXXFLAGS  := -DG4MULTITHREADED -DG4USE_STD11 -ftls-model=global-dynamic -pthread

ALL_TOOLS      += geant4vis
geant4vis_EX_LIB := G4FR G4modeling G4RayTracer G4Tree G4visHepRep G4vis_management G4visXXX G4VRML G4GMocren G4zlib
geant4vis_EX_USE := geant4core

ALL_TOOLS      += giflib
giflib_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/giflib/4.2.3/include
giflib_EX_LIB := gif
giflib_EX_USE := root_cxxdefaults

ALL_TOOLS      += git

ALL_TOOLS      += glibc

ALL_TOOLS      += glimpse

ALL_TOOLS      += gmake

ALL_TOOLS      += gmp
gmp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gmp-static/6.1.0/include
gmp_EX_LIB := gmp
gmp_EX_USE := mpfr

ALL_TOOLS      += gnuplot

ALL_TOOLS      += graphviz
graphviz_EX_USE := expat zlib libjpeg-turbo libpng

ALL_TOOLS      += gsl
gsl_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/2.2.1/include
gsl_EX_LIB := gsl gslcblas
gsl_EX_USE := root_cxxdefaults

ALL_TOOLS      += hdf5
hdf5_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/hdf5/1.8.17/include
hdf5_EX_LIB := hdf5

ALL_TOOLS      += hector
hector_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/hector/1.3.4_patch1-giojec14/include
hector_EX_LIB := Hector
hector_EX_USE := root_cxxdefaults

ALL_TOOLS      += hepmc_headers
hepmc_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/hepmc/2.06.07/include
hepmc_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += hepmc
hepmc_EX_LIB := HepMCfio HepMC
hepmc_EX_USE := hepmc_headers

ALL_TOOLS      += heppdt
heppdt_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/heppdt/3.03.00-giojec/include
heppdt_EX_LIB := HepPDT HepPID
heppdt_EX_USE := root_cxxdefaults

ALL_TOOLS      += herwig
herwig_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/herwig/6.521-giojec7/include
herwig_EX_LIB := herwig herwig_dummy
herwig_EX_USE := root_cxxdefaults f77compiler lhapdf photos

ALL_TOOLS      += herwigpp
herwigpp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/herwigpp/7.0.2-giojec4/include/Herwig
herwigpp_EX_USE := root_cxxdefaults lhapdf thepeg madgraph5amcatnlo openloops

ALL_TOOLS      += histfactory
histfactory_EX_LIB := HistFactory
histfactory_EX_USE := roofitcore roofit roostats rootcore roothistmatrix rootgpad rootxml rootfoam

ALL_TOOLS      += histogrammar
histogrammar_EX_USE := python py2-numpy py2-pandas

ALL_TOOLS      += igprof

ALL_TOOLS      += iwyu-cxxcompiler
iwyu-cxxcompiler_EX_USE := llvm-cxxcompiler

ALL_TOOLS      += jemalloc
jemalloc_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/jemalloc/4.2.1/include
jemalloc_EX_LIB := jemalloc
jemalloc_EX_USE := root_cxxdefaults

ALL_TOOLS      += jimmy_headers
jimmy_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/jimmy/4.2-giojec7/include
jimmy_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += jimmy
jimmy_EX_LIB := jimmy
jimmy_EX_USE := f77compiler herwig jimmy_headers

ALL_TOOLS      += ktjet
ktjet_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/ktjet/1.06-giojec7/include
ktjet_EX_LIB := KtEvent
ktjet_EX_USE := root_cxxdefaults
ktjet_EX_FLAGS_CPPDEFINES  := -DKTDOUBLEPRECISION

ALL_TOOLS      += lapack

ALL_TOOLS      += lcov

ALL_TOOLS      += lhapdf
lhapdf_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6-giojec7/include
lhapdf_EX_LIB := LHAPDF
lhapdf_EX_USE := yaml-cpp root_cxxdefaults

ALL_TOOLS      += libffi
libffi_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libffi/3.2.1/include
libffi_EX_LIB := ffi

ALL_TOOLS      += libhepml
libhepml_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libhepml/0.2.1/interface
libhepml_EX_LIB := hepml
libhepml_EX_USE := root_cxxdefaults

ALL_TOOLS      += libjpeg-turbo
libjpeg-turbo_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libjpeg-turbo/1.3.1/include
libjpeg-turbo_EX_LIB := jpeg turbojpeg
libjpeg-turbo_EX_USE := root_cxxdefaults

ALL_TOOLS      += libpng
libpng_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libpng/1.6.16/include
libpng_EX_LIB := png
libpng_EX_USE := root_cxxdefaults zlib

ALL_TOOLS      += libtiff
libtiff_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libtiff/4.0.3/include
libtiff_EX_LIB := tiff
libtiff_EX_USE := root_cxxdefaults libjpeg-turbo zlib

ALL_TOOLS      += libungif
libungif_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libungif/4.1.4/include
libungif_EX_LIB := ungif
libungif_EX_USE := root_cxxdefaults zlib

ALL_TOOLS      += libuuid
libuuid_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libuuid/2.22.2/include
libuuid_EX_LIB := uuid
libuuid_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += libxml2
libxml2_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libxml2/2.9.1/include/libxml2
libxml2_EX_LIB := xml2
libxml2_EX_USE := root_cxxdefaults

ALL_TOOLS      += libxslt
libxslt_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libxslt/1.1.28-giojec4/include/libxslt
libxslt_EX_LIB := xslt

ALL_TOOLS      += llvm-analyzer-ccompiler
llvm-analyzer-ccompiler_EX_USE := llvm-ccompiler

ALL_TOOLS      += llvm-analyzer-cxxcompiler
llvm-analyzer-cxxcompiler_EX_USE := llvm-cxxcompiler

ALL_TOOLS      += llvm-ccompiler
llvm-ccompiler_EX_USE := gcc-ccompiler

ALL_TOOLS      += llvm-cxxcompiler
llvm-cxxcompiler_EX_USE := gcc-cxxcompiler
llvm-cxxcompiler_EX_FLAGS_CXXFLAGS  := -Wno-c99-extensions -Wno-c++11-narrowing -D__STRICT_ANSI__ -Wno-unused-private-field -Wno-unknown-pragmas -Wno-unused-command-line-argument -ftemplate-depth=512 -Wno-error=potentially-evaluated-expression
llvm-cxxcompiler_EX_FLAGS_REM_CXXFLAGS  := -Wno-non-template-friend -Werror=format-contains-nul -Werror=maybe-uninitialized -Werror=unused-but-set-variable -Werror=return-local-addr -fipa-pta -frounding-math -mrecip -Wno-psabi

ALL_TOOLS      += llvm-f77compiler
llvm-f77compiler_EX_USE := gcc-f77compiler

ALL_TOOLS      += llvm
llvm_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.8.0-giojec4/include
llvm_EX_LIB := clang
llvm_EX_FLAGS_LDFLAGS  := -Wl,-undefined -Wl,suppress
llvm_EX_FLAGS_CXXFLAGS  := -D_DEBUG -D_GNU_SOURCE -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -O3 -fomit-frame-pointer -fPIC -Wno-enum-compare -Wno-strict-aliasing -fno-rtti

ALL_TOOLS      += madgraph5amcatnlo
madgraph5amcatnlo_EX_LIB := MadGraph5aMCatNLO
madgraph5amcatnlo_EX_USE := root_cxxdefaults

ALL_TOOLS      += mcdb
mcdb_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/mcdb/1.0.3/interface
mcdb_EX_LIB := mcdb
mcdb_EX_USE := root_cxxdefaults xerces-c

ALL_TOOLS      += mctester
mctester_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/mctester/1.25.0a-giojec14/include
mctester_EX_LIB := HEPEvent HepMCEvent MCTester
mctester_EX_USE := root_cxxdefaults root hepmc

ALL_TOOLS      += md5
md5_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/md5/1.0.0/include
md5_EX_LIB := cms-md5

ALL_TOOLS      += meschach
meschach_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/meschach/1.2.pCMS1/include
meschach_EX_LIB := meschach
meschach_EX_USE := root_cxxdefaults

ALL_TOOLS      += millepede
millepede_EX_USE := sockets pcre zlib

ALL_TOOLS      += mpfr
mpfr_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/mpfr-static/3.1.3/include
mpfr_EX_LIB := mpfr

ALL_TOOLS      += opengl
opengl_EX_LIB := GL GLU
opengl_EX_USE := x11

ALL_TOOLS      += openldap
openldap_EX_USE := openssl db6

ALL_TOOLS      += openloops

ALL_TOOLS      += openssl
openssl_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/openssl/1.0.2d/include
openssl_EX_LIB := ssl crypto
openssl_EX_USE := root_cxxdefaults

ALL_TOOLS      += oracle
oracle_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/oracle/12.1.0.2.0/include
oracle_EX_LIB := clntsh
oracle_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += oracleocci
oracleocci_EX_LIB := occi
oracleocci_EX_USE := oracle

ALL_TOOLS      += pacparser
pacparser_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pacparser/1.3.5/include
pacparser_EX_LIB := pacparser
pacparser_EX_USE := root_cxxdefaults

ALL_TOOLS      += pcre
pcre_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pcre/8.37/include
pcre_EX_LIB := pcre
pcre_EX_USE := root_cxxdefaults zlib bz2lib

ALL_TOOLS      += photos_headers
photos_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/photos/215.5/include
photos_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += photos
photos_EX_LIB := photos
photos_EX_USE := photos_headers f77compiler

ALL_TOOLS      += photospp
photospp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/photospp/3.61/include
photospp_EX_LIB := Photospp PhotosppHepMC PhotosppHEPEVT
photospp_EX_USE := hepmc

ALL_TOOLS      += professor

ALL_TOOLS      += protobuf
protobuf_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/protobuf/2.6.1/include
protobuf_EX_LIB := protobuf
protobuf_EX_USE := root_cxxdefaults

ALL_TOOLS      += py2-cjson
py2-cjson_EX_USE := python

ALL_TOOLS      += py2-cx-oracle
py2-cx-oracle_EX_USE := python oracle

ALL_TOOLS      += py2-cycler

ALL_TOOLS      += py2-dablooms

ALL_TOOLS      += py2-docopt

ALL_TOOLS      += py2-dxr
py2-dxr_EX_USE := python sqlite py2-futures py2-jinja py2-markupsafe py2-ordereddict py2-parsimonious py2-pygments py2-pysqlite

ALL_TOOLS      += py2-futures

ALL_TOOLS      += py2-h5py

ALL_TOOLS      += py2-ipython
py2-ipython_EX_USE := python

ALL_TOOLS      += py2-jinja

ALL_TOOLS      += py2-lint
py2-lint_EX_USE := python

ALL_TOOLS      += py2-markupsafe

ALL_TOOLS      += py2-matplotlib
py2-matplotlib_EX_USE := python zlib libpng py2-numpy py2-python-dateutil py2-cycler

ALL_TOOLS      += py2-networkx

ALL_TOOLS      += py2-numpy
py2-numpy_EX_USE := python zlib lapack

ALL_TOOLS      += py2-ordereddict

ALL_TOOLS      += py2-pandas
py2-pandas_EX_USE := python py2-numpy py2-python-dateutil py2-pytz

ALL_TOOLS      += py2-parsimonious

ALL_TOOLS      += py2-prettytable

ALL_TOOLS      += py2-pycurl
py2-pycurl_EX_USE := python

ALL_TOOLS      += py2-pygithub

ALL_TOOLS      += py2-pygments

ALL_TOOLS      += py2-pyparsing

ALL_TOOLS      += py2-pysqlite

ALL_TOOLS      += py2-python-dateutil
py2-python-dateutil_EX_USE := python

ALL_TOOLS      += py2-pytz
py2-pytz_EX_USE := python

ALL_TOOLS      += py2-pyyaml

ALL_TOOLS      += py2-requests

ALL_TOOLS      += py2-root_numpy

ALL_TOOLS      += py2-schema

ALL_TOOLS      += py2-scikit-learn
py2-scikit-learn_EX_USE := python py2-numpy py2-matplotlib py2-scipy 

ALL_TOOLS      += py2-scipy
py2-scipy_EX_USE := python py2-numpy lapack

ALL_TOOLS      += py2-setuptools
py2-setuptools_EX_USE := python

ALL_TOOLS      += py2-six

ALL_TOOLS      += py2-sqlalchemy
py2-sqlalchemy_EX_USE := python

ALL_TOOLS      += pyclang
pyclang_EX_USE := python

ALL_TOOLS      += pydata
pydata_EX_FLAGS_LDFLAGS  := $(PYDATA_BASE)/lib/pydata.o
pydata_EX_FLAGS_NO_RECURSIVE_EXPORT  := 1

ALL_TOOLS      += pyminuit2

ALL_TOOLS      += pyqt
pyqt_EX_USE := python qt sip

ALL_TOOLS      += pythia6_headers
pythia6_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pythia6/426/include
pythia6_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += pythia6
pythia6_EX_LIB := pythia6 pythia6_dummy pythia6_pdfdummy
pythia6_EX_USE := pythia6_headers f77compiler

ALL_TOOLS      += pythia8
pythia8_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pythia8/219-giojec5/include
pythia8_EX_LIB := pythia8
pythia8_EX_USE := root_cxxdefaults cxxcompiler hepmc lhapdf

ALL_TOOLS      += python-ldap
python-ldap_EX_USE := openssl openldap python

ALL_TOOLS      += python
python_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/python/2.7.11-giojec4/include/python2.7
python_EX_LIB := python2.7
python_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += qd
qd_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qd/2.3.13/include
qd_EX_LIB := qd_f_main qdmod qd

ALL_TOOLS      += qt3support
qt3support_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/Qt3Support
qt3support_EX_LIB := Qt3Support
qt3support_EX_USE := root_cxxdefaults qtbase
qt3support_EX_FLAGS_CPPDEFINES  := -DQT3_SUPPORT

ALL_TOOLS      += qtbase
qtbase_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/Qt /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/QtCore /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/QtXml
qtbase_EX_LIB := QtCore QtXml
qtbase_EX_USE := root_cxxdefaults zlib
qtbase_EX_FLAGS_CPPDEFINES  := -DQT_ALTERNATE_QTSMANIP -DQT_CLEAN_NAMESPACE -DQT_THREAD_SUPPORT

ALL_TOOLS      += qtdesigner
qtdesigner_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/QtDesigner
qtdesigner_EX_LIB := QtDesigner
qtdesigner_EX_USE := root_cxxdefaults qtbase qt

ALL_TOOLS      += qt
qt_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/QtOpenGL /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/qt/4.8.7/include/QtGui
qt_EX_LIB := QtOpenGL QtGui
qt_EX_USE := root_cxxdefaults qtbase qt3support x11 opengl

ALL_TOOLS      += rivet
rivet_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/rivet/2.5.2/include
rivet_EX_LIB := Rivet
rivet_EX_USE := python cython hepmc fastjet gsl yoda

ALL_TOOLS      += roofitcore
roofitcore_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.08-giojec3/include
roofitcore_EX_LIB := RooFitCore
roofitcore_EX_USE := rootcore roothistmatrix rootgpad rootminuit root_cxxdefaults

ALL_TOOLS      += roofit
roofit_EX_LIB := RooFit
roofit_EX_USE := roofitcore rootcore rootmath roothistmatrix

ALL_TOOLS      += roostats
roostats_EX_LIB := RooStats
roostats_EX_USE := roofitcore roofit rootcore roothistmatrix rootgpad

ALL_TOOLS      += rootcling
rootcling_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.08-giojec3/include
rootcling_EX_LIB := Core
rootcling_EX_USE := root_interface sockets pcre zlib

ALL_TOOLS      += rootcore
rootcore_EX_LIB := Tree Net
rootcore_EX_USE := rootmathcore rootthread

ALL_TOOLS      += root_cxxdefaults

ALL_TOOLS      += rooteg
rooteg_EX_LIB := EG
rooteg_EX_USE := rootgraphics

ALL_TOOLS      += rooteve
rooteve_EX_LIB := Eve
rooteve_EX_USE := rootgeompainter rootrgl

ALL_TOOLS      += rootfoam
rootfoam_EX_LIB := Foam
rootfoam_EX_USE := roothistmatrix

ALL_TOOLS      += rootgeom
rootgeom_EX_LIB := Geom
rootgeom_EX_USE := rootrio rootmathcore

ALL_TOOLS      += rootgeompainter
rootgeompainter_EX_LIB := GeomPainter
rootgeompainter_EX_USE := rootgeom rootgraphics

ALL_TOOLS      += rootgpad
rootgpad_EX_LIB := Gpad Graf
rootgpad_EX_USE := roothistmatrix

ALL_TOOLS      += rootgraphics
rootgraphics_EX_LIB := TreePlayer Graf3d Postscript
rootgraphics_EX_USE := rootgpad

ALL_TOOLS      += rootguihtml
rootguihtml_EX_LIB := GuiHtml
rootguihtml_EX_USE := rootinteractive

ALL_TOOLS      += roothistmatrix
roothistmatrix_EX_LIB := Hist Matrix
roothistmatrix_EX_USE := rootcore

ALL_TOOLS      += roothtml
roothtml_EX_LIB := Html
roothtml_EX_USE := rootgpad

ALL_TOOLS      += rootinteractive
rootinteractive_EX_LIB := Gui
rootinteractive_EX_USE := libjpeg-turbo libpng rootgpad rootrint

ALL_TOOLS      += root_interface
root_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.08-giojec3/include
root_interface_EX_USE := root_cxxdefaults

ALL_TOOLS      += rootmathcore
rootmathcore_EX_LIB := MathCore
rootmathcore_EX_USE := rootcling

ALL_TOOLS      += rootmath
rootmath_EX_LIB := GenVector MathMore
rootmath_EX_USE := rootcore gsl

ALL_TOOLS      += rootminuit2
rootminuit2_EX_LIB := Minuit2
rootminuit2_EX_USE := rootgpad

ALL_TOOLS      += rootminuit
rootminuit_EX_LIB := Minuit
rootminuit_EX_USE := rootgpad

ALL_TOOLS      += root
root_EX_USE := rootphysics
root_EX_FLAGS_NO_CAPABILITIES  := yes

ALL_TOOLS      += rootmlp
rootmlp_EX_LIB := MLP
rootmlp_EX_USE := rootgraphics

ALL_TOOLS      += rootphysics
rootphysics_EX_LIB := Physics
rootphysics_EX_USE := roothistmatrix

ALL_TOOLS      += rootpy
rootpy_EX_LIB := PyROOT
rootpy_EX_USE := rootgraphics

ALL_TOOLS      += rootrflx
rootrflx_EX_USE := root_interface rootcling
rootrflx_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DCMS_DICT_IMPL -D_REENTRANT -DGNUSOURCE -D__STRICT_ANSI__
rootrflx_EX_FLAGS_GENREFLEX_GCCXMLOPT  := -m64
rootrflx_EX_FLAGS_GENREFLEX_ARGS  := --deep

ALL_TOOLS      += rootrgl
rootrgl_EX_LIB := RGL
rootrgl_EX_USE := rootinteractive rootgraphics

ALL_TOOLS      += rootrint
rootrint_EX_LIB := Rint
rootrint_EX_USE := rootcling

ALL_TOOLS      += rootrio
rootrio_EX_LIB := RIO
rootrio_EX_USE := rootcling

ALL_TOOLS      += rootspectrum
rootspectrum_EX_LIB := Spectrum
rootspectrum_EX_USE := roothistmatrix

ALL_TOOLS      += rootthread
rootthread_EX_LIB := Thread
rootthread_EX_USE := rootrio

ALL_TOOLS      += roottmva
roottmva_EX_LIB := TMVA
roottmva_EX_USE := rootmlp rootminuit

ALL_TOOLS      += rootxmlio
rootxmlio_EX_LIB := XMLIO
rootxmlio_EX_USE := rootrio

ALL_TOOLS      += rootxml
rootxml_EX_LIB := XMLParser
rootxml_EX_USE := rootcore libxml2

ALL_TOOLS      += scons

ALL_TOOLS      += self
self_EX_INCLUDE := /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/combinetool/CMSSW_8_1_0/src /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/combinetool/CMSSW_8_1_0/include/LCG /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_1_0/src /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_1_0/include/LCG
self_EX_LIBDIR := /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/combinetool/CMSSW_8_1_0/lib/slc6_amd64_gcc530 /lustre/home/adrianodif/jpsiphi/2018/data_2018/analysis/utilities/combinetool/CMSSW_8_1_0/external/slc6_amd64_gcc530/lib /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_1_0/lib/slc6_amd64_gcc530 /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_1_0/external/slc6_amd64_gcc530/lib
self_EX_FLAGS_SYMLINK_DEPTH_CMSSW_SEARCH_PATH  := 2
self_EX_FLAGS_LLVM_ANALYZER  := llvm-analyzer
self_EX_FLAGS_SKIP_TOOLS_SYMLINK  := cxxcompiler ccompiler f77compiler gcc-cxxcompiler gcc-ccompiler gcc-f77compiler llvm-cxxcompiler llvm-ccompiler llvm-f77compiler llvm-analyzer-cxxcompiler llvm-analyzer-ccompiler icc-cxxcompiler icc-ccompiler icc-f77compiler x11 dpm
self_EX_FLAGS_DEFAULT_COMPILER  := gcc
self_EX_FLAGS_EXTERNAL_SYMLINK  := PATH LIBDIR CMSSW_SEARCH_PATH
self_EX_FLAGS_NO_EXTERNAL_RUNTIME  := LD_LIBRARY_PATH PATH CMSSW_SEARCH_PATH
self_EX_FLAGS_OVERRIDABLE_FLAGS  := CPPDEFINES CXXFLAGS FFLAGS CFLAGS CPPFLAGS LDFLAGS
self_ORDER := 20000
IS_PATCH:=

ALL_TOOLS      += sherpa
sherpa_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/sherpa/2.2.1-giojec5/include/SHERPA-MC
sherpa_EX_LIB := SherpaMain ToolsMath ToolsOrg
sherpa_EX_USE := root_cxxdefaults hepmc lhapdf qd blackhat fastjet sqlite openloops

ALL_TOOLS      += sigcpp
sigcpp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/sigcpp/2.6.2/include/sigc++-2.0
sigcpp_EX_LIB := sigc-2.0
sigcpp_EX_USE := root_cxxdefaults

ALL_TOOLS      += sip
sip_EX_USE := python

ALL_TOOLS      += sloccount

ALL_TOOLS      += sockets
sockets_EX_LIB := nsl crypt dl rt

ALL_TOOLS      += sqlite
sqlite_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/sqlite/3.12.2-giojec2/include
sqlite_EX_LIB := sqlite3
sqlite_EX_USE := root_cxxdefaults

ALL_TOOLS      += starlight
starlight_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/starlight/r193-giojec7/include
starlight_EX_LIB := Starlib
starlight_EX_USE := root_cxxdefaults clhep

ALL_TOOLS      += tauola_headers
tauola_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tauola/27.121.5/include
tauola_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += tauola
tauola_EX_LIB := pretauola tauola
tauola_EX_USE := f77compiler tauola_headers

ALL_TOOLS      += tauolapp
tauolapp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tauolapp/1.1.5-giojec7/include
tauolapp_EX_LIB := TauolaCxxInterface TauolaFortran TauolaTauSpinner
tauolapp_EX_USE := root_cxxdefaults hepmc f77compiler pythia8 lhapdf

ALL_TOOLS      += tbb
tbb_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tbb/44_20151115oss/include
tbb_EX_LIB := tbb
tbb_EX_USE := root_cxxdefaults

ALL_TOOLS      += tcmalloc_minimal
tcmalloc_minimal_EX_LIB := tcmalloc_minimal

ALL_TOOLS      += tcmalloc
tcmalloc_EX_LIB := tcmalloc

ALL_TOOLS      += thepeg
thepeg_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/thepeg/2.0.2-giojec3/include/ThePEG
thepeg_EX_LIB := ThePEG LesHouches
thepeg_EX_USE := root_cxxdefaults lhapdf gsl

ALL_TOOLS      += tinyxml
tinyxml_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tinyxml/2.5.3-giojec3/include
tinyxml_EX_LIB := tinyxml

ALL_TOOLS      += tkonlineswdb
tkonlineswdb_EX_LIB := DeviceDescriptions Fed9UDeviceFactory
tkonlineswdb_EX_USE := tkonlinesw oracle oracleocci

ALL_TOOLS      += tkonlinesw
tkonlinesw_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tkonlinesw/4.1.0-1-giojec7/include
tkonlinesw_EX_LIB := ICUtils Fed9UUtils
tkonlinesw_EX_USE := root_cxxdefaults xerces-c
tkonlinesw_EX_FLAGS_CXXFLAGS  := -DCMS_TK_64BITS

ALL_TOOLS      += toprex_headers
toprex_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/toprex/4.23/include
toprex_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += toprex
toprex_EX_LIB := toprex
toprex_EX_USE := toprex_headers f77compiler

ALL_TOOLS      += utm
utm_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/utm/r47119-xsd330-patch-giojec2/include
utm_EX_LIB := tmeventsetup tmtable tmxsd tmgrammar tmutil
utm_EX_USE := root_cxxdefaults

ALL_TOOLS      += valgrind
valgrind_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/valgrind/3.11.0/include
valgrind_EX_USE := root_cxxdefaults

ALL_TOOLS      += vdt_headers
vdt_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/vdt/v0.3.2-giojec2/include
vdt_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += vdt
vdt_EX_LIB := vdt
vdt_EX_USE := vdt_headers

ALL_TOOLS      += x11
x11_EX_USE := sockets

ALL_TOOLS      += xerces-c
xerces-c_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/xerces-c/3.1.3/include
xerces-c_EX_LIB := xerces-c
xerces-c_EX_USE := root_cxxdefaults

ALL_TOOLS      += xrootd
xrootd_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/xrootd/4.4.1/include/xrootd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/xrootd/4.4.1/include/xrootd/private
xrootd_EX_LIB := XrdUtils XrdClient
xrootd_EX_USE := root_cxxdefaults

ALL_TOOLS      += xz
xz_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/xz/5.2.1/include
xz_EX_LIB := lzma
xz_EX_USE := root_cxxdefaults

ALL_TOOLS      += yaml-cpp
yaml-cpp_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/yaml-cpp/0.5.1-giojec6/include
yaml-cpp_EX_LIB := yaml-cpp
yaml-cpp_EX_USE := boost

ALL_TOOLS      += yoda
yoda_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/yoda/1.6.5/include
yoda_EX_LIB := YODA
yoda_EX_USE := python cython root_cxxdefaults

ALL_TOOLS      += zlib
zlib_EX_INCLUDE := /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/zlib/1.2.8/include
zlib_EX_LIB := z
zlib_EX_USE := root_cxxdefaults

