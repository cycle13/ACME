#!/bin/bash 
set -e
HOMME=~/codes/acme-dev/components/homme
wdir=~/scratch1/sweqx
MACH=$HOMME/cmake/machineFiles/redsky.cmake
NCPU=0


#
#  problem setup
#
tstep=30 
nu=1.2e-6 #3.5e-8 for hv_scaling=3.2, 1.2e-6 for hv_scaling=4.0
hvscaling=4
test_case=swtc2
name=${test_case}-tensor
mesh=~/codes/mapping/grids/mountain_10_x2.g

#
# run/build directoryies
#
input=$HOMME/test/sw_conservative
builddir=$wdir/bld
rundir=$wdir/$name
mkdir -p $builddir
mkdir -p $rundir


cd $builddir
if ! ([ -e  CMakeCache.txt ]) then
  rm -rf CMakeFiles CMakeCache.txt
  cmake -C $MACH -DSWEQX_PLEV=1  -DSWEQX_NP=4 $HOMME
  exit
fi


#make clean 
make -j4 sweqx
exe=$builddir/src/sweqx/sweqx


mkdir -p $rundir/movies
cd $rundir
rsync $mesh .
let sfreq=3600
sfreq=`echo "$sfreq / $tstep" | bc`

sed s/tstep.\*/"tstep = $tstep"/  $input/swtc2-tensor-hv.nl |\
sed s/hypervis_scaling.\*/"hypervis_scaling = $hvscaling"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
    > input.nl

mpirun -np $NCPU  $exe  < input.nl | tee  sweq.out

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out
mv -f movies/swtc21.nc movies/$name.nc



