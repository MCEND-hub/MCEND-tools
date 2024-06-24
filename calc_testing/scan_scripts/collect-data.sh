#/bin/bash

# this collects misc MCEND data files for you
zipfiles='no'
dohhg='no'
dohhgx='no'
hhg=''
doacf='no'
Rn='1'
prof='no'
fftw_folder='hhg_spec_gen'

command_usage()
{
    echo ' '
    echo 'Hi! This collects MCEND data for you!'
    echo ' '
    echo '   collect_data.sh [storagename] -[additional options]'
    echo '                           -gz | -cz (gzip .chk files)'
    echo '                           -acf   (generate autocorrelation function spectrum)'
    echo '                           -hhg   (generate hhg spectrum for z direction)'
    echo '                           -hhgx  (generate hhg spectrum for x direction)'
    echo '                           -z     (take FFT of <z>)'
    echo '                           -x     (take FFT of <x>)'
    echo '                           -R     (take FFT of <R>)'
    echo '                           -p     (profile time with gprof)'
    echo '                           -h     (display usage)'
    echo '   '
    echo ' example: ./collect-data.sh LiH_421_dz_hhg -hhg -gz   '
    echo '   '
   exit
}


if [ $# -eq 0 ] || [ ${1} == -h ] || [ ${1} == --help ]; then
   command_usage
# else the JobName is taken as the first argument given.
else
   Rn=${1}
#   shift
   if [[ ${Rn} == [\-][a-z]* ]]; then
       echo "setting run directory to 1"
       Rn='1'
   else
       shift
   fi

# Parse command line options
# while there are more flags, check for each one if it is one of the following
   while [ "$1" != "" ]; do
       case $1 in
           -cz | --chkgzip | -gz )  zipfiles='zip'
                              ;;
           -acf | --autocf ) doacf='yes'
                              ;;
           -hhg | --hhg )    dohhg='yes'
                             ;;
           -hhgx | --hhgx )  dohhgx='yes'
                             ;;
           -z | --z )        direction='z'
                             ;;
           -R | --R )        direction='R'
                             ;;
           -x | --x )        direction='x'
                             ;;
           -h  | --help )    command_usage
                             exit
                             ;;
           -p  | --profile ) prof='yes'
                             ;;
           * )               echo
                             echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                             echo "* You put something in wrong, (${1}) isn't an option *"
                             echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                             echo
                             exit 1
       esac
       shift
   done

fi

Rtitle="Run"
#DATE=`date +%d.%m.%y`
DATE=`date +%m.%d.%y`
foldloc="${Rtitle}.${DATE}-${Rn}"
mkdir -p $foldloc

total_time=`grep "Wall time" mcend*out | awk '{print$4}'`
ttime_min=`echo "scale=2;${total_time}/60" | bc -l`

if [ $prof == 'yes' ]; then
    if [ -e gmon.out ]; then
        echo "Profiling with gprof requested..."
        gprof mcend gmon.out > "time-stats.out"
        echo "Done."
        mv "time-stats.out" $foldloc
        rm -f gmon.out
    else
        echo "Cannot profile without gmon.out file!!!"
    fi
fi

rm -f gmon.out

if [ $doacf == 'yes' ]; then
    # there is no missing $ for acf, it's feeding that string directly into the program
    ./${fftw_folder}/FT_transform acf &> acf.log
#    rm -f count.txt
elif [ $dohhg == 'yes' ]; then
    ./${fftw_folder}/FT_transform hhg &> hhg.log
elif [ $dohhgx == 'yes' ]; then
    ./${fftw_folder}/FT_transform hhgx &> hhg.log
else
    ./${fftw_folder}/FT_transform $direction  &> fft.log
#    rm -f count.txt
fi


mv *\.t      $foldloc
mv mcend*out $foldloc

if [ -e finalpsi ]; then
    cp finalpsi startpsi
    mv finalpsi $foldloc
fi

if [ -e finalrho ]; then
    mv finalrho $foldloc
fi

if [ -e ft_omega.dat ]; then
    mv ft_omega.dat $foldloc
fi


if [ $zipfiles == 'zip' ]; then
    tar -zcf chk-tar-bomb.tgz *.chk
    mv chk-tar-bomb.tgz ${foldloc}
    rm *\.chk
 else
    mv *\.chk $foldloc
fi

#python rdm_1_rho_check.py
mv ssqmat.diag $foldloc
mv *\.ij $foldloc
mv *\.Rt $foldloc
mv *\.log $foldloc
mv *\.list $foldloc
#mv 'rdm_mel_finalstep.pdf' $foldloc

rm -f *.o[0-9][0-9]*
echo "Results stored in ${foldloc}"
ttime_hr=`echo "scale=2;${ttime_min}/60" | bc -l`
echo "Calculation took ${total_time} s/${ttime_min} mins/${ttime_hr} hrs to run"
