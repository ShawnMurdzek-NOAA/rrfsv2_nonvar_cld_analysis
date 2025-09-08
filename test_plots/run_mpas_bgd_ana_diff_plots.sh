
# This script can create MPAS background, analysis, or difference (analysis - background) plots
# The code comes from the py_scripts repo (https://github.com/ShawnMurdzek-NOAA/py_scripts/tree/main)

# HPC machine
machine='gaeaC6'

# MPAS input files
invariant_file=../../test_data/south3.5km.invariant.nc_L60_RAP
bgd_file=../../test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc
ana_file=../run/mpasout.nc

# Model levels to plot
lvls=(9 10 11 12 13 14 15)

# Plotting domain
minlat=25
maxlat=42
minlon=-105
maxlon=-79

#-------------------------------------------------------------------------------
# You shouldn't need to change anything below this line

case ${machine} in

  'gaeaC6')
    module use /usw/conda/modulefiles
    module load miniforge
    conda activate /ncrc/home2/Shawn.S.Murdzek/.conda/envs/my_py
    export PYTHONPATH=$PYTHONPATH:/ncrc/home2/Shawn.S.Murdzek/src/
    script_dir=~/src/py_scripts/mpas
    ;;

  *)
    echo "Unknown machine option ${machine}"
    ;;

esac

fields=('qc'   'qi')
vmin=(0         0)
vmax=(0.00005  0.00005)

fields_diff=('theta'    'qv'   'qc'    'qi')
vmin_diff=(-0.25        -0.001 -0.0002 -0.0002)
vmax_diff=(0.25          0.001  0.0002  0.0002)

for l in ${lvls[@]}; do
  echo
  echo ${l}

  for j in "${!fields[@]}"; do
    python ${script_dir}/plot_mpas_hcrsxn.py ${bgd_file} \
                                             ${invariant_file} \
                                             --outtag bgd \
                                             -l ${l} \
                                             -f ${fields[j]} \
                                             --vmin ${vmin[j]} \
                                             --vmax ${vmax[j]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

    python ${script_dir}/plot_mpas_hcrsxn.py ${ana_file} \
                                             ${invariant_file} \
                                             --outtag ana \
                                             -l ${l} \
                                             -f ${fields[j]} \
                                             --vmin ${vmin[j]} \
                                             --vmax ${vmax[j]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

  done
  for j in "${!fields_diff[@]}"; do

    python ${script_dir}/plot_mpas_hcrsxn.py ${ana_file} \
                                             ${invariant_file} \
                                             --outtag diff \
                                             -l ${l} \
                                             -f ${fields_diff[j]} \
                                             --file2 ${bgd_file} \
                                             -c bwr \
                                             --vmin ${vmin_diff[j]} \
                                             --vmax ${vmax_diff[j]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

  done
done
