
# This script can create MPAS background, analysis, or difference (analysis - background) plots
# The code comes from the py_scripts repo (https://github.com/ShawnMurdzek-NOAA/py_scripts/tree/main)

# HPC machine
machine='gaeaC6'

# MPAS input files
invariant_file=../../test_data/south3.5km.invariant.nc_L60_RAP
bgd_file=../../test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc
ana_file=../run/mpasout.nc

# Model levels to plot
lvls=(2 5 10)

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

fields=('theta'    'qv'   'qc'   'qi'   'qr'   'qs'   'qg'   'nr'  'ni'  'nc'  'cldfrac')
vmin_diff=(-2      -0.002 -0.001 -0.001 -0.001 -0.001 -0.001 -5000 -5000 -5000 -1)
vmax_diff=(2        0.002  0.001  0.001  0.001  0.001  0.001  5000  5000  5000  1)
for i in "${!fields[@]}"; do
  echo
  echo ${fields[i]}
  for l in ${lvls[@]}; do
    echo ${l}
    python ${script_dir}/plot_mpas_hcrsxn.py ${bgd_file} \
                                             ${invariant_file} \
                                             --outtag bgd \
                                             -l ${l} \
                                             -f ${fields[i]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

    python ${script_dir}/plot_mpas_hcrsxn.py ${ana_file} \
                                             ${invariant_file} \
                                             --outtag ana \
                                             -l ${l} \
                                             -f ${fields[i]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

    python ${script_dir}/plot_mpas_hcrsxn.py ${ana_file} \
                                             ${invariant_file} \
                                             --outtag diff \
                                             -l ${l} \
                                             -f ${fields[i]} \
                                             --file2 ${bgd_file} \
                                             -c bwr \
                                             --vmin ${vmin_diff[i]} \
                                             --vmax ${vmax_diff[i]} \
                                             --minlat ${minlat} \
                                             --maxlat ${maxlat} \
                                             --minlon ${minlon} \
                                             --maxlon ${maxlon}

  done
done
