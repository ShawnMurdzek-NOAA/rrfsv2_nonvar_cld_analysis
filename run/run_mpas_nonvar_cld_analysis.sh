#!/bin/sh

#SBATCH -M c6
#SBATCH -A bil-pmp
#SBATCH -t 00:10:00
#SBATCH --ntasks=48
#SBATCH --qos=debug

# Script to run all components of the nonvar cloud analysis for MPAS for a single time

# Input Parameters
# ================

machine='gaeaC6'
env_dir='../env/'
fix_dir='../fix/'
obs_dir="../../../../test_data/obs"
mpas_invariant_file="../../../../test_data/south3.5km.invariant.nc_L60_RAP"
mpasout_file="../../../../test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc"
dx=3500.0
valid=2024050813
proj_name='CONUS'

# Option to run each component
run_larccld=1
run_lightning=1
run_metarcld=1
run_mosaic=1
run_cloudanalysis=1


# Prepare to Run Program
# ======================

# Load environment
echo "Configuring environment"
source ${env_dir}/${machine}.env
module list
echo

# Get input data
cp ${mpas_invariant_file} invariant.nc
cp ${mpas_invariant_file} mesh.nc
cp ${mpasout_file} mpasout.nc
cp ${fix_dir}/prepobs_prep_RAP.bufrtable prepobs_prep.bufrtable
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.lgycld.tm00.bufr_d lgycld.bufr_d
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.lghtng.tm00.bufr_d lghtngbufr
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.prepbufr.tm00 prepbufr
cp ${obs_dir}/reflectivity/MergedReflectivityQC*${valid::8}-${valid:8:2}0* .
ls MergedReflectivityQC* > filelist_mrms

# Cleaning
clean_files=( 'lightning_interp.txt'  
	      'lightning_raw_1.txt'  
	      'nasa_larc_obs_interp.txt'  
	      'nasa_larc_obs_raw.txt'  
	      'processed_metar_obs_mpas.txt' )
for f in ${clean_files[@]}; do
  if [[ -f ${f} ]]; then
    rm ${f}
  fi
done


# Run larccld.fd
# ==============

if [[ ${run_larccld} -gt 0 ]]; then

cp ../larccld.fd/process_larccld.exe .

cat << EOF > namelist.nasalarc
 &setup
  analysis_time = ${valid},
  bufrfile='NASALaRCCloudInGSI_bufr.bufr',
  npts_rad=3,
  ioption = 2,
  userDX = ${dx},
  proj_name = '${proj_name}',
  debug=1,
 /
EOF

echo
echo "*****************************************************"
echo "Starting larccld.fd"
date
srun --export=ALL process_larccld.exe
date
echo "Completed larccld.fd"

fi


# Run lightning.fd
# ================

if [[ ${run_lightning} -gt 0 ]]; then

cp ../lightning.fd/process_Lightning.exe .

cat << EOF > namelist.lightning
 &setup
  analysis_time = ${valid},
  minute=00,
  trange_start=-10,
  trange_end=10,
  obs_type = "bufr",
  proj_name = '${proj_name}',
  debug=1
 /

EOF

echo
echo "*****************************************************"
echo "Starting lightning.fd"
date
srun --export=ALL process_Lightning.exe
date
echo "Completed lightning.fd"

fi


# Run metarcld.fd
# ===============

if [[ ${run_metarcld} -gt 0 ]]; then

cp ../metarcld.fd/process_metarcld.exe .

cat << EOF > namelist.metarcld
 &setup
  analysis_time = ${valid},
  prepbufrfile='prepbufr'
  twindin=0.5,
  metar_impact_radius=17,
  proj_name = '${proj_name}',
  debug=1,
 /
EOF

echo
echo "*****************************************************"
echo "Starting metarcld.fd"
date
srun --export=ALL process_metarcld.exe
date
echo "Completed metarcld.fd"

fi


# Run mosaic.fd
# =============

if [[ ${run_mosaic} -gt 0 ]]; then

cp ../mosaic.fd/process_NSSL_mosaic.exe .

cat << EOF > namelist.mosaic
   &setup
    tversion=1,
    analysis_time = ${valid},
    dataPath = './',
   /
EOF

cat << EOF > namelist.mosaic_netcdf
   &setup_netcdf
    output_netcdf = .true.,
    max_height = 11001.0,
    use_clear_air_type = .true.,
    precip_dbz_thresh = 10.0,
    clear_air_dbz_thresh = 5.0,
    clear_air_dbz_value = 0.0,
    precip_dbz_horiz_skip = 0,
    precip_dbz_vert_skip = 0,
    clear_air_dbz_horiz_skip = 0,
    clear_air_dbz_vert_skip = 0,
    remove_bdy = .false.,
   /
EOF

echo
echo "*****************************************************"
echo "Starting mosaic.fd"
date
srun --export=ALL process_NSSL_mosaic.exe
date
echo "Completed mosaic.fd"

fi


# Run cloudanalysis.fd
# ====================

if [[ ${run_cloudanalysis} -gt 0 ]]; then

cp ../cloudanalysis.fd/mpas_nonvarcldana.exe .

cat << EOF > gsiparm.anl
  
 &SETUP
  iyear=${valid::4},
  imonth=${valid:4:2},
  iday=${valid:6:2},
  ihour=${valid:8:2},
  iminute=00,
 /
 &RAPIDREFRESH_CLDSURF
   dfi_radar_latent_heat_time_period=20.0,
   metar_impact_radius=10.0,
   metar_impact_radius_lowCloud=4.0,
   l_pw_hgt_adjust=.true.,
   l_limit_pw_innov=.true.,
   max_innov_pct=0.1,
   l_cleanSnow_WarmTs=.true.,
   r_cleanSnow_WarmTs_threshold=5.0,
   l_conserve_thetaV=.true.,
   i_conserve_thetaV_iternum=3,
   l_cld_bld=.true.,
   l_numconc=.true.,
   cld_bld_hgt=1200.0,
   l_precip_clear_only=.false.,
   build_cloud_frac_p=0.50,
   clear_cloud_frac_p=0.10,
   iclean_hydro_withRef_allcol=1,
   i_gsdcldanal_type=6,
   i_gsdsfc_uselist=1,
   i_lightpcp=1,
   i_gsdqc=2,
   l_saturate_bkCloud=.true.,
   l_qnr_from_qr=.false.,
   n0_rain=100000000.0
   i_T_Q_adjust=1,
   l_rtma3d=.false.,
   i_precip_vertical_check=0,
   l_cld_uncertainty=.false.,
 /
EOF

echo
echo "*****************************************************"
echo "Starting cloudanalysis.fd"
date
srun --export=ALL mpas_nonvarcldana.exe
date
echo "Completed cloudanalysis.fd"

fi
