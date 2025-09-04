#!/bin/sh

#SBATCH -M c6
#SBATCH -A bil-pmp
#SBATCH -t 00:10:00
#SBATCH --ntasks=10
#SBATCH --qos=debug

# Script to run mpas_nonvarcldana.exe for MPAS for a single time

# Input Parameters
# ================
machine='gaeaC6'
env_dir='../env/'
fix_dir='../fix/'
obs_dir="../../test_data/obs"
mpas_invariant_file="../../test_data/south3.5km.invariant.nc_L60_RAP"
mpasout_file="../../test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc"
dx=3500.0
valid=2024050813


# Run Program
# =========================

# Load environment
echo "Configuring environment"
source ${env_dir}/${machine}.env
module list
echo

# Get input data
cp ${mpas_invariant_file} invariant.nc
cp ${mpasout_file} mpasout.nc
cp ../metarcld.fd/mpas_metarcloud.bin .
cp ../mosaic.fd/RefInGSI3D.dat .
cp ../lightning.fd/LightningInMPAS.dat .

# Create namelist
# Options come from RRFSv1 version of rrfs-workflow (FV3-based)
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

# Clean old text files
#clean_files=( 'NASALaRC_cloud4mpas.bin' 'nasa_larc_obs_interp.txt' 'nasa_larc_obs_raw.txt' )
#for f in ${clean_files[@]}; do
#  if [[ -f ${f} ]]; then
#    rm ${f}
#  fi
#done

# Run program
date
echo
srun --export=ALL mpas_nonvarcldana.exe
echo
date
