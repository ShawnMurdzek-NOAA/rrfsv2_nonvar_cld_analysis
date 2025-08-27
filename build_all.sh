
# Build all nonvar cloud analysis components

# User-defined arguments

machine='gaeaC6'

# You shouldn't need to modify anything beyond this point

source ./env/${machine}.env
module list
echo

programs=( wps_map_utils
	   larccld.fd
	   lightning.fd
	   metarcld.fd
	   mosaic.fd )

top_dir=`pwd`
for p in ${programs[@]}; do
  cd ${p}
  echo
  echo ${p}
  make clean
  make
  cd ${top_dir}
done
