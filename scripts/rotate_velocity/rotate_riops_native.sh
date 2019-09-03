# Script to rotate riops_native forecast files
# Nancy Soontiens, March 2018


# RIOPS_NATIVE specific variables and directories
GRID=/data/hdd/grids/riops_native/grid_angle.nc
VARKEEP="nav_lat,nav_lon"

# Global attributes to remove and comment to add
ATTS_TO_REMOVE="history"
COMMENT="DERIVED PRODUCT - created with script ${SCRIPTNAME}"

# Directory structure
OUT=/data/hdd/riops_native/riopsf/rotated/
FILES=/data/hdd/riops_native/riopsf/*.nc

for f in $FILES; do
    xvel=$f
    yvel=$f
    basename=$(basename $f .nc)
    savedir=$OUT
    if [[ ! -e $savedir ]] ; then
	mkdir -p $savedir
    fi
    # Ocean current first
    savefile=$savedir/${basename}_cardinal_velocity.nc
    echo $GRID $xvel $yvel $savefile
    bash rotate_velocity.sh \
	 $GRID \
	 $xvel \
	 $yvel \
	 $savefile \
	 $VARKEEP \
         "$ATTS_TO_REMOVE"
done
