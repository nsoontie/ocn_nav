# Script to rotate salish sea model files
# Nancy Soontiens, March 2018

# Salish Sea specific variables and directories
GRID=/data/hdd/grids/salishsea/grid_angle_201702.nc
VARKEEP="nav_lon,nav_lat,time_counter,depth"
# Global attributes to remove
ATTS_TO_REMOVE="history"
# Directory structure
BASEDIR=/data/hdd/salishsea/model_Agrid/
subdirs=$(find $BASEDIR -type d -maxdepth 1 -mindepth 1)

for d in $subdirs; do
    echo $d 
    xvel=$d/*grid_U*
    yvel=$d/*grid_V*
    basename=$(basename $xvel .nc)
    dirname=$(dirname $xvel)
    subdir=$(basename $dirname)
    savedir=$BASEDIR/$subdir/
    savefile=$(echo $savedir/$basename | sed "s/grid_U_//")
    savefile=${savefile}_cardinal_velocity.nc
    echo $GRID $xvel $yvel $savefile
    bash rotate_velocity.sh \
	 $GRID \
	 $xvel \
	 $yvel \
	 $savefile \
    	 $VARKEEP \
	 "$ATTS_TO_REMOVE"
    # Clean up global attributes
    ncatted -h \
	    -a comment1,global,o,c,'velocities unstaggered to T grid' \
	    -a comment2,global,o,c,"created with script ${SCRIPTNAME}" \
	    -a description,global,o,c,'ocean velocities in cardinal directions' \
	    -a title,global,o,c,'ocean velocities in cardinal directions' \
	    $savefile    
done
