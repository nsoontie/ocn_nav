# Script to subset giops data as requested by a user on github
# Issue 759:
# I am looking for data from GIOPS Historical Daily dataset between August 18-22, 2018. The variables are vosaline (salinity) and votemper (temperature) for all depth levels. The bounding spatial limits are "max_range":"70.000, -124.1",
# "min_range":"69.3, -124.4"

# indices for subsetting were determined elsewhere

src='/home/soontiensn/remote2/hank/GIOPS/daily/201808'
out='/home/soontiensn/giops_tmps'
for i in {17..22}; do
  f=$src/giops_201808${i}00_024.nc
  echo $f
  basename=$(basename $f)
  ncks -d x,566,573 -d y,934,942 $f $out/$basename 
done
