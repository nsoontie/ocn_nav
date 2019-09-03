# Script to remove metadata from files

data_dir=/data/hdd/riops/riopsf/rotated/
atts_to_remove="institution source product_version contact history"
comment="DERIVED PRODUCT - created with script rotate_riops_forecast.sh"

for file in $data_dir/*.nc; do
    echo $file
    for att in $atts_to_remove; do
	ncatted -h -a $att,global,d,, $file
    done
    ncatted -h -a comment,global,m,c,"$comment" $file
done
