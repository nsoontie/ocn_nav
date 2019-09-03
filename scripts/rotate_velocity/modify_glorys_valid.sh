basedir=/data/hdd/glorys/v4/v4
for file in $basedir/*/*/*ice_cardinal_velocity.nc; do
echo $file
basename=$(basename $file)
out=tmp_$basename

cp $file $out
vars=$(ncdump -h $file)
var_list="air_ice_stress_east air_ice_stress_north ice_east_vel ice_north_vel iwind_stress_east iwind_stress_north"

for var in $var_list; do
echo $var
fill=$(grep -Po "${var}:_FillValue = \K.*(?=f ;)" <<< "${vars}")
min=$(grep -Po "${var}:valid_min = \K.*(?=f ;)" <<< "${vars}")
max=$(grep -Po "${var}:valid_max = \K.*(?=f ;)" <<< "${vars}")
echo $fill $min $max
ncap2 -O -h -s "where(${var} < ${min}) ${var}=${fill}" $out $out
ncap2 -O -h -s "where(${var} > ${max}) ${var}=${fill}" $out $out

done

mv $out ${file}
done
