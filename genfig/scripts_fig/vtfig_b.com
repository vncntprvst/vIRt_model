#!/bin/sh
#Plotting the membrane potentials of several I and F neurons,
#rastergrams of these populations and the "vibrissa angle".

suffix_a=dir_here/pb_v_rast_no_intra
gsyn_ext_a=0.5
suffix_b=dir_here/pb_v_rast_intra   
gsyn_ext_b=0.5
ipar=$1
irepeat=$2
tmax=3100.0

if [ $# -eq 0 ]
then
    col_file_a=$suffix_a/irt.col
    wsk_file_a=$suffix_a/irt.wsk
else
    col_file_a=$suffix_a/irt.col.$ipar.$irepeat
    wsk_file_a=$suffix_a/irt.wsk.$ipar.$irepeat
fi

if [ $# -eq 0 ]
then
    col_file_b=$suffix_b/irt.col
    wsk_file_b=$suffix_b/irt.wsk
else
    col_file_b=$suffix_b/irt.col.$ipar.$irepeat
    wsk_file_b=$suffix_b/irt.wsk.$ipar.$irepeat
fi
    
cat > line.xx <<EOF
2600.0 -84.0
2800.0 -84.0

3030.0 -20.0
3030.0  10.0
EOF

cat > linew.xx <<EOF
3030.0 100.0
3030.0 183.3333333
EOF

cat > frame.xx <<EOF
1950.0 0.03
1950.0 0.94

2025.0 0.03
2025.0 0.94
EOF

awk '{print $1 - 1250.0, '$gsyn_ext_a' * $2}' $suffix_a/irt.sex > $suffix_a/irt.sex.xx

# IRT
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1250.0 -0.01)) print $1 - 1250.0, $6}' $col_file_a > $suffix_a/irt.col.I1.1.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1250.0 -0.01)) print $1 - 1250.0, $7}' $col_file_a > $suffix_a/irt.col.I1.2.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1250.0 -0.01)) print $1 - 1250.0, $12}' $col_file_a > $suffix_a/irt.col.I2.1.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1250.0 -0.01)) print $1 - 1250.0, $13}' $col_file_a > $suffix_a/irt.col.I2.2.xx

# FN
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > -0.01)) print $1 - 1250.0, $18}' $col_file_a > $suffix_a/irt.col.F.1.xx

# wsk
awk '{
if (NF >= 2)
{
  if (NF == 2)
  {
    if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01))
    { 
      print $1 - 1250.0, $2 * 0.12
    }
  }
}
}' dir_here/pb_v_rast_no_intra/irt.wsk > dir_here/pb_v_rast_no_intra/irt.wsk.xx


awk '{print $1 - 1200.0, '$gsyn_ext_b' * $2}' $suffix_b/irt.sex > $suffix_b/irt.sex.xx

# IRT
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01)) print $1 - 1200.0, $6}' $col_file_b > $suffix_b/irt.col.I1.1.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01)) print $1 - 1200.0, $7}' $col_file_b > $suffix_b/irt.col.I1.2.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01)) print $1 - 1200.0, $12}' $col_file_b > $suffix_b/irt.col.I2.1.xx
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01)) print $1 - 1200.0, $13}' $col_file_b > $suffix_b/irt.col.I2.2.xx

# FN
awk '{if (($1 <= '$tmax' + 0.01) && ($1 > -0.01)) print $1 - 1200.0, $18}' $col_file_b > $suffix_b/irt.col.F.1.xx

# wsk
awk '{if (NF >= 2)
{
  if (($1 <= '$tmax' + 0.01) && ($1 > 1200.0 -0.01))
  { 
    print $1 - 1200.0, $2 * 0.12
  }
}
else
{
  print "   "
}
}' $wsk_file_b > $suffix_b/irt.wsk.xx
	
#gracebat \
xmgrace \
        -graph 0 $suffix_a/irt.sex.xx \
	-graph 1 $suffix_a/irt.sex.xx \
        -graph 2 $suffix_a/irt.col.I1.1.xx \
        -graph 3 $suffix_a/irt.col.I1.2.xx \
        -graph 4 $suffix_a/irt.col.I2.1.xx \
        -graph 5 $suffix_a/irt.col.I2.2.xx \
        -graph 6 $suffix_a/irt.col.F.1.xx line.xx \
        -graph 7 $suffix_a/irt.wsk.xx linew.xx \
        -graph 8 $suffix_b/irt.sex.xx \
	-graph 9 $suffix_b/irt.sex.xx \
        -graph 10 $suffix_b/irt.col.I1.1.xx \
        -graph 11 $suffix_b/irt.col.I1.2.xx \
        -graph 12 $suffix_b/irt.col.I2.1.xx \
        -graph 13 $suffix_b/irt.col.I2.2.xx \
        -graph 14 $suffix_b/irt.col.F.1.xx \
        -graph 15 $suffix_b/irt.wsk.xx \
        -hdevice EPS -p ../genfig/scripts_fig/vtfig_b.gr \
	-printfile $suffix_a/vtfig_b.eps

#linew.xx \

/bin/rm $suffix_a/irt.col.I*.xx $suffix_a/irt.col.F*.xx
#$suffix_a/irt.wsk.xx
/bin/rm line.xx frame.xx
/bin/rm $suffix_b/irt.col.I*.xx $suffix_b/irt.col.F*.xx $suffix_b/irt.wsk.xx
