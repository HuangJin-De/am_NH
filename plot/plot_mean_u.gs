"reinit"
"set display color white"
"c"

"set lwid 13 16"
"set lwid 14 26"
"set lwid 15 10"
"set lwid 16 6"
"set annot 1 15"
"set strsiz 0.19"
"set xlopts 1 15 0.15"
"set ylopts 1 15 0.15"
"set clopts 1 15 0.15"
"set rgb 201 100 100 100 220"
"set grid on 3 201 12"
*"set grid off"

"set font 11 file /home/der0318/.grads/Helvetica.ttf"

"open ../gs_ctl/era5_daily.ctl"

"set lat 10 80"
"set lev 1000 100"

"set grads off"
"set parea 1.2 10 1.2 6.5"
"set ylevs 1000 925 850 700 500 250 100"
"color -80 80 5 -kind blue->cyan->white->yellow->red"
"d mean(emc,t=1,t=5082)"

"xcbar 10.15 10.25 2.2 5.5 -edge triangle -fs 2 -ft 15 -fh 0.15 -fw 0.15"

"define um=mean(u,t=1,t=5082)"

"set gxout contour"
"set clevs 5 10 15 20 25 30"
"set cthick 16"
"set clab masked"
"d um"

"set gxout contour"
"set clevs 0"
"set cthick 15"
"set rgb 152 180 180 180"
"set ccolor 152"
"set clab masked"
"d um"

"set gxout contour"
"set clevs -10 -5"
"set cthick 16"
"set ccolor 1"
"set cstyle 2"
"set clab masked"
"d um"

"close 1"

"open ../gs_ctl/int.ctl"

"set z 1"
"set t 1"
"set lat 10 80"

"set xlabs |||||||"
"set parea 1.2 10 6.6 8"
"set ylpos 0 l"
"set ylopts 1 0.15 15"
"set vrange -2 3"
"set ylevs -1 0 1 2"
"set cmark 0"
"set ccolor 1"
"set cthick 15"
"d const(u,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d mean(mean(emc*p,z=1,z=121),t=1,t=42)"

"set ylpos 0 r"
"set ylopts 11 0.15 15"
"set vrange -5 20"
"set ylevs  0 5 10 15"
"set cmark 0"
"set cthick 16"
"set ccolor 11"
"d mean(mean(u,z=1,z=121),t=1,t=42)"


"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 0.7 latitude"

"set string 1 c 15 90"
"set strsiz 0.18"
"draw string 0.5 4.35 pressure [hPa]"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 8.2 Mean zonal wind and eddy flux [DJFM]"

"printim ./figure/mean.png x2048 y1536"
"c"


