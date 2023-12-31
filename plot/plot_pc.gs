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

"open ../gs_ctl/era5_pc.ctl"

"set lat 10 80"
"set z 1"

"set grads off"
"set parea 1.2 10 2.2 6.5"
"set vrange -2 2"
"set ylint 0.5"
"set cmark 0"
"set ccolor 1"
"set cthick 15"
"d const(u,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d u(t=1)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d u(t=2)"
"set cmark 0"
"set cthick 16"
"set ccolor 11"
"d u(t=3)"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 1.7 latitude"

"set string 1 c 15 90"
"set strsiz 0.18"
"draw string 0.5 4.35 wind speed [m s`a-1`n]"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 6.7 EOF1 (red,52.4%), EOF2 (blue,27.7%), and EOF3 (cyan,10.7%) [DJFM]"

"printim ./figure/pc.png x2048 y1536"
*"c"


