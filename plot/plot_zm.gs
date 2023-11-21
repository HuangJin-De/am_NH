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

"open ../gs_ctl/era5_zm.ctl"

"set tlsupp year"

n=36 
"set t "n*121-120" "n*121""
"set z 1"

"set grads off"
"set parea 1.2 10 1.2 4.15"
"set vrange -40 40"
"set ylint 10"
"set cmark 0"
"set ccolor 1"
"set cthick 15"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d z(z=2)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d m(z=2)"

"set grads off"
"set parea 1.2 10 4.55 7.5"
"set vrange -40 40"
"set ylint 10"
"set cmark 0"
"set ccolor 1"
"set cthick 15"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d z(z=1)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d m(z=1)"


"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 0.7 time"

"set string 1 c 15 90"
"set strsiz 0.18"
"draw string 0.5 2.675 index"

"set string 1 c 15 90"
"set strsiz 0.18"
"draw string 0.5 6.025 index"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 7.8 Z1, M1 (upper) and Z2, M2 (lower) [Z:red, M:blue]"

"printim ./figure/zm_series.png x2048 y1536"
*"c"


