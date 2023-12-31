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

"open ../gs_ctl/era5_spectrum.ctl"

"set mproj off"
"set t 1"
"set z 1"
"set lat -30 30"

"set grads off"
"set parea 1.2 10 1.2 4.15"
"set vrange -0.2 1"
"set ylint 0.2"
"set cmark 0"
"set cthick 16"
"set ccolor 1"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d mean(m(z=1),t=1,t=42)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d mean(m(z=2),t=1,t=42)"

"set xlabs ||||||||||||" 

"set parea 1.2 10 4.45 7.5"
"set vrange -0.2 1"
"set ylint 0.2"
"set cmark 0"
"set cthick 16"
"set ccolor 1"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d mean(z(z=1),t=1,t=42)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d mean(z(z=2),t=1,t=42)"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 0.7 lag [days]"
"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 7.8 Autocorrelation of z1, z2 and m1, m2"
"printim ./figure/auto_zm.png x2048 y1536"
"c"




"set grads off"
"set parea 1.2 10 1.2 4.15"
"set vrange -0.4 0.8"
"set ylint 0.2"
"set cmark 0"
"set cthick 16"
"set ccolor 1"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d mean(z1m2(z=1),t=1,t=42)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d mean(z2m1(z=1),t=1,t=42)"

"set xlabs ||||||||||||" 

"set parea 1.2 10 4.45 7.5"
"set vrange -0.4 0.8"
"set ylint 0.2"
"set cmark 0"
"set cthick 16"
"set ccolor 1"
"d const(z,0,-a)"
"set cmark 0"
"set cthick 16"
"set ccolor 4"
"d mean(zm(z=1),t=1,t=42)"
"set cmark 0"
"set cthick 16"
"set ccolor 2"
"d mean(zm(z=2),t=1,t=42)"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 0.7 lag [days]"
"set string 1 c 15 0"
"set strsiz 0.22"
"draw string 5.6 8.2 Cross-correlation" 
"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 7.8 z1m1:b, z2m2:r (upper) and z1m2:b, z2m1:r (lower)"
"printim ./figure/cross_zm.png x2048 y1536"
*"c"




