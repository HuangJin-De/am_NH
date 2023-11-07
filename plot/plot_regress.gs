"reinit"
"set display color white"
"c"

"set lwid 13 16"
"set lwid 14 26"
"set lwid 15 10"
"set lwid 16 6"
"set annot 1 16"
"set strsiz 0.19"
"set xlopts 1 15 0.15"
"set ylopts 1 15 0.15"
"set clopts 1 15 0.15"
"set rgb 201 100 100 100 220"
"set grid on 3 201 12"
"set grid off"

"set font 11 file /home/der0318/.grads/Helvetica.ttf"

var1.1="z.1(z=1)"
var1.2="z.1(z=2)"
var1.3="m.1(z=1)"
var1.4="m.1(z=2)"

vname.1="z1"
vname.2="z2"
vname.3="m1"
vname.4="m2"

var2.1="z1.2(z=1)"
var2.2="z1.2(z=2)"
var2.3="m1.2(z=1)"
var2.4="m1.2(z=2)"

"open ../gs_ctl/era5_zm.ctl"
*"open ../gs_ctl/regress_zm.ctl"


n=1

j=1
while (j<=4)

i=1
while (i<=4)

"set z 1"
"set x 1"
"set t 1 5082"

"set grads off"
"set parea 1.2 10 1.2 7.5"
"set gxout scatter"
"set vrange -50 50"
"set vrange2 -50 50"
"set digsiz 0.1"
"set cmark 3"
"set rgb 124 170 170 170"
"set ccolor 124"
"d "var1.j";"var1.i""

"set t 1"
"set gxout grfill"
"define mm=mean("var1.i",t=1,t=5082)"
"define zm=mean("var1.j",t=1,t=5082)"
"d sum(("var1.i"-mm)*("var1.j"-zm),t=1,t=5082)/sum(pow("var1.j"-zm,2),t=1,t=5082)"
bb=sublin(result,3)
a=subwrd(bb,4)
say a
"d mm-zm*"a""
b=subwrd(result,4)
say b

"set xlab off"
"set ylab off"
"set x -50 50"
"set parea 1.2 10 1.2 7.5"
"set gxout line"
"set vrange -50 50"
"set ccolor 1"
"set cmark 0"
"set cthick 16"
"d lon/0.625*"a"-"b""

"set xlab on"
"set ylab on"

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 0.7 "vname.j""

"set string 1 c 15 90"
"set strsiz 0.18"
"draw string 0.5 4.35 "vname.i""

"set string 1 c 15 0"
"set strsiz 0.18"
"draw string 5.6 7.8 slope = "math_format('%5.3f',a)", intersection = "math_format('%5.3f',b)""

"printim ./figure/"vname.j"_"vname.i"_reg.png x2048 y1536"
"c"

n=n+1
say n

i=i+1
endwhile

j=j+1
endwhile


