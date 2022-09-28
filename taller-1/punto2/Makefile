min=0
max=3
step=1

all: cfigura.pdf

cfigura.pdf: cfigura.gp lambda3.dat
	gnuplot cfigura.gp


cfigura.gp: lambda3.dat
	echo "set term pdf" >>$@
	echo "set out 'cfigura.pdf'" >>$@
	echo "set title 'Modos normales' " >>$@
	echo "set xlabel 'r'" >>$@
	echo "set ylabel 'R(r,lambda)'" >>$@
	echo "plot 'lambda0.dat' w l t 'R0(r,2.40527)',  'lambda1.dat' w l t 'R0(r,5.52245)','lambda2.dat' w l t 'R0(r,8.65954)', 'lambda3.dat' w l t 'R0(r,11.8022)'">>$@

lambda3.dat: a.out
	./a.out 1 0 >$llambda0.dat;
	./a.out 1 1 >$llambda1.dat;
	./a.out 1 2 >$llambda2.dat;
	./a.out 1 3 >$llambda3.dat;

a.out: puntoc.cpp
	g++ puntoc.cpp


clean:
	rm -rf *.txt *~ *.out  *.gp *.cpp# *.dat
