min=0
max=3
step=1

puntoc:3.dat

3.dat: a.out
	for i in $$(seq $(min) $(step) $(max)); do ./a.out $$i>$$i.dat;

a.out: puntoc.cpp
	g++ puntoc.cpp


clean:
	rm -rf *.txt *~ *.out  *.gp *.cpp# *.dat
