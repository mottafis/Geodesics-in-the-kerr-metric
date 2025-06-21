.PHONY: clear all
all: Burlisch-1dim.png Burlisch1dim.png burlisch-13d.png burlisch13d.png ergosfera.png rk4vsburlisch.png geodesics_comparison.pdf kerr_geodesics.pdf 
Burlisch-1dim.png: fotongeneralbur_-1.dat fotongeneral.dat fotongeneralplots.py
	python3 fotongeneralplots.py
Burlisch1dim.png: fotongeneralbur_1.dat fotongeneral.dat  fotongeneralplots.py
	python3 fotongeneralplots.py
burlisch-13d.png: fotongeneralbur_-1.dat fotongeneral.dat fotongeneralplots.py
	python3 fotongeneralplots.py
burlisch13d.png: fotongeneralbur_1.dat fotongeneral.dat fotongeneralplots.py
	python3 fotongeneralplots.py
ergosfera.png: fotongeneralbur_-1.dat fotongeneral.dat fotongeneralplots.py
	python3 fotongeneralplots.py
rk4vsburlisch.png: fotongeneralbur_-1.dat fotongeneralbur_1.dat fotongeneral.dat fotongeneralplots.py
	python3 fotongeneralplots.py
geodesics_comparison.pdf: direct_orbit.txt retrograde_orbit.txt plot_first_proff.py
	python3 plot_first_proff.py
kerr_geodesics.pdf: direct_orbit.txt retrograde_orbit.txt plot.py
	python3 plot.py 
fotongeneralbur_-1.dat fotongeneralbur_1.dat: con1y-1burlisch
	./con1y-1burlisch
con1y-1burlisch: con1y-1burlisch.cpp
	g++ con1y-1burlisch.cpp -o con1y-1burlisch 
fotongeneral.dat: fotongeneral
	./fotongeneral
fotongeneral: fotongeneral.cpp
	g++ fotongeneral.cpp -o fotongeneral
direct_orbit.txt retrograde_orbit.txt: Proyect
	./Proyect
proyect: Proyect.cpp
	g++ Proyect.cpp -o Proyect

