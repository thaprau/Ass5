galsim: galsim.c
	gcc -g -O3 galsim.c graphics/graphics.c graphics/graphics.h -o galsim -lm -lX11 -lpthread
run:
	./galsim 10000 input_data/ellipse_N_10000.gal 200 0.00001 0.252 0 2

clean:
	rm galsim