galsim: galsim.c
	gcc -g galsim.c graphics/graphics.c graphics/graphics.h -o galsim -lm -lX11 -lpthread
run:
	./galsim 500 input_data/ellipse_N_00500.gal 200 0.00001 0 0 8

clean:
	rm galsim