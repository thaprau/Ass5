galsim: galsim.c
	gcc -g galsim.c graphics/graphics.c graphics/graphics.h -o galsim -lm -lX11 -lpthread
run:
	./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0.252 0 16

clean:
	rm galsim