test:
	
	gcc -g -Wall -I/home/gsl/include -c main.c functions.c read_fn.c invert.c func2.c
	
	gcc -L/home/gsl/lib main.o functions.o read_fn.o invert.o func2.o -lgsl -lgslcblas -lumfpack -lm -o test

clean:
	rm test main.o func2.o functions.o read_fn.o invert.o print_complex_mat.txt print_mat_octave.txt 0_CALCULAT_FILE.txt PLOT_SAVED_BUS.csv PLOT_MY_ALGO.csv

clean_txt:
	rm write.txt

#gcc -g main.c read_fn.c functions.c -o test -lm