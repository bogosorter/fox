fox: fox.c
	@mpicc fox.c -o fox -lm -Wall -O2

clean:
	@rm fox