fox: fox.c
	@mpicc fox.c -o fox -lm -Wall

clean:
	@rm fox