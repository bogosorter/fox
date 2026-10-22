fox: fox.c
	@mpicc fox.c -o fox -lm -Wall -O2

fox_serial: fox_serial.c
	@mpicc fox_serial.c -o fox_serial -lm -Wall -O2

clean:
	@rm fox