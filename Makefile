CC 	= gcc
TARGET 	= ecc
SRC 	= ecc.c
FLAGS 	= -Wall -Wextra -O3 -march=native -lgmp
all:
	$(CC) -o $(TARGET) $(SRC) $(FLAGS)
clean:
	rm $(TARGET)
