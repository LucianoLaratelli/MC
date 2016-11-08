CC = g++
CFLAGS = -g -Wall
TARGET = grand
all : $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

clean:
	$(RM) $(TARGET) 
	$(RM) energies.dat
	$(RM) positions.xyz
	$(RM) qsts.dat
