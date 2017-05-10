CC = gcc  
XX = g++  
CFLAGS = -Wall -O –g  
TARGET = ./SpectraMatch  
  
SOURCES = $(wildcard *.cpp)  
OBJS = $(patsubst %.cpp, %.o, $(SOURCES))  
  
$(TARGET):$(OBJS)  
	$(XX) $(OBJS) -lpthread -o $(TARGET)  
	chmod a+x $(TARGET)  
  
clean:  
	rm -rf *.o ./SpectraMatch
