CC := g++

CFLAGS := -std=c++11 -Wno-unused-command-line-argument -O2 
CCLIBS  := -lz -lpthread

TARGET := geva

SUBDIRS := $(shell find ./src/ -type d)
INCDIRS := $(addprefix -I, $(SUBDIRS))

SOURCES := $(shell find . -name '*.cpp')
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))


all: $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCDIRS) -c $< -o $@ $(CCLIBS)

$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $(CFLAGS) $(OBJECTS) $(CCLIBS)


.PHONY: clean

clean:
	@rm -f $(OBJECTS) $(TARGET)
