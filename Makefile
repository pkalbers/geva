CC := g++

CFLAGS := -std=c++11 -Wno-unused-command-line-argument -O2 -march=native
CCLIBS  := -lz -lpthread

TARGET := geva_v1beta

SUBDIRS := $(shell find ./src/ -type d)
INCDIRS := $(addprefix -I, $(SUBDIRS))

SOURCES := $(shell find . -name '*.cpp')
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))


all: $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -flto $(INCDIRS) -c $< -o $@ $(CCLIBS)

$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) -flto $(CFLAGS) $(OBJECTS) $(CCLIBS)


.PHONY: clean

clean:
	@rm -f $(OBJECTS) $(TARGET)


