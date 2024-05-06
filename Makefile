SRC=src/
CAMERA_SRC=$(SRC)camera/
FILTER_SRC=$(SRC)filter/
SERIAL_SRC=$(SRC)serial/
BUILD=build/
INSTALL=bin/
LIBS=$(INSTALL)libs/

COMMON_FLAGS = -Wall -Wextra -g -O0
CFLAGS = $(COMMON_FLAGS) -std=gnu99  
CPPFLAGS = $(COMMON_FLAGS) -fPIC
LFLAGS = -L $(BUILD) -lm -pthread
INCLUDES_FLAGS = -I $(CAMERA_SRC) -I $(FILTER_SRC) -I $(SERIAL_SRC) -I $(SRC)

MAIN_C = $(SRC)main.c
MAIN_H = $(SRC)config.h
MAIN_O = $(patsubst $(SRC)%.c,$(BUILD)%-c.o,$(MAIN_C))

LIBCAMERA_CPP = $(CAMERA_SRC)camera_interface_C.cpp $(CAMERA_SRC)camera_vision.cpp 
LIBCAMERA_HPP = $(CAMERA_SRC)camera_vision.hpp $(CAMERA_SRC)camera.h 
LIBCAMERA_CPP_O = $(patsubst $(CAMERA_SRC)%.cpp,$(BUILD)%-cpp.o,$(LIBCAMERA_CPP))

LIBCAMERA_C = $(CAMERA_SRC)aruko_detect.c
LIBCAMERA_H = $(CAMERA_SRC)aruko_detect.h
LIBCAMERA_C_O = $(patsubst $(CAMERA_SRC)%.c,$(BUILD)%-c.o,$(LIBCAMERA_C))

LIBFILTER_CPP = $(FILTER_SRC)filter_C_interface.cpp
LIBFILTER_HPP = $(FILTER_SRC)filter.h $(FILTER_SRC)filter.hpp $(FILTER_SRC)queue.hpp $(FILTER_SRC)filter_in_use.h
LIBFILTER_CPP_O = $(patsubst $(FILTER_SRC)%.cpp,$(BUILD)%-cpp.o,$(LIBFILTER_CPP))

LIBFILTER_C = $(FILTER_SRC)timer.c
LIBFILTER_H = $(FILTER_SRC)timer.h
LIBFILTER_C_O = $(patsubst $(FILTER_SRC)%.c,$(BUILD)%-c.o,$(LIBFILTER_C))

SERIAL_C = $(SERIAL_SRC)SerialQueue.c $(SERIAL_SRC)socket.c
SERIAL_H = $(SERIAL_SRC)SerialQueue.h $(SERIAL_SRC)socket.h
SERIAL_O = $(patsubst $(SERIAL_SRC)%.c,$(BUILD)%-c.o,$(SERIAL_C))

ALL_FILES = $(MIAN_C) $(MAIN_H) $(LIBCAMERA_CPP) $(LIBCAMERA_HPP) $(LIBCAMERA_C) $(LIBCAMERA_H) $(LIBFILTER_CPP) $(LIBFILTER_C) $(LIBFILTER_H) $(LIBFILTER_HPP) $(SERIAL_C) $(SERIAL_H)

all: install

.depend : $(ALL_FILES)
	gcc -MM $(MAIN_C) $(SERIAL_C) $(LIBCAMERA_C) $(LIBFILTER_C) $(INCLUDES_FLAGS) | sed 's/^[^ ]/$(patsubst %/,%\/,$(BUILD))&/' | awk '{ gsub(/\.o/,"-c.o"); print $0}' > .depend
	g++ -MM $(LIBCAMERA_CPP) $(LIBFILTER_CPP) $(INCLUDES_FLAGS) | sed 's/^[^ ]/$(patsubst %/,%\/,$(BUILD))&/' | awk '{ gsub(/\.o/,"-cpp.o"); print $0}' >> .depend 

-include .depend

$(BUILD):
	@mkdir -p $(BUILD)

do-build: $(BUILD) $(BUILD)motorControl

install : do-build $(LIBS)
	@cp $(BUILD)libcamera.so $(LIBS)
	@cp $(BUILD)libfilter.so $(LIBS)
	@cp $(BUILD)motorControl $(INSTALL)

$(LIBS):$(INSTALL)
	@mkdir -p $(LIBS)

$(INSTALL) :
	@mkdir $(INSTALL) 

%c.o :
	gcc -c $(CFLAGS) $< $(INCLUDES_FLAGS) $(LFLAGS) -o $@

%cpp.o :
	g++ -c $(CPPFLAGS) $< `pkg-config --cflags --libs opencv4` -lyaml-cpp $(INCLUDES_FLAGS) $(LFLAGS) -o $@

$(BUILD)libcamera.so: $(LIBCAMERA_CPP_O)
	g++ $(CPPFLAGS) $(INCLUDES_FLAGS) -shared $^ `pkg-config --cflags --libs opencv4` -lyaml-cpp $(LFLAGS) -o $@

$(BUILD)libfilter.so : $(LIBFILTER_CPP_O)
	g++ $(CPPFLAGS) $(INCLUDES_FLAGS) -shared $^ $(LFLAGS) -o $@

$(BUILD)motorControl : $(BUILD)libcamera.so $(BUILD)libfilter.so $(SERIAL_O) $(MAIN_O) $(SERIAL_O) $(LIBCAMERA_C_O) $(LIBFILTER_C_O)
	gcc $(CFLAGS) $(MAIN_O) $(SERIAL_O) $(LIBCAMERA_C_O) $(LIBFILTER_C_O) $(LFLAGS) -Wl,-rpath=$(LIBS) -L ./ -lcamera -lfilter -o $@

pio:
	pio run -t upload -d arduino/

cu:
	cu -s 115200 -l /dev/ttyACM0

filter:
	python3 filters/filter_calculus.py y

no-autofocus:
	v4l2-ctl --device=/dev/video2 --set-ctrl=focus_automatic_continuous=0

autofocus:
	v4l2-ctl --device=/dev/video2 --set-ctrl=focus_automatic_continuous=1

.PHONY: all clean cu filter pio install do-build

clean:
	@rm -rf $(BUILD) *~ .depend $(INSTALL)