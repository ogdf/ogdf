OGDF_INCLUDE = ../../include/
CC = g++
CXX = g++
CPPFLAGS = -I$(OGDF_INCLUDE) -std=c++11
ifdef DEBUG
  OGDF_LIBS = ../../_debug/
  CPPFLAGS += -W -Wall -DOGDF_DEBUG -Wno-unused-parameter -O0 -g -ggdb
else
  OGDF_LIBS = ../../_release/
  CPPFLAGS += -O3 -DNDEBUG
endif
LDFLAGS = -pthread -L$(OGDF_LIBS) -L../..
LDLIBS = -lOGDF -lCOIN

all: $(OUTPUT) $(OUTPUTS)

$(OUTPUT): $(OUTPUT:%=%.o) $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%o: %cpp

clean:
	$(RM) $(OUTPUTS) $(OUTPUTS:%=%.o) $(OUTPUT) $(OUTPUT:%=%.o) $(OBJ)

.PHONY: clean
