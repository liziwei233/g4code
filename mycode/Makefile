OBJDIR = Object
SRCDIR = Source
INCDIR = Include

HEADERS = $(wildcard $(INCDIR)/*.h)
OBJECTS = $(subst $(INCDIR), $(OBJDIR), $(HEADERS:.h=.o))

CC = g++
CFLAGS = -O2 `root-config --libs --cflags` \
		 -I$(INCDIR)

LDFLAGS = `root-config --glibs` 

all: $(OBJECTS) draw
#testbeam: $(OBJECTS) 

#draw: $(OBJECTS) draw.cpp 
#	@echo "Building executable 'draw'..."
#	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


$(OBJDIR)/DrawMyfunc.o: $(SRCDIR)/DrawMyfunc.cpp $(INCDIR)/DrawMyfunc.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm $(OBJDIR)/*.o
