LIBS := -lglfw3 -lpthread -lX11 -lXrandr -lXxf86vm -lXi -lGL 

all: muFLIP

muFLIP: muflip.cpp
	g++ -Wall muflip.cpp $(LIBS) -o muFLIP

clean:
	rm muFLIP
