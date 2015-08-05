CFLAGS = `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`

% : %.cpp
	g++ -I/usr/local/OpenCV/include -L/usr/local/OpenCV/lib $(CFLAGS) $(LIBS) -lopencv_nonfree -o $@ $<
	

