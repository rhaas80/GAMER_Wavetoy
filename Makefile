.PHONY: clean

all: wavetoy
	./wavetoy >out.asc
	bash plot_data.sh

wavetoy: wavetoy.cc main.cc
	bash main.cc

clean:
	rm -f out.asc wavetoy frame*.png *~
