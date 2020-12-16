ROOTlibs=$(shell root-config --libs)
ROOTflags=$(shell root-config --cflags)

INCs=-Iinc -I/usr/include/eigen3

INCfiles=inc/chebFitter.h inc/nodes.h

main: obj/chebFitter.o  obj/nodes.o  obj/main.o 
	g++ -O3  ${INCs} $^   -o $@ ${ROOTlibs}  -fopenmp



obj/%.o: src/%.cc  ${INCfiles}
	g++ -O3 -c ${INCs} ${ROOTflags}  $<   -o $@ -fopenmp


clean:
	rm main obj/*.o
