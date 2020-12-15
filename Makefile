ROOTlibs=$(shell root-config --libs)
ROOTflags=$(shell root-config --cflags)

INCs=inc

INCfiles=inc/chebFitter.h inc/nodes.h

main: obj/chebFitter.o  obj/nodes.o  obj/main.o 
	g++ -O3 -I${INCs} $^   -o $@ ${ROOTlibs}  -fopenmp



obj/%.o: src/%.cc  ${INCfiles}
	g++ -O3 -c -I${INCs} ${ROOTflags}  $<   -o $@ -fopenmp


clean:
	rm main obj/*.o
