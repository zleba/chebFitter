ROOTlibs=$(shell root-config --libs)
ROOTflags=$(shell root-config --cflags)

INCs=-Iinc -I/usr/include/eigen3

INCfiles=inc/chebFitter.h  inc/chebFitter2D.h inc/nodes.h

tagV:  obj/tagV.o obj/chebFitter.a  
	g++ -O3   $^   -o $@     ${ROOTlibs}  -fopenmp


testFit:  obj/testFit.o obj/chebFitter.a  
	g++ -O3   $^   -o $@     ${ROOTlibs}  -fopenmp

obj/chebFitter.a: obj/chebFitter.o  obj/chebFitter2D.o  obj/nodes.o
	ar rvs $@   $^


obj/%.o: src/%.cc  ${INCfiles}
	g++ -O3 -c ${INCs} ${ROOTflags}  $<   -o $@ -fopenmp


clean:
	rm testFit obj/*.o
