output: TestStraws.o
	g++ TestStraws.o -o TestStraws `root-config --libs` `root-config --cflags`

TestStraws.o: TestStraws.C
	g++ -g -c TestStraws.C `root-config --libs` `root-config --cflags`
	
clean:
	rm -f *.o
