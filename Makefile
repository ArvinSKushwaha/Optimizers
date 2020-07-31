main: main.cpp
	${CXX} $< -o $@ -O3

clean:
	rm ./main

run:
	./main