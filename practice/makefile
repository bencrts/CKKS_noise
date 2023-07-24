textbook:
	rm -rf textbook
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ textbook-ckks.cpp -o textbook -L/opt/homebrew/lib ../HEAAN/HEAAN/lib/libHEAAN.a  -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '
	./textbook

rns:
	rm -rf rns
	g++ -std=c++11 -pthread -o "rns" rns-ckks.cpp ../FullRNS-HEAAN/lib/libFRNSHEAAN.a -I/../src/
	@echo 'Finished building target: $@'
	@echo ' '
	./rns

clean:
	rm -rf textbook
	rm -rf rns
