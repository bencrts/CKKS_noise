# This sample should be run after make libHEAAN.a file in lib folder

# All Target
encryption: clean HEAAN-encryption run
complex: clean HEAAN-complex run
nr: clean HEAAN-nr run

# Tool invocations

HEAAN-encryption:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++   encryption/main.cpp  -o HEAAN -L/opt/homebrew/lib ../lib/libHEAAN.a -I/../src -lntl -lgmp -lm -O2 -std=c++11 -pthread
	@echo 'Finished building target: $@'
	@echo ' '

HEAAN-nr:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ NR/main.cpp -o HEAAN -L/opt/homebrew/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '

HEAAN-complex:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ complex/main.cpp -o HEAAN -L/opt/homebrew/lib ../lib/libHEAAN.a  -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	rm -rf HEAAN

run:
	./HEAAN
