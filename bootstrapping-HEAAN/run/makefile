# This sample should be run after make libHEAAN.a file in lib folder

# All Target
all: clean BOOTS

# Tool invocations
BOOTS:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ boots.cpp -o BOOTS -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	rm -rf BOOTS
