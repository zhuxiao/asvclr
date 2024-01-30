#!/bin/sh
	if [ ! -d "bin" ]; then 
		mkdir bin
	fi
	cd src/
	make
	mv -f asvclr ../bin
	cd ../

