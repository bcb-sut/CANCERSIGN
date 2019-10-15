 #!/bin/bash

cd "${0%/*}"

if [ ! -d temp ]; then
	mkdir temp
fi

cd temp
if [ -f pid.txt ]; then
	p_id=`cat pid.txt`
	kill -9 $p_id  
fi
rm -r *
cd ..

cd src
Rscript -e 'library(methods); shiny::runApp(launch.browser=TRUE)'
