VERSION=1.1
default:
	- rm -rf build dist
	- export VERSION=$(VERSION) ;\
	  ./setup.py --quiet bdist_rpm  
	- sudo rpm -e resultsFile
	- sudo rpm -hiv dist/resultsFile*.noarch.rpm
