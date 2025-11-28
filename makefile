.env:
	echo "ROOT_DIR=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))" | tee .env
