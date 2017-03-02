tests:
	nosetests
dev_test:
	fswatch -o -0 -e '__.*' NCBITK/tests | xargs -0 -n 1 -I {} make tests
	# fswatch -o NCBITK/tests | (while read; do make tests; done)
	# fswatch -o -e '__.*' NCBITK/tests | xargs -n 1 make tests
	# fswatch -o NCBITK/tests | xargs -n 1 make tests
	# fswatch -0 NCBITK/tests | xargs -0 -n 1 make tests
	# fswatch -o NCBITK/tests | xargs -n 1 make tests
	# fswatch -0 NCBITK/tests | xargs -0 -n 1 make tests
	# fswatch -0 NCBITK/tests | xargs -0 -n 1 make tests
	# fswatch -0 NCBITK/tests | xargs -0 -n 1 make tests
