tests:
	nosetests
dev_test:
	fswatch -0 NCBITK/tests | xargs -0 -n 1 make tests
