spaghetti.r
=========

Скрипт для построения диагностических графиков в [модели уточнения Бушмаринова](https://github.com/dmitrienka/metarefine).


````
$ spaghetti.r  -h
Usage: spaghetti.r [options] directory


Options:
	-W WIDTH, --width=WIDTH
		Width of the output file, [default 12]

	-H HEIGHT, --height=HEIGHT
		Height of the output file, [default 7]

	-p, --Plot
		Use existing _rtable.dat for plotting, [default FALSE]

	-e, --Errors
		Use bond errors in outlier detection, [default FALSE]

	-k KPAR, --Kpar=KPAR
		Use k in the outlier detection function, negatives for table values, [default -1]

	-a, --Angles
		Use angles instead of bonds

	-h, --help
		Show this help message and exit

````
