FLAGS = -std=c11 -Wall -Werror -Wextra

all: test s21_matrix.a gcov_report

s21_matrix.a: *.c
	gcc $(FLAGS) -c *.c
	ar rc s21_matrix.a *.o
	ranlib s21_matrix.a
	rm -rf *.o

test: tests/*.check s21_matrix.a
	checkmk clean_mode=1 tests/*.check > test.c
	gcc $(FLAGS) test.c s21_matrix.a -lcheck -lm -lpthread -o test
	./test

gcov_report: test
	gcc $(FLAGS) -coverage *.c -o gcovreport -lcheck -lm
	./gcovreport
	lcov -t "gcovreport" -o gcovreport.info -c -d .
	genhtml -o report gcovreport.info
	open report/./index.html

clean: 
	rm -rf test.c *.o *.a test *.gcda *.gcno gcovreport* report