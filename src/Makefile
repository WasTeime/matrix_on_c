CC := gcc
CFLAGS := -std=c11 -Wall -pedantic -Wextra -Wfloat-equal -Wfloat-conversion -Wvla
UNITFLAGS := -lcheck
COV_FLAG :=
OS:=$(shell uname -s)
ifeq ($(OS), Darwin)
	UNITFLAGS += $(shell pkg-config --cflags --libs check)
else
	UNITFLAGS += -lsubunit -lrt -lm -lpthread
endif

# Объектные файлы
OBJS := ./out/s21_matrix.o
UNITOBJS := ./out/s21_matrix_test.o ./out/s21_allocate_free_test.o ./out/s21_eq_test.o ./out/s21_sum_sub_test.o ./out/s21_mul_test.o ./out/s21_transpose_test.o ./out/s21_det_test.o ./out/s21_complements.o ./out/s21_minor_test.o ./out/s21_inverse_test.o

.PHONY : all clean gcov_report del_lib

test : unit_tests.exe
	./unit_tests.exe || echo

all : COV_FLAG += --coverage
all : clean gcov_report

unit_tests.exe : $(OBJS) $(UNITOBJS) s21_matrix.a
	$(CC) -o unit_tests.exe $(UNITOBJS) $(UNITFLAGS) -L. -ls21_matrix -fprofile-arcs -ftest-coverage

gcov_report : COV_FLAG += --coverage
gcov_report : del_lib test
	lcov -t "coverage" -o coverage.info -c -d . 
	genhtml -o report coverage.info
	rm -rf ./*.g*

del_lib :
	rm -rf ./*.a
	rm -rf ./s21_matrix.o
	rm -rf ./out/s21_matrix.o
	
s21_matrix.a : $(OBJS)
	ar rc $@ $^
	ranlib $@
	cp $@ lib$@

./out/%.o : ./unit_tests/%.c
	mkdir -p ./out
	$(CC) $(CFLAGS) $(UNITFLAGS) -c $<
	mv ./*.o ./out/

./out/%.o : ./%.c ./s21_matrix.h
	mkdir -p ./out
	$(CC) $(CFLAGS) $(COV_FLAG) -c $<
	mv ./*.o ./out/

clean :
	rm -f ./*.exe
	rm -rf ./out
	rm -rf ./*.a
	rm -rf ./*.g*
	rm -rf ./report
	rm -rf ./*.info
	
memcheck: unit_tests.exe
	valgrind --leak-check=summary --show-leak-kinds=definite --track-origins=yes ./unit_tests.exe || echo "Memory check completed"
	@if [ $$? -eq 0 ]; then \
		echo "No memory leaks detected."; \
	else \
		echo "Memory leaks found."; \
	fi