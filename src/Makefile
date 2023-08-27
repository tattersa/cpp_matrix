CC = g++
CPPFLAGS = -lgtest -std=c++17
CHECKFLAGS = -Wall -Werror -Wextra
SRC = s21_matrix_oop.cc
TEST = tests.cc

all: test gcov_report

test: clean
	@${CC} ${TEST} ${SRC} ${CHECKFLAGS} ${CPPFLAGS} -o test
	@./test

%.o: %.cc
	@${CC} -std=c++17 ${CHECKFLAGS} -c $^

%.a: %.o
	@ar rcs $@ $^
	@ranlib $@

clean:
	@/bin/rm -rf *.o *.a test *.gcno *gcda report *.info test *.out

gcov_report: ${TEST} ${SRC}
	@$(CC) --coverage $^ ${CPPFLAGS} ${CHECKFLAGS} -o test
	@./test
	@lcov -t "test" -o test.info --no-external -c -d ./
	@genhtml -o report test.info
	@open report/index.html

style:
	clang-format --style=google -i *.cc *.h

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem . *.cc


valgrind:
	valgrind --tool=memcheck --leak-check=full ./test
