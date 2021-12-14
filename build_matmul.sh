#!/bin/bash

clang++ -o matmul -std=c++11 -O3 -fsanitize=address matmul.cpp

time ./matmul
