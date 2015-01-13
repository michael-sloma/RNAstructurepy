CC = gcc
CFLAGS = -DSHORT -O3 -Wall -D disablecoax
CUDA_CC = nvcc
CUDA_CFLAGS = -DSHORT -O3 -use_fast_math -D disablecoax 
DAT = ../data_tables
BIN = ../exe
