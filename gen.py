#!/bin/python
import os
import sys
import math
import argparse
import subprocess
import time
import multiprocessing
from multiprocessing import Process, Queue
from subprocess import Popen, PIPE, STDOUT
#import scipy.special


# cmake target to build and execute
CMAKE_TARGET            = "code"
CMAKE_TARGET_INCLUDE    = "main.h"
CMAKE_TARGET_DIR        = ""
CMAKE_TARGET_FLAG       = ""
CMAKE_LOGGING_FILE      = "out.log"

def wait_timeout(proc, seconds):
	"""Wait for a process to finish, or raise exception after timeout"""
	start = time.time()
	end = start + seconds
	interval = min(seconds / 1000.0, .25)
	while True:
		result = proc.poll()
		if result is not None:
			return result
		if time.time() >= end:
			#os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
			proc.kill()
			raise RuntimeError("Process timed out")

		time.sleep(interval)


def rebuild():
	p = Popen(["make", CMAKE_TARGET, "-j1"], stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True, cwd= "./cmake-build-release")
	p.wait()
	if p.returncode != 0:
		print("ERROR Build", p.returncode, p.stdout.read().decode("utf-8"))
	return p.returncode


def run(seconds=-1):
	start = time.time()
	f = open(CMAKE_LOGGING_FILE, "w")

	#stdout=f, stderr=f,
	p = Popen(["./" + CMAKE_TARGET_DIR+CMAKE_TARGET, CMAKE_TARGET_FLAG], preexec_fn=os.setsid, cwd= "./cmake-build-release")
	if seconds != -1:
		try:
			c = wait_timeout(p, seconds)
			t = time.time() - start
			# print("runtime: ", t)
			return (c, t)
		except:
			# print("runtime is over")
			return -1, -1
	else:
		p.wait()
		t = time.time() - start
		# print("runtime: ", t)
		return p.returncode, t

def write_config(r:int, m: int, s: int, t: int, threads: int, ITERS=20):
	""""""

	with open(CMAKE_TARGET_INCLUDE, "w") as f:
		f.write("""
#ifndef CONFIG_SET
#define CONFIG_SET

#include <iostream>

#define BENCHMODE
""")

		f.write("constexpr uint32_t n=" + str(1) + ";\n")
		f.write("constexpr uint32_t r=" + str(r) + ";\n")
		f.write("constexpr uint32_t m=" + str(m) + ";\n")
		f.write("constexpr uint32_t threads=" + str(threads) + ";\n")
		f.write("constexpr uint32_t s= uint32_t(1) << " + str(s) + ";\n")
		f.write("constexpr uint32_t t= uint32_t(1) << " + str(t) + ";\n")
		f.write("constexpr uint32_t BENCHES=" + str(ITERS) + ";\n")
		f.write("constexpr uint32_t BENCHESMULT=" + str(ITERS*10000) + ";\n")
		f.write("#endif //CONFIG_SET")


def bench(args):	
	global CMAKE_LOGGING_FILE
	for r in range(15, args.r):
		for m in [2**(r//3)]: #, 2**r//3
			s = math.ceil(r/3)
			t = math.ceil(r/3)
			CMAKE_LOGGING_FILE  = "r" + str(r) + "_m" + str(m) + "_t" + str(t) + "_s" + str(s) + "_iters" + str(args.iterations) + ".log"
			print("\nCurrent Set: ", CMAKE_LOGGING_FILE)
			write_config(r, m, s, t, args.threads, args.iterations)

			if rebuild() != 0:
				print("ERROR Build")
				continue

			run(args.seconds)


def bench_inner(args):
	global CMAKE_TARGET
	global CMAKE_LOGGING_FILE
	threads=1

	for r in range(15, args.r):
		CMAKE_TARGET="code_inner"
		CMAKE_LOGGING_FILE  = "inner_r" + str(r) + "_m0_t0_s0_iters" + str(args.iterations) + ".log"
		write_config(r, r, r, r, args.threads, args.iterations)

		if rebuild() != 0:
			print("ERROR Build")
			continue

		run(args.seconds)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Skip some pages.')
	parser.add_argument('-r', help='log scale', default=24, type=int, required=False)
	parser.add_argument('-m', help='.', default=1, type=int, required=False)
	parser.add_argument('-s', help='log scale', default=8, type=int, required=False)
	parser.add_argument('-t', help='log scale', default=8, type=int, required=False)
	parser.add_argument('--threads', help='', default=1, type=int, required=False)
	parser.add_argument('--iterations', help='', default=5, type=int, required=False)
	parser.add_argument('--bench', action='store_true')
	parser.add_argument('--inner', action='store_true')
	parser.add_argument('--seconds', help='Kill a programm after x seconds. Default -1: run infinite', default=-1, type=int, required=False)

	parser.add_argument('--executable', help='', default="code", type=str, required=False)
	parser.add_argument('--include', help='', default="main.h", type=str, required=False)

	args = parser.parse_args()

	CMAKE_TARGET = args.executable
	CMAKE_TARGET_INCLUDE = args.include

	if args.bench:
		if args.inner:
			bench_inner(args)
		else:
			bench(args)
		exit()

	write_config(args.r, args.m, args.s, args.t, args.threads, args.iterations)

	if rebuild() != 0:
		print("ERROR Build")
		exit()

	print("Build Finished")
	retcode, time = run(args.seconds)
	if retcode == -1:
		print("ERROR Runtime over after", args.seconds, "seconds")
