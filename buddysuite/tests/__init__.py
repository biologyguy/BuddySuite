# coding=utf-8
import sys
import os
main_dir = os.path.dirname(os.path.abspath(__file__)).split("/")[:-1]
main_dir = "/".join(main_dir)
sys.path.insert(0, main_dir)